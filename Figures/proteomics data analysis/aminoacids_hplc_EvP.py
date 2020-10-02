__author__ = "Sieze Douwenga"

# If used in ipython, WD should be usr\surfdrive\Shared\!Labjournal\Chemostats\05Chemostats_WT_ccpa

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
DILUTIONRATE = 0.5 # h-1

pd.set_option('display.max_columns', 100)
# Import metadata ###
metadata = pd.read_csv("data/raw_fluxdata/metadata.txt",
                       delimiter = '\t')
# Import amino acid hplc data ###
hplcdata = pd.read_csv("data/raw_fluxdata/hplc_aa/20190314_05chemostat_aa_data_complete.txt",
                      delimiter = '\t')
# Import manually curated hplcdata
manualentry = pd.read_csv("data/raw_fluxdata/hplc_aa/manual_concentrations_run1.txt",
                          delimiter = '\t')
# Import amino acid concentrations in medium kby design
mediumdesign = pd.read_csv("data/raw_fluxdata/medium_bydesign/aminoacids.txt",
                          delimiter = '\t')

# Overwrite incorrect entries in automated hplcdataset with manually curated data
for sample,compound,concentration,area in\
        zip(manualentry.Sample, manualentry.Compound, manualentry.Concentrations, manualentry.Area):
    sampleindex = hplcdata[hplcdata['Sample name'] == sample].index
    compoundindeces = hplcdata.filter(like=compound).columns # output: [concentration, area]
    hplcdata.at[sampleindex,compoundindeces[0]] = concentration
    hplcdata.at[sampleindex,compoundindeces[1]] = area

hplcdata.rename(columns = {'Sample name':'hplcid'}, inplace=True)
# extract hplcrun information from hplcid and as separate column
runs = []
replicates = []
hplcreps = []
reactornum = []
strain = []
for name in hplcdata['hplcid']:
    try:
        runs.append(re.findall(r'run[0-9]',name)[0])
    except(IndexError):
        runs.append('run1')
    try:
        replicates.append(re.findall(r't[0-9]',name)[0][-1])
    except(IndexError):
        replicates.append(1)
    try:
        hplcreps.append(re.findall(r'rep[0-9]of[0-9]',name)[0])
    except(IndexError):
        hplcreps.append('rep1of1')
    try:
        reactornum.append(re.findall(r'R[0-9]',name)[0])
    except(IndexError):
        try:
            reactornum.append(re.findall(r'b[0-9]',name)[0])
        except(IndexError):
            reactornum.append(0)
    try:
        strain.append(re.findall(r'WT',name)[0])
    except(IndexError):
        try:
            strain.append(re.findall(r'ccpa',name)[0])
        except(IndexError):
            strain.append('nostrain')


hplcdata.insert(1,'run',np.array(runs))
hplcdata.insert(2,'vesselreplicate',np.array(replicates))
hplcdata.insert(3,'hplcrunreps',np.array(hplcreps))
hplcdata.insert(2,'reactornum',np.array(reactornum))
hplcdata.insert(2,'strain',np.array(strain))

# Remove chemostat data from 1st run if it was run again in the second run (some overlap is present)
chemid = metadata[metadata['vesseltype']=='chemostat']['id']
overlappedhplcid = hplcdata[[(idval in list(chemid) and runval== 'run2') for idval,runval in zip(hplcdata.id,hplcdata.run)]]['id']
hplcdata1 = hplcdata[-hplcdata['id'].isin(overlappedhplcid)]
hplcdata2 = hplcdata[hplcdata['run']=='run2']
hplcdata2 = hplcdata2[hplcdata2['id'].isin(overlappedhplcid)]
hplcdata = pd.concat([hplcdata1,hplcdata2])

# remove non concentration data
hplcdata = hplcdata.drop([colname for colname in hplcdata.columns if 'Area' in colname],axis=1)

# remove outliers: calibration 100% run 1; and ccpa R3 t1
calib100id = metadata[(metadata['vesseltype']=='calibration') & (metadata['strength'] == '100%')]['id'].values[0]
hplcdata = hplcdata[-((hplcdata['id']==calib100id) & (hplcdata['run']=='run1'))]
hplcdata = hplcdata[-((hplcdata['strain']=='ccpa') & (hplcdata['reactornum'] == 'R3') &\
                     (hplcdata['vesselreplicate']=='1'))]

hplcdata.loc[:,hplcdata.filter(like='Conc').columns] = hplcdata.loc[:,hplcdata.filter(like='Conc').columns]/10

# Import dry weight data, annotate correctly using metadata id's and normalise ###
dwdata_ccpa= pd.read_csv("data/raw_fluxdata/dryweights/20181128_dry_weights_ccpa.txt",
                     delimiter = '\t')
dwdata_ccpa['strain'] = 'ccpa'
dwdata_wt= pd.read_csv("data/raw_fluxdata/dryweights/20181207_dry_weights_WT.txt",
                         delimiter = '\t')
dwdata_wt['strain'] = 'wt'
dwdata = pd.concat([dwdata_ccpa,dwdata_wt],ignore_index=True)
dwdata['reactor'] = [float(re.findall(r'[0-9]+',stringentry)[0]) for stringentry in dwdata['reactor']]
dwdata = dwdata.melt(id_vars=['strain','reactor','measurement'],value_name="value", var_name='replicate')
dwdata = dwdata.pivot_table(index=["strain","reactor","replicate"],columns="measurement",values="value")
dwdata.reset_index(inplace=True)
dwdata.rename(columns={'reactor':'vesselnum'},inplace=True)
dwdata.columns.name = None

# add correct metadata ids to dryweight data
dwids = []
for id,strain,vesselnum,vesseltype in zip(metadata.id,metadata.strain,metadata.vesselnum,metadata.vesseltype):
    if vesseltype!= 'chemostat':
        continue
    for dwstrain,dwvesselnum in zip(dwdata.strain,dwdata.vesselnum):
        if dwstrain == strain and dwvesselnum == vesselnum:
            dwids.append(id)
dwdata.insert(0,'id',np.array(dwids))
dwdata.drop(['strain','vesselnum'],axis=1,inplace=True)

# calculate the dryweight with sdev #####
dwdata['dw_gpL'] = (dwdata['post_mg'] - dwdata['pre_mg']) / dwdata['volume_mL']
meandw = pd.Series(dwdata.groupby('id')['dw_gpL'].mean(), name='mean')
sdevdw = pd.Series(dwdata.groupby('id')['dw_gpL'].std(), name='sdev')
meandw = pd.concat([meandw,sdevdw],axis=1)

# calculate mean aa hplc values for technical reps
hplcmeanreps = hplcdata.groupby(['id','strain','reactornum','vesselreplicate','run']).mean()
hplcmstd = hplcdata.groupby(['id','strain','reactornum','vesselreplicate','run']).sem()

# split data into medium and chemostat aas
chemoids = list(metadata[metadata['vesseltype']=='chemostat']['id'])
chemohplc = hplcmeanreps.loc[hplcmeanreps.index.get_level_values('id').isin(chemoids),:]
mediumid = list(metadata[metadata['vesseltype']=='medium']['id'])
mediumhplc = hplcmeanreps.loc[hplcmeanreps.index.get_level_values('id').isin(mediumid),:]

# relative difference between design and medium measurement in the hplc
flipped = mediumhplc.T.sort_index()
difference = flipped.subtract(mediumdesign['Concentration'].values, axis=0)
reldifference = difference.divide((mediumdesign['Concentration'].values),axis=0)*100

# plot raw hplc data ####
#mediumhplc.plot(marker='o', colormap='tab20',linestyle='')
#plotchemo = chemohplc.groupby(['strain','reactornum']).plot(marker='o',colormap='tab20')
# plot the calibration curves (mean values if there were replicates)
calibids = list(metadata[metadata['vesseltype']=='calibration']['id'])
calibhplc = hplcmeanreps.loc[hplcmeanreps.index.get_level_values('id').isin(calibids),:]
#calibhplc.groupby('run').plot(marker='o', colormap='tab20')

# flux calculations ###
# aa concentrations
meanmedium = mediumhplc.groupby(['id','strain','reactornum']).mean()

# select the medium vessels that were connected to the chemostat during sampling (multiple vessels were used in some cases
# WT only had 1 bottle, so the 2nd bottle measurement in the second run was likely mislabeled
b1wtid = metadata.loc[(metadata['vesseltype'] == 'medium') & (metadata['strain'] == 'wt') & (metadata['vesselnum'] == 1),'id'].values[0]
b3ccpaid = metadata.loc[(metadata['vesseltype'] == 'medium') & (metadata['strain'] == 'ccpa') & (metadata['vesselnum'] == 3),'id'].values[0]
samplemediumids = [b1wtid,b3ccpaid]
mediumduringsampling = meanmedium.loc[meanmedium.index.get_level_values('id').isin(samplemediumids),:]

grpchemo = chemohplc.groupby(['id','strain','reactornum'])
meanchemo = grpchemo.mean()
semchemo = chemohplc.groupby(['id','strain','reactornum']).sem()
#meanchemo.to_csv("data/CcpA_chemostat/AA_residual_conc.txt",sep='\t',index=None)
#semchemo.to_csv("data/CcpA_chemostat/AA_residual_conc_sem.txt",sep='\t',index=None)
# dw concentrations
meandw = dwdata.groupby(['id'])['dw_gpL'].mean()
semdw = dwdata.groupby(['id'])['dw_gpL'].sem()

#flux
sorteddf = meanchemo.sort_index(axis=1)
concentrationdiff = sorteddf.subtract(mediumdesign['Concentration'].values, axis=1)
flux = concentrationdiff.divide(meandw.values, axis=0) * DILUTIONRATE # mmol/gDW h-1
#flux.to_csv("data/CcpA_chemostat/Fluxdata.txt",sep='\t',index=None)

#fluxerror
relerrchemosquared = (semchemo/meanchemo)**2
relerrchemosquared = relerrchemosquared.fillna(0)
relfluxerror = np.sqrt(relerrchemosquared.add((semdw/meandw).values**2,axis=0)) # as fraction
fluxerror = relfluxerror * abs(flux)
#fluxerror.to_csv("data/CcpA_chemostat/Fluxerror.txt",sep='\t',index=None)