import pandas as pd
import matplotlib as mpl
mpl.use('Qt4Agg')  # an upgrade to the regular interactive plots
import matplotlib.pyplot as plt
import seaborn as sns

#TODO: compare aa concentrations to medium

MU = 0.5 # h-1
#biomass_MW = 24.82 #g/mol   b.subtilis, source: https://bionumbers.hms.harvard.edu/bionumber.aspx?id=111949&ver=3&trm=elemental+biomass+formula&org=
biomass_MW = 27.4 #g/mol uit het hoofd van Rinke/Eunice lactis IL...

# Exportfilename
exportname = 'concentrations_chemostat05_ccpaWT.txt'

########### Import data ##################
# hplcdata
hplcfile = "data/raw_fluxdata/hplc_oa/20191030_internalstandard_chemostats05_ccpa_WT_MG1363_clean.tsv"
hplcdata_raw = pd.read_csv(hplcfile, sep='\t')
hplcdata_raw = hplcdata_raw.set_index('Sample name', drop=True)
hplcdata_raw.index.name = 'samplename'

# chemical formulae
chemformfile = "data/raw_fluxdata/hplc_oa/hplc_compoundformulae.txt"
chemform = pd.read_csv(chemformfile, sep='\t')

# dilution rates
dilutionratesfile = "data/raw_fluxdata/dilution_rates/dilutionrates.txt"
dilutionrates = pd.read_csv(dilutionratesfile, sep='\t')
dilutionrates.set_index(['strain','reactor'],inplace=True)

# dryweights
dryweightfile_ccpa = "data/raw_fluxdata/dryweights/20181128_dry_weights_ccpa.txt"
dryweightfile_wt   = "data/raw_fluxdata/dryweights/20181207_dry_weights_WT.txt"
dw_ccpa = pd.read_csv(dryweightfile_ccpa, index_col=0, sep='\t')
dw_wt = pd.read_csv(dryweightfile_wt, index_col=0, sep='\t')

# convert dryweights to long format and merge wt and ccpa
dw_ccpa.reset_index(inplace=True)
dw_ccpa = dw_ccpa.melt(id_vars=['reactor','measurement'],var_name='timepoint',value_name='weight_mg')
dw_ccpa = dw_ccpa.pivot_table(index=['reactor','timepoint'],columns='measurement',values='weight_mg').reset_index()
dw_ccpa.columns.name=None
dw_ccpa.insert(0,'strain',['ccpa']*len(dw_ccpa))
dw_wt.reset_index(inplace=True)
dw_wt = dw_wt.melt(id_vars=['reactor','measurement'],var_name='timepoint',value_name='weight_mg')
dw_wt = dw_wt.pivot_table(index=['reactor','timepoint'],columns='measurement',values='weight_mg').reset_index()
dw_wt.columns.name=None
dw_wt.insert(0,'strain',['WT']*len(dw_wt))
dw_raw = pd.merge(dw_ccpa,dw_wt,how='outer')

dw_raw['dryweight_gpl'] = (dw_raw['post_mg'] - dw_raw['pre_mg'])/dw_raw['volume_mL']
dw = dw_raw.groupby(['strain','reactor'])['dryweight_gpl'].describe()

## Take only concetrations columns
concentration_locs = [colname for colname in hplcdata_raw.columns if "_Conc" in colname]
concentrations = hplcdata_raw.loc[:,concentration_locs]

indexnames = concentrations.index.name
concentrations.reset_index(inplace=True)
concentrations = concentrations.melt(id_vars=indexnames,var_name='Compound',value_name='concentration')
concentrations['Compound'].\
    replace(regex=True,inplace=True,to_replace='_Conc',value='')

# drop all columns containing MQ
concentrations = concentrations[~concentrations.samplename.str.contains("MQ")]


## Extract metadata in filename to data in separate columns
reactorcol = []
timepointcol = []
straincol = []
hplcrepcol = []
ismediumcol=[]
for name in concentrations.samplename:
    splitname = name.split('_')
    if 'R' in name:
        reactorcol.append(splitname[-2])
        timepointcol.append(splitname[-1])
        straincol.append(splitname[1])
        hplcrepcol.append('m1')
        ismediumcol.append(False)
    else:
        if 'CDMpc' in name:
            reactorcol.append('medium_'+splitname[-2])
            hplcrepcol.append(splitname[-1])
            ismediumcol.append(True)
        else:
            reactorcol.append('n/a')
            hplcrepcol.append('n/a')
            ismediumcol.append(False)
        timepointcol.append('t0')
        straincol.append(splitname[1])

concentrations.insert(1,'hplcrep',hplcrepcol)
concentrations.insert(1,'timepoint',timepointcol)
concentrations.insert(1,'reactor',reactorcol)
concentrations.insert(1,'strain',straincol)
concentrations['timepoint'] = timepointcol

# For OAA HPLC only glucose is present in medium - i.e. remove other compounds
concentrations = concentrations.loc[-((concentrations.reactor.isin(['medium_bottle1','medium_bottle2','medium_chembottle','medium_onbottle'])) &\
                   (concentrations.loc[:,'Compound'] != 'Glucose')),:]

for compound,cnum in zip(chemform.Compound, chemform.c):
    concentrations.loc[concentrations.Compound==compound, 'mM_C'] = concentrations.loc[concentrations.Compound==compound, 'concentration'] * cnum

# Add the mean dry weight as a mM_C for final carbon balance dataframe
dw2 = dw.drop('count',axis=1).iloc[:,0:1].reset_index()
dw2['Biomass'] = dw2['mean']/ biomass_MW * 1000
dw2.drop('mean',inplace=True,axis=1)

carbonbalance = concentrations.pivot_table(values='mM_C', index=['samplename','strain','reactor'],columns='Compound').reset_index()
carbonbalance = pd.merge(carbonbalance,dw2,how='outer')
carbonbalance.loc[carbonbalance.strain=='WT','C_in'] =\
    carbonbalance.loc[carbonbalance.samplename=='20190110_WT_CDMpc_chembottle_m1','Glucose'].values[0]
carbonbalance.loc[carbonbalance.strain=='ccpa','C_in'] = \
    carbonbalance.loc[carbonbalance.samplename=='20190110_ccpa_CDMpc_bottle2_m1','Glucose'].values[0]

carbonbalance.set_index(['samplename','strain','reactor'],inplace=True)

carbonbalance_nodw = carbonbalance.drop('Biomass',axis=1)
carbonbalance_nodw['rel_C_balance'] = (carbonbalance_nodw.C_in - carbonbalance_nodw.iloc[:,0:carbonbalance_nodw.shape[1]-1].sum(axis=1))/carbonbalance_nodw.C_in
carbonbalance['rel_C_balance'] = (carbonbalance.C_in - carbonbalance.iloc[:,0:carbonbalance.shape[1]-1].sum(axis=1))/carbonbalance.C_in

# Determine fluxes
fluxdata = concentrations.pivot_table(values='concentration', index=['samplename','strain','reactor'],\
                                      columns='Compound').reset_index()
fluxdata.loc[fluxdata.strain=='WT','Glucose'] = \
    fluxdata.loc[fluxdata.samplename=='20190110_WT_CDMpc_chembottle_m1','Glucose'].values[0]*-1
fluxdata.loc[fluxdata.strain=='ccpa','Glucose'] = \
    fluxdata.loc[fluxdata.samplename=='20190110_ccpa_CDMpc_bottle2_m1','Glucose'].values[0]*-1
fluxdata = fluxdata.groupby(['strain','reactor']).describe()
fluxdata = fluxdata.loc[:,fluxdata.columns.get_level_values(None).isin(['mean'])]
fluxdata.columns = fluxdata.columns.droplevel(None)
fluxdata = fluxdata.loc[fluxdata.index.get_level_values('reactor').isin(['R1','R2','R3','R4']),:]
dw = dw.loc[dw.index.get_level_values('reactor').isin(['R1','R2','R3','R4']),:]
dw = dw.loc[:,['mean','std','count']]
dw['rel_sd_squared'] = (dw['std']/dw['mean'])**2
fluxdata = fluxdata.divide(dw['mean'],axis=0)
fluxdata = fluxdata.multiply(dilutionrates['dilutionrate'],axis=0)
fluxdata['unit'] = 'mmol/gDW/h'
fluxdata.reset_index(inplace=True)
#fluxdata.to_csv("data/CcpA_chemostat/2018_05chemostats_wtccpa_fluxdata_organic_acids.txt",sep='\t',index=None)
fluxdata = fluxdata.melt(id_vars=['strain','reactor','unit'],value_name='flux')


## Plot concentrations
sns.set_style('darkgrid',{'axes.grid':True})
sns.set_context(context='talk')
g = sns.FacetGrid(data=fluxdata,col='strain')
g.map(sns.stripplot,"Compound","flux",alpha=0.5,zorder=1,jitter=True,size=10,color='orange')
g.map(sns.pointplot,"Compound","flux",linestyles="",ci=95,errwidth=1,capsize=0.1,scale=0.3,markers='x')
#sns.pointplot('condition','mumax',linestyles='',data=mumaxparams.sort_values('mumax',ascending=False),errwidth=3,ci=95,scale=1,markers=['_']*len(mumaxparams.condition.unique()))
#sns.stripplot(x='Comound',y='flux',hue='reactor',data=fluxdata)
## Plot fluxdata
#sns.catplot(data=fluxdata,x="Compound",y='flux',hue='reactor',col='strain',size=20)



#
### Split ccpa and wt experiment into different dataframes
#ccpa_locs = [rowname for rowname in hplcdata_raw.index if "ccpa" in rowname]
#
#wt_locs = [rowname for rowname in hplcdata_raw.index if "WT" in rowname]
#
#concentrations_ccpa = concentrations.loc[ccpa_locs,:]
#concentrations_ccpa = concentrations_ccpa.loc[concentrations_ccpa.index.get_level_values('reactor')!='n/a']
#concentrations_ccpa.index = concentrations_ccpa.index.droplevel('Sample name')
#concentrations_wt = concentrations.loc[wt_locs,:]
#concentrations_wt = concentrations_wt.loc[concentrations_wt.index.get_level_values('reactor')!='n/a']
#concentrations_wt.index = concentrations_wt.index.droplevel('Sample name')
#
#
#ccpa_reactors = concentrations_ccpa.groupby('reactor')
#ccpa_obsnum = len(concentrations_ccpa.index.get_level_values('timepoint').unique())
#ccpa_stderr = ccpa_reactors.std()/np.sqrt(ccpa_obsnum)
#ccpa_CI = t.ppf(1-0.025, (ccpa_obsnum-1)) * ccpa_stderr
##ccpa_reactors.mean().plot(kind='bar',yerr=ccpa_CI,title='ccpa',ylim=[-1,60])
#
#wt_reactors = concentrations_wt.groupby('reactor')
#wt_obsnum = len(concentrations_wt.index.get_level_values('timepoint').unique())
#wt_stderr = wt_reactors.std()/np.sqrt(wt_obsnum)
#wt_CI = t.ppf(1-0.025, (wt_obsnum-1)) * wt_stderr
##wt_reactors.mean().plot(kind='bar',yerr=wt_CI,title='wt',ylim=[-1,60])
#
#### CONVERT CONCENTRATIONS INTO FLUXES
#flux_ccpa = ccpa_reactors.mean().divide(dw_ccpa,axis=0) * MU  # mM/h/gDW
#flux_ccpa.plot(kind='bar', title='ccpa')#,ylim=[-1,35])
#flux_wt = wt_reactors.mean().divide(dw_wt,axis=0) * MU # mM/h/gDW
#flux_wt.plot(kind='bar', title='wt')#,ylim=[-1,35])
#
#plt.show()
