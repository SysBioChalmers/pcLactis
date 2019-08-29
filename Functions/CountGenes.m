%% CountGenes
function Gene_list = CountGenes(M_Model_total)

Gene_list = struct();

% Count M model genes.
Gene_list.M_gene = M_Model_total.genes;

% Count E model genes.
Gene_list.E_protein_gene = importdata('E_model_protein_gene_list.txt');
Gene_list.E_rRNA_gene = importdata('E_model_rRNA_list.txt');
Gene_list.E_tRNA_gene = importdata('E_model_tRNA_list.txt');
Gene_list.E_total_gene = unique([Gene_list.E_protein_gene;
                               Gene_list.E_rRNA_gene;
                               Gene_list.E_tRNA_gene]);

% Total genes.
Gene_list.total_gene = unique([Gene_list.M_gene;Gene_list.E_total_gene]);

% Total RNA genes.
Gene_list.E_RNA_gene = unique(cat(1,Gene_list.E_rRNA_gene,Gene_list.E_tRNA_gene));

% Total protein genes, will be used in protein stoichiometry.
Gene_list.total_protein_gene = unique([Gene_list.M_gene;Gene_list.E_protein_gene]);

