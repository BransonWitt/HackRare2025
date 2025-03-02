from filter_sufficiency import *
from  finalGeneSearch import *
import pandas as pd

gff_ref = "genomic.gff"
fna_ref = "GCF_000001405.13_GRCh37_genomic.fna"

sample_mut = "modified_proteins.fna"

def get_median_mut(gene, df):
    column = df[df["gene_ID1"] == gene]['pct_change']
    return column.median()

def get_average_mut(gene, df):
    column = df[df["gene_ID1"] == gene]['pct_change']
    
    return column.mean()

def matrix(disease, gff_mut, fna_mut, gff_reference = gff_ref, fna_reference = fna_ref):
    building_matrix = build_matrix_for_disease_with_sop(disease)
    
    building_matrix = pd.DataFrame(building_matrix)
    building_matrix = building_matrix.drop_duplicates().reset_index(drop=True)
    
    targeted_genes = building_matrix['Gene'].dropna().unique()
    
    coded_to_target = genes_name_dict(targeted_genes, gff_reference)
    
    large_matrix = compare_genes(gff_reference, fna_reference, gff_mut, fna_mut, genes_to_look_at=targeted_genes)
    
    large_matrix = large_matrix.replace({"gene_ID1":coded_to_target})
        
    building_matrix["Median Mut"] = building_matrix['Gene'].apply(lambda x: get_median_mut(x, large_matrix))
    building_matrix["Mean Mut"] = building_matrix['Gene'].apply(lambda x: get_average_mut(x, large_matrix))
    
    return(building_matrix)
    
    
matrix('developmental and epileptic encephalopathy', gff_ref, sample_mut)