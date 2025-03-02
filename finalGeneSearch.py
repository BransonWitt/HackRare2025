from BCBio import GFF
import pandas as pd
from Bio import SeqIO
import re

def get_gene_id(specificString:str) -> int:
    """takes in a string from a gff file 

    Args:
        specificString (str): string from gff file

    Returns:
        int: gene id from the string of description from the gff file
    """
    #Regex match
    match = re.search(r'GeneID:(\d+)', specificString)

    #Applying regex
    if match:
        gene_id_entry = match.group(1)  # Extract the matched value (the number after GeneID:)
        return int(gene_id_entry)
    else:
        return None

def get_protein_specifics(specificString:str):
    #Gets the id for the protein
    data = specificString.split(';')
    
    protein_id_entry = next((entry for entry in data if entry.startswith("protein_id=")), None)
    
    #None handling, I might've messed up either way
    if protein_id_entry != None:
        protein_id_entry = protein_id_entry[11:]
    
    return protein_id_entry


def get_symbol(specificString):
    #Gets whats called the symbol from the NIH for the proteins
    data = specificString.split(';')
    
    symbol_id_entry = next((entry for entry in data if entry.startswith("gene=")), None)
    
    #None handling
    if symbol_id_entry != None:
        symbol_id_entry = symbol_id_entry[5:]
        
    return symbol_id_entry

def get_gene_name(specificString:str) -> str:
    #Gets the gene id for transferring between different formats
    data = specificString.split(';')
    
    gene_id_entry = next((entry for entry in data if entry.startswith("ID=")), None)
    
    #None handling
    if gene_id_entry != None:
        gene_id_entry = gene_id_entry[8:]
    
    return str(gene_id_entry)


def extract_CDS_from_gff(gff_file, fasta_file, requestedGeneID):
    """Extracts all CDS entries from a GFF file

    Args:
        gff_file (gff_file): gff file location
        fasta_file (fasta_file): fasta file location
        requestedGeneID (list): list of strings which correspond to geneIDs
    """
    # Load the genome sequence from the FASTA file
    genome_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    
    #Reading the GFF file into a pandas dataframe
    gff_df = pd.read_csv(gff_file, comment="#", delimiter="\t", header=None)
    
    #Searching through the dataframe to sort down to only entries that have to do with CDS 
    protein_coding = gff_df[gff_df[2] == 'CDS']
    gff_df = gff_df.reindex() #reindexing cut down version

    #Filling more columns
    protein_coding["protein_id"] = protein_coding[8].apply(get_protein_specifics)
    protein_coding["gene_ID"] = protein_coding[8].apply(get_gene_id)
    protein_coding["Symbol"] = protein_coding[8].apply(get_symbol)

    #Shortening down dataframe even further to the proteins within the specific gene ids requested
    shortened = protein_coding[protein_coding['gene_ID'].isin(requestedGeneID)]
    
    columns = ["Gene Location", "Protein Id", "Start", "End", "Sequence"]
    gene_information = []
    
    #Running through each entry in the severly cut down data frame into the refrenced position in the 
    for index, row in shortened.iterrows():
        start = row[3]
        end = row[4]
        gene = row[0]
        strand = row[6]
        protein_id = row["protein_id"]
        
        sequence = genome_dict[gene].seq[start:end]
        
        if strand == '-':
            sequence = sequence.reverse_complement()
        
        gene_information.append([gene, protein_id, start, end, str(sequence)])
    
    gene_df = pd.DataFrame(gene_information, columns=columns)
    
    return gene_df
        
def extract_genes_from_gff(geneCommonName:list, gff_file) -> list:
    """Takes the gff file and converts a list of geneCommonNames into the gene ids for the NIH

    Args:
        geneCommonName (list): list of the gommon internationally recognized gene names
        gff_file (_type_): gff file location

    Returns:
        list: list of 
    """
    #Reads the gff file into a pandas data frame
    gff_df = pd.read_csv(gff_file, comment="#", delimiter="\t", header=None)
    
    #Filtering out only the gene entries
    gene_coding = gff_df[gff_df[2] == 'gene']
    gene_coding = gene_coding.reindex() #reindex
    
    #Populating other columns using functions
    gene_coding['Name'] = gene_coding[8].apply(get_gene_name)
    gene_coding["gene_ID"] = gene_coding[8].apply(get_gene_id)
    
    #Filtering only the gene common names desired 
    gene_coding = gene_coding[gene_coding['Name'].isin(geneCommonName)]
    
    #returning a list of the NIH gene IDs for the common gene names desires when input into the function originally
    requestedGeneIDColumn = gene_coding.get('gene_ID')
    requestedGeneID = requestedGeneIDColumn.to_list()
    
    return requestedGeneID
            
# Example usage
gff_file = "C:\\Users\\brans\\Downloads\\Harvard_Hackathon\\genomic.gff"  # Your GFF file
fasta_file = "C:\\Users\\brans\\Downloads\\Harvard_Hackathon\\GCF_000001405.13_GRCh37_genomic.fna"  # Your FASTA file with the genome sequence
mutated_fasta_file = "C:\\Users\\brans\\Downloads\\Harvard_Hackathon\\modified_proteins.fna"

test_genes = ['CTSB', 'SCN1A', 'ABAT', 'DNM1']
#rGeneID = extract_genes_from_gff(test_genes, gff_file)

#print(rGeneID)
#extract_CDS_from_gff(gff_file, fasta_file, rGeneID)

def get_DNA_pct_change(DNA1, DNA2):
    #Assuming the same length of each DNA and checking for similarity or difference
    assert(len(DNA1) == len(DNA2))
    
    print(DNA1)
    print(DNA2)
    print( )
    
    mutations = 0
    
    for i in range(len(DNA1)):
        if DNA1[i] != DNA2[i]:
            
            mutations += 1
    
    
    return (mutations / len(DNA1))


def compare_genes(genes_to_look_at, gff_file1, fasta_file1, gff_file2, fasta_file2):
    
    #Transforming genes from international code like "ABAT" to a NIH Code like 1759
    requested_genes = extract_genes_from_gff(genes_to_look_at, gff_file1)
    
    proteinDNA1 = extract_CDS_from_gff(gff_file1, fasta_file1, requested_genes)
    proteinDNA2 = extract_CDS_from_gff(gff_file2, fasta_file2, requested_genes)
    
    new_columns = ["Gene Location Reference", "Protein ID Reference", "Start Reference", "End Reference", "Sequence Reference", "Gene Location Mutant", "Protein ID Mutant", "Start Mutant", "End Mutant", "Sequence Mutant"]
    final_df =  pd.concat([proteinDNA1, proteinDNA2], axis=1).reindex(proteinDNA1.index)
    
    final_df = final_df.set_axis(new_columns, axis=1)
    
    final_df['pct_change'] = final_df.apply(
        lambda row: get_DNA_pct_change(row['Sequence Reference'], row['Sequence Mutant']),
        axis=1
    )
    
    print(final_df[['Sequence Mutant', 'Sequence Reference', 'pct_change']].describe())
    return final_df


compare_genes(test_genes, gff_file, fasta_file, gff_file, mutated_fasta_file)



#Sample Output
"""
Feature: NP_001005484.1
Location: 69091 - 70008
Sequence: TGGTGACTGAATTCATTTTTCTGGGTCTCTCTGATTCTCAGGAACTCCAGACCTTcctatttatgttgttttttgta
TTCTATGGAGGAATCGTGTTTGGAAACCTTCTTATTGTCATAACAGTGGTATCTGACTCCCACCTTCACTCTCCCATGTACTTCCTG
CTAGCCAACCTCTCACTCATTGATCTGTCTCTGTCTTCAGTCACAGCCCCCAAGATGATTACTGACTTTTTCAGCCAGCGCAAAGTC
ATCTCTTTCAAGGGCTGCCTTGTTCagatatttctccttcacttctttgGTGGGAGTGAGATGGTGATCCTCATAGCCATGGGCTTT
GACAGATATATAGCAATATGCAAGCCCCTACACTACACTACAATTATGTGTGGCAACGCATGTGTCGGCATTATGGCTGTCACATGG
GGAATTGGCTTTCTCCATTCGGTGAGCCAGTTGGCGTTTGCCGTGCACTTACTCTTCTGTGGTCCCAATGAGGTCGATAGTTTTTAT
TGTGACCTTCCTAGGGTAATCAAACTTGCCTGTACAGATACCTACAGGCTAGATATTATGGTCATTGCTAACAGTGGTGTGCTCACT
GTGTGTTCTTTTGTTCTTCTAATCATCTCATACACTATCATCCTAATGACCATCCAGCATCGCCCTTTAGATAAGTCGTCCAAAGCT
CTGTCCACTTTGACTGCTCACATTACAGTAGTTCTTTTGTTCTTTGGACCATGTGTCTTTATTTATGCCTGGCCATTCCCCATCAAG
TCATTAGATAAATTCCTTGCTGTATTTTATTCTGTGATCACCCCTCTCTTGAACCCAATTATATACACACTGAGGAACAAAGACATG
AAGACGGCAATAAGACAGCTGAGAAAATGGGATGCACATTCTAGTGTAAAGTTTTAG
"""