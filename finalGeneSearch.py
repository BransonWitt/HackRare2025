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
    
    epilepsy_genes = [['Gene', 'Start', 'End']]
    
    #Running through each entry in the severly cut down data frame into the refrenced position in the 
    for index, row in shortened.iterrows():
        start = row[3]
        end = row[4]
        gene = row[0]
        strand = row[6]
        protein_id = row["protein_id"]
        
        sequence = genome_dict[gene].seq[start:end]
        
        epilepsy_genes.append([gene, start, end])
        
        if strand == '-':
            sequence = sequence.reverse_complement()
        
        epilepsy_df = pd.DataFrame(epilepsy_genes)
        epilepsy_df.to_csv("C:\\Users\\brans\\Downloads\\ncbi_dataset\\ncbi_dataset\\data\\GCF_000001405.13\\epilepsy.csv")
        
        # Print or return the sequence
        print(f"Feature: {gene}")
        print(f"Location: {start} - {end}")
        print(f"Sequence: {sequence}")
        print()  # Newline for separation
        
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
gff_file = "genomic.gff"  # Your GFF file
fasta_file = "GCF_000001405.13_GRCh37_genomic.fna"  # Your FASTA file with the genome sequence

test_genes = ['CTSB', 'SCN1A', 'ABAT', 'DNM1']
rGeneID = extract_genes_from_gff(test_genes, gff_file)

print(rGeneID)
extract_CDS_from_gff(gff_file, fasta_file, rGeneID)













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