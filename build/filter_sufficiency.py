import pandas as pd

file_path = "C:\\Users\\solom\\Downloads\\Matched_Clingen_Data.csv"
df = pd.read_csv(file_path)

# Standardize column names
df.columns = df.columns.str.strip()

# Extract chromosome number from GRCh38 column
df['Chromosome'] = df['GRCh38'].astype(str).str.extract(r'chr(\d+)')

# Define chromosomes of interest
chromosomes_of_interest = {"2", "13", "3", "12", "4", "6", "14", "17", "11", "19"}

# Filter for Gain-of-Function (GoF) genes
GoF = df[df["Chromosome"].isin(chromosomes_of_interest)]

# Filter for Loss-of-Function (LoF) genes
LoF = df[
    (df["MOI_x"].str.contains("AR", na=False)) |
    (df["MOI_y"].str.contains("X-linked inheritance", na=False)) |
    (df["Haploinsufficiency"].str.contains('Gene Associated with Autosomal Recessive Phenotype', na=False))
]
LoF = LoF[LoF["Haploinsufficiency"] != "No Evidence for Haploinsufficiency"]

# Identify indices for LoF and GoF
lof_indices = set(LoF.index)
gof_indices = set(GoF.index)

# Find Dominant Negative (DN) by excluding GoF and LoF
dn_indices = set(df.index) - lof_indices - gof_indices
DN = df.loc[list(dn_indices)]


# Function to classify a gene symbol
def classify_gene(gene_symbol):
    gene_symbol = gene_symbol.upper()  # Normalize input to uppercase

    if gene_symbol in GoF["Gene Symbol"].values:
        return "GoF"
    elif gene_symbol in LoF["Gene Symbol"].values:
        return "LoF"
    elif gene_symbol in DN["Gene Symbol"].values:
        return "DN"
    else:
        return "Not classified"

#Function to assign classification values based on predefined rules
def assign_classification(classification_value):
    classification_map = {
        "Definitive": 1.0,
        "Strong": 0.9,
        "Medium": 0.5,
        "Weak": 0.3,
        "No Known Disease": 0.05,
        "Disputed": 0.05,
        "Refuted": 0.1
    }
    return classification_map.get(classification_value, None)

# Function to assign SOP values based on SOP category
def assign_sop(sop_value):
    if sop_value in ["SOP1", "SOP2", "SOP3", "SOP4"]:
        return 0.02
    elif sop_value in ["SOP5", "SOP6", "SOP7", "SOP8"]:
        return 0.05
    elif sop_value in ["SOP9", "SOP10", "SOP11"]:
        return 0.15
    else:
        return None  # Default case

# Function to build and update the matrix based on user-input disease
def build_matrix_for_disease_with_sop(disease_name):
    print(f"\nProcessing Disease: {disease_name}")

    # Check if required columns exist
    required_columns = ["Disease Label", "Gene Symbol", "SOP","Classification"]
    for col in required_columns:
        if col not in df.columns:
            print(f"Error: Missing column {col}")
            return

    # Filter dataset for the given disease
    disease_genes = df[df["Disease Label"].str.contains(disease_name, case=False, na=False)][["Gene Symbol", "SOP","Classification"]]

    if disease_genes.empty:
        print(f"No matching genes found for disease: {disease_name}")
        return

    # Initialize matrix list
    matrix_data = []

    # Iterate over genes and classify them
    for _, row in disease_genes.iterrows():
        gene = row["Gene Symbol"]
        sop_value = row["SOP"]
        classification_value = row["Classification"]


        # Assign SOP based on the SOP category
        classification = classify_gene(gene)
        sop_score = assign_sop(sop_value)
        classification_score = assign_classification(classification_value)

        # Define row structure
        matrix_row = {
            "Gene": gene,
            "Coefficient": classification_score if classification_score else 0.1,
            "SOP": sop_score,
            "LoF": "*" if classification == "LoF" else "-",
            "GoF": "*" if classification == "GoF" else "-",
            "DN": "*" if classification == "DN" else "-",
        }

        matrix_data.append(matrix_row)

    # Convert to DataFrame
    matrix_df = pd.DataFrame(matrix_data)
    
    return matrix_df

# Test with example disease names
def user_input_disease_classification():
    user_disease = input("Enter a disease name: ").strip()
    build_matrix_for_disease_with_sop(user_disease)

