# **Gen-IE**

## **Overview**
**Gen-ie** is a revolutionary genomic sequencing tool designed to facilitate early diagnosis of rare diseases. By analyzing **500 rare disease genes**, it identifies mutations in nitrogenous bases within DNA that contribute to disease development. **Gen-ie** processes the entire human genome, which consists of approximately **3 billion base pairs**, delivering precise insights into the genetic basis of rare diseases. This enables faster and more accurate diagnoses, significantly reducing the lengthy and challenging journey patients often face in obtaining a diagnosis.

## **Background & Research**
Whole-genome sequencing (**WGS**) has demonstrated transformative potential in rare disease research. A national health system study involving **13,037 participants** identified **95 Mendelian gene-disease associations**, confirming **79 as causal**. The study also uncovered **novel non-coding variants** that contribute to disease by disrupting critical genes such as *ARPC1B, GATA1, LRBA,* and *MPL*. These findings emphasize the power of **WGS** in advancing diagnostic and etiological discovery in clinical settings (**Turro et al., 2020**). 

## **Data Analysis & Scoring Matrix**
**Gen-ie** utilizes the **ClinGen API** for individualized research on proteins and mutated genes. Key analytical parameters include:
- **Mutation location**: Specific chromosome regions
- **Haploinsufficiency & triplosensitivity**: Indicators of gene dosage sensitivity
- **Hereditary inheritance patterns**: *X-linked, autosomal dominant/recessive*

To quantify relevance, a **scoring matrix (1-10)** is implemented with color gradients based on the **ClinGen API** data. The relevance metric follows **standard operating procedures (SOPs):**
- **Definitive (1.0)**, **Strong (0.9)**, **Moderate (0.5)**, **Weak (0.3)**
- **No known disease (0.05)**, **Disputed (0.05)**, **Refuted (0.1)**

Additionally, gene mutations are categorized based on their mechanisms:
- **Loss of Function (LOF)**
- **Gain of Function (GOF)**
- **Dominant Negative (DN)**

### **Key Findings**
Through **UniProt** analysis:
- **LOF mutations** correlate with **positive haploinsufficiency**, *X-linked recessive,* and *autosomal recessive* inheritance.
- **GOF mutations** are associated with *autosomal dominant inheritance*, **positive triplosensitivity**, and mutations in chromosomal regions *2, 3, 6, 12, 13, 14, 17, and 19*.
- **DN mutations** were assigned in cases where clinical research was insufficient.

## **Functionality**
**Gen-ie** enables users to input a genome sequence in **FASTA format**. The system then:
1. **Maps the sequence against a database** to detect mutations in **500 rare disease genes**.
2. **Classifies genes based on mutation type** (*GOF, LOF, or DN*).
3. **Generates a relevance score** based on published research data.

Currently, the system operates at an estimated **70% efficiency** in sorting gene data functions, with the remaining **30% limited by data scarcity**. Future enhancements will focus on expanding research on understudied genes and incorporating predictive tools like **PolyPhen-2** to assess amino acid substitution effects on protein function.

## **Technical Stack**
**Gen-ie** is built using:
- **Python** for coding and data analysis
- **ClinGen API** for gene relevance and inheritance research
- **Reactome** for pathway identification
- **PubMed & Google Scholar** for academic research
- **UniProt & NIH resources** for protein and gene function analysis
- **Excel** for gene mutation categorization

## **Future Improvements**
- Enhancing data processing efficiency through broader gene research
- Implementing **PolyPhen-2** to improve mutation impact prediction
- Refining the user interface for better accessibility and efficiency

## **Conclusion**
**Gen-ie** represents a groundbreaking step in rare disease diagnosis, merging advanced genomic sequencing with a powerful analytical framework. By providing an **efficient, user-friendly, and research-oriented solution**, **Gen-ie** accelerates medical advancements and enhances patient care. This project underscores the importance of **collaboration**, leveraging diverse tools and research methodologies to drive innovation in genetic diagnostics.

