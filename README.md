# Major features of parasite adaptation revealed by genomes of Plasmodium falciparum population samples archived for over 50 years

## Overview
This repository contains the data processing and analysis workflows for a genomic study of \*Plasmodium falciparum\* malaria parasites from archived placental blood samples collected in The Gambia between 1966 and 1971. The project investigates the evolutionary history of the parasite prior to widespread antimalarial drug pressure, identifying signatures of selection and temporal genomic changes over a ~50-year period.

## **Key Findings:**
- No drug resistance alleles were detected in the historic (1966–1971) parasite population.
- Strong signatures of positive selection were identified at multiple loci, particularly in genes encoding surface proteins involved in erythrocyte binding and antibody recognition.
- Major allele frequency changes over 50 years were observed at loci including *gdv1* (sexual conversion regulation), *Pfsa1*, and *Pfsa3* (associated with sickle-cell trait infection).
- Historic infections showed high genomic complexity and were largely genetically unrelated.

## Project Structure

 ```
├── data/ genomic data (VCFs, filtered SNPs)
    └── metadata
    └── processed
    └── raw
├── scripts/                        # Analysis scripts (R, shell)
├── results/                        # Output figures, tables, and statistical summaries
    └── Figures
    └── Tables
├── docs/                           # Supplementary documentation and methods details
└── README.md                       # This file
```
## Methods Summary

### 1. Sample Identification and DNA Extraction
- **Source:** 61 lyophilised placental blood samples from routine deliveries in Greater Banjul, The Gambia (1966–1971).
- **Ethics:** Approved by the Joint Ethics Committee of the MRC Unit and the Gambian Government; samples were anonymised.
- **Extraction:** DNA was extracted using the QIAamp® Blood Extraction Midi kit (QIAGEN).

### 2. Sequencing and Variant Calling
- **Library Prep:** Selective whole genome amplification to enrich *P. falciparum* DNA.
- **Sequencing:** Illumina HiSeq, 150bp paired-end reads (Wellcome Sanger Institute, MalariaGEN pipeline).
- **Variant Calling:** MalariaGEN Pf7 pipeline; 49,565 coding SNPs detected, filtered to 35,998 high-quality bi-allelic SNPs.
- **Final Dataset:** 54 high-coverage historic samples retained for analysis.

### 3. Comparative Temporal Samples
- **2015 Samples:** 89 *P. falciparum* -infected blood samples from Brikama Health Centre (The Gambia).
- **Historic Comparators:** Previously sequenced samples from 1984 and 2001 from the same region.
- **Uniform Processing:** All samples processed and variant-called using the same pipeline.

### 4. Statistical & Population Genetic Analyses

#### A. Genomic Complexity of Infections
- **Within-infection fixation index (FWS):** Calculated using `moimix` R package (3,927 SNPs, MAF ≥5%).
- **Complexity of Infection (COI):** Estimated via the REAL McCOIL method (`McCOILR`).

#### B. Population Structure and Relatedness
- **SNP Pruning:** 6,267 low-LD SNPs selected (sliding window, r² > 0.1).
- **Multidimensional Scaling (MDS):** Performed using `plink` and visualized in R.
- **Pairwise Relatedness:** Identity-by-state (IBS) distances displayed as heatmaps (`Pheatmap`).

#### C. Identifying Genomic Regions Under Selection (1966–1971)
- **Balancing Selection:** Beta score (β) via `BetaScan`; Tajima’s D for 3,038 genes.
- **Recent Positive Selection:** Integrated Haplotype Score (iHS) calculated with `Rehh` (12,416 SNPs, MAF ≥0.03).

#### D. Temporal Genomic Changes (1966–1971 vs. 2015)
- **Temporal Differentiation:** FST (Nei’s GST) calculated for 12,416 SNPs (`mmod` R package).
- **Visualization:** Manhattan and zoomed locus plots (`qqman`, `locuszoom`).
- **Extended Timeline:** Allele frequency trajectories analysed across 1984, 2001, and 2015 samples.

## Data Availability
- **Sequence Data:** Available at the European Nucleotide Archive (ENA). Accession numbers and coverage details are provided in Supplementary Table S1.
- **VCF Files:** Processed variant calls for historic and contemporary samples are available upon request.
- **Scripts:** Custom R and Python scripts for reproduction of analyses are provided in the `scripts/` directory.

## Software & Dependencies
Key tools and packages used in this analysis:
- **Variant Calling & Filtering:** MalariaGEN Pf7 pipeline, VCFtools
- **Population Genetics:** `plink`, `moimix`, `Rehh`, `BetaScan`, `mmod`
- **Visualization:** `ggplot2`, `Pheatmap`, `qqman`, `locuszoom`
- **General Analysis:** R (≥4.0), Python (≥3.8)

Installation instructions for specific packages are provided within the individual scripts.

## Reproducibility
All analysis steps are documented in the `scripts/` directory with commented code. Key parameters and filtering thresholds are detailed in the Methods section above and within the scripts. Sample sizes and SNP sets for each analysis are explicitly noted in the relevant script headers.

## Citation
If you use this data or code, please cite the corresponding publication:
[Publication details will be added upon acceptance]

## Acknowledgements
- Sample collection and archiving: MRC Unit, The Gambia.
- Sequencing and data processing: Wellcome Sanger Institute (MalariaGEN pipeline).

## - Funding: [Details to be added].

## Contact
For questions regarding data or analysis, please contact:
- **Principal Investigator:** [Alfred Amambua-Ngwa, MRC Unit The Gambia at London School of Hygiene and Tropical Medicine, Banjul, The Gambia, alfred.ngwa@lshtm.ac.uk]
- **Principal Investigator:** [David J. Conway, Department of Infection Biology, London School of Hygiene and Tropical Medicine, London, UK, david.conway@lshtm.ac.uk]
- **Bioinformatics Lead:** [Mouhamadou Fadel Diop, MRC Unit The Gambia at London School of Hygiene and Tropical Medicine, Banjul, The Gambia, mdiop@mrc.gm]

## License
This project is licensed under the [MIT License] – see the LICENSE file for details.

