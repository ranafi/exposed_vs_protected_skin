# Circadian Rhythms in Sun-Protected and Sun-Exposed Skin

This repo contains the R code associated with the manuscript *Comparative circadian transcriptome analysis reveals dampened and phase-advanced rhythms in sun-exposed human skin.* The data and results files are currently omitted due to privacy concerns, but will be added prior to publication of the paper.

The goal of this project was to detect and model circadian rhythms in gene expression, using time course bulk RNA-seq data from 20 subject, with skin biopsies taken from both the forearm (the *EXPOSED* site) and the buttocks (the *PROTECTED* site), at times 12 noon, 6PM, 12 midnight, and 6AM.

The scripts should be run in the following order:

1. **process_data.R** -- processes raw counts bulk RNA-seq data and stores linear-scale CPM expression matrices and sample information as RDS objects. Again, the data is currently excluded from this repo because of privacy concerns, but will be uploaded prior to publication.
2. **first_pass_cosinor.R** -- implements Cosinor regression analysis, applied to each condition (*PROTECTED* and *EXPOSED*) separately. Identifies rhythmic genes in each condition, and stores results (including fit amplitude and acrophase estimates for each rhythmic gene in each condition) as both RDS objects and .txt files in the *Results* folder.
3. **differential_cosinor.R** -- implements a differential Cosinor analysis to identify genes that were significantly differentially rhythmicm between the two conditions.

4. **make_pathway_analysis_inputs.R** -- prepares input files for PSEA, GSEA, and Enrichr pathway analyses.

5. **supplementary_tables.R** -- prepares supplementary tables, which are lists of genes and their rhythmic parameters that fall into several different categories. Combined supplementary table Excel spreadsheet is available in ```Results/supplementary_tables/supplementary_tables.xlsx```.



If you have any questions, please feel free to contact me by email at mikest (at) udel.edu.



