# MScR-scRNA-seq-analysis

Scripts:

1. main.R: The main script for data analysis using the no-pulsatile RDS file
2. old_Main.R: The main script for data analysis usign the orginial RDS file (containing pulsatile)
3. annotation.R: For annotation purposes
4. cell-typing.R: FeaturePlot + DotPlots using validated marker genes (from Fidanza et al. 2020. Blood)
5. main_rescaled: The main script for data analysis usign the orginial RDS file (containing pulsatile), remove the omitted features from the data and rescaled. Do not use, it does not make sense
6. add-quotation-marks.R: Script that will add quotation around the gene names in xsl files
7. trajectory-analysis.R: Trajecotry analysis for EHT + haematopoietic (cluster 7 + 6), including functions for plotting features over pseudotime (all coloured/ two-colour only)
8. determine-resolution.R: Script to determine which resolution to use
9. metabolic-gene-expression.R: Script to visualise metabolic-assocaited gene lists obtained from literature review
10. 240614-merge-rep1,rep2-myself.R: Script to perform harmony on rep1 + rep2 RDS prepared by AF (no pulsatile)
11. automated-removing-omitted-features.R: Remove gene list from the csv files and save the new file accordingly
12. subset-clusters-for-further-analysis.R: Script to subset the clusters
13. heatmap-on-top10-filtered-out-omitted.R: Heatmap visualisation
14. MiloR-speeded.R: Updated speed MiloR analysis by the dev
15. go-TF+membrane-marker.R: Uses GO term to identify potential TFs and MMs
