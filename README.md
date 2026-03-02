
# Differential gene expression analysis of _S. cerevisiae_ at different stages of velum formation 
## Using data from Mardanov et al. 2020

### for BINF 6110 Assignment 2
### March 1st 2026
### by Eva Burguete-Innocente


## Introduction

The fermentation of sherry wine is facilitated by flor yeasts _(Saccharomyces cerevisiae)_ that form a biofilm or velum that covers the surface of the wine (Mardanov et al., 2023). These yeasts are able to tolerate ethanol and contain many specific genetic variations (Mardanov et al., 2023). The velum formation occurs when nitrogen and sugar resources are diminished after alcoholic fermentation, so the yeast switch to oxidative metabolism (Legras et al., 2016). The cells aggregate together and rise to the surface, eventually building up to form the velum (Legras et al., 2016). 

RNA-sequencing data analysis is highly important across biological contexts (Rosati et al., 2024). One of the primary steps is aligning the RNA-seq reads to a reference genome or reference transcriptome, followed by counting of the sequenced reads, which can be achieved with pseudoalignment tools (Corchete et al., 2020; Patro, Mount, and Kingsford, 2014). Pseudoaligners count k-mers occurring in reads rather than aligning each base to a reference (Corchete et al., 2020; Patro, Mount, and Kingsford, 2014). Some common pseudoaligners are Kallisto (Bray et al., 2016), Sailfish (Patro, Mount, and Kingsford, 2014), and Salmon (Patro et al., 2017). All three methods vastly improved the speed of RNA-seq alignment and quantification compared to traditional methods, but Salmon is the most accurate of the three (Patro et al., 2017) and the most precise (Corchete et al., 2020). Salmon takes into account sample-specific bias and was shown to have a better performance in terms of correctly quantifying differential expressed genes (Patro et al., 2017). For this reason, I will use Salmon to align and quantify the RNA-seq reads. 

After read alignment and quantification, differential gene expression analysis is often performed, a statistical quantification of differentially expressed genes between experimental conditions (Corchete et al., 2020). A wide variety of tools exist for this purpose, such as edgeR (Robinson, McCarthy, & Smyth, 2010), DEseq2 (Love, Hubers, & Anders, 2014), and limma (Ritchie et al., 2015), to name only a few of the most popular. Several reviews have shown that DESeq2 is able to accurate identify differential expressed genes, and is particularly useful for small datasets (Love, Hubers, & Anders, 2014; Rapaport et al., 2013; Seyednasrollah, Laiho, & Elo, 2015). However, authors emphasise that no single tool will be perfectly suited to any analysis (Rapaport et al., 2013; Seyednasrollah, Laiho, & Elo, 2015). Given its good performance, ease of use (a quickly installed R package), and the fact that there are many informative and detailed vignettes available for DEseq2 (such as Love, Hubers, & Anders, 2025), I will use it in my analysis.

Functional enrichment analysis of the differentially expressed genes is necessary to identify any significantly enriched biological pathways and assess the biological importance of the significant genes (Rosati et al., 2024). This can be done via over-representation analysis (ORA) or gene set enrichment analysis (GSEA), of which GSEA has a slight advantage, as rather than applying a strict cutoff, the list of genes is ranked (Goeman & Bühlmann, 2007). Gene Ontology (GO) and KEGG databases contain comprehensive information about gene sets, and are commonly used (Geistlinger et al., 2021; Mardanov et al., 2020). Again, many tools exist for this purpose, like simplifyEnrichment (Gu & Hübschmann, 2023), GeneFEAST (Taylor et al., 2025) and clusterProfiler (Yu et al., 2012). ClusterProfiler is a widely used R package which allows GSEA to be perfomed with GO and KEGG annotations (Wu et al., 2021). As a comprehensive tutorial exists for this package, I used clusterProfiler to perform functional enrichment analysis (Khalfan, 2019). 

Mardanov et al. (2020) compared gene expression levels of _S. cerevisiae_ at three stages of velum formation (early, thin, and mature) using transcriptome data from experimental conditions simulating biological wine aging. They characterised gene expression changes between stages and identified which functional categories were most enriched. I will use their RNA-sequencing reads to perform differential gene expression analysis using the best-performing tools of my choice as specified above.

## Methods

I downloaded the raw RNA-seq reads from NCBI SRA (accessions: SRR10551657-SRR10551665; Mardanov et al., 2020). I used FastQC v0.12.1 to visualise read quality (Andrews, 2010). I used Salmon v1.10.3 to build an index for the reference transcriptome and quantify the transcript counts for each sample (Patro et al., 2017). I specified -k 31 for building the index, as this builds the index over k-mers of length 31, which Patro et al. (2017) found works well for reads that are 75bp or longer. However, I later realised the reads are only 50 base pairs long, and I did not have time before the submission deadline to go back and change this parameter- this is a limitation of the results and in future work, I would reduce the -k flag to reflect the shorter read length, as per the [Salmon documentation](https://salmon.readthedocs.io/en/latest/salmon.html). For quantifying the transcript counts, I used four threads and the -gcBias flag as recommended by Love, Anders, and Huber in the DEseq2 vignette (2025). This flag estimates a correction factor for common biases in RNA-seq data (Love, Hogenesch, and Irizarry, 2016; Patro et al., 2017). Both FastQC and Salmon were installed in conda environments and analyses run on a virtual machine running Ubuntu 25.10 with 8 cores and 12 GB of RAM. Code used for these analyses can be found [here.](./01_Salmon_analyses).

I used the _S. cerevisiae_ strain S288C transcriptome (GCF_000146045.2) as the reference for building the index, not _S. cerevisiae_ strain I-329 (GenBank PTER00000000) as Mardanov et al. (2020) used. This is due to the fact that I was having trouble matching the transcript annotations when analysing differential gene expression with DEseq2, so I had to use a different reference transcriptome. 
For analysing the differential gene expression between different stages of velum, I used DEseq2 (Love, Hubers, & Anders, 2014). I imported the count data with the tximport v1.34.0 (Soneson, Love, & Robinson, 2016). To assign gene IDs to the transcripts, I used the R packages txdbmaker v1.2.1 (Pagès et al., 2025) and AnnotationDbi v1.68.0 (Pagès et al., 2025b) to make a database of yeast gene IDs from the reference genome .gff file. I performed differential gene expression analysis following the DEseq2 vignette by Love, Hubers, and Anders (2025). I prefiltered for genes with low counts as recommended, and computed results for each contrast between stages of velum formation (i.e. Mature vs. Thin, Mature vs. Early, Thin vs. Early). I performed log-fold shrinkage for ease of visualisation and plotted the significantly differentially expressed genes for each contrast. I also applied the variance stabilising transformation as per Anders and Huber (2010), with the blind parameter set to True, as experimental design is not expected to influence the counts. I plotted heatmaps of the counts, performed a PCA, and assessed dispersion of the counts. 

I then performed a GO term and KEGG pathway gene set enrichment analysis using the packages clusterProfiler v4.14.6 (Yu et al., 2012), DOSE v4.0.1 (Yu et al., 2015) and enrichplot v1.26.6 (Yu, 2025). I used the _S. cerevisiae_ genome annotation package to provide the GO and KEGG annotations (Carlson, 2024). I ran the GSEA with the following parameters: search all GO/KEGG categories, gene name as the key type, minimum genes in set = 2, maximum genes in set = 1000, p-value cutoff = 0.05, database as yeast reference genome database, and the Benjamini-Hochberg correction for multiple testing. Both the diffrenetial expression analyses and gene set enrichment were run in RStudio v4.4.2 (RStudio team, 2025). All code for these analyses is [here.](./R/02_DE_and_GSEA.Rmd)

## Results

Overall sequence was quality was good according to FastQC. There were elevated levels of overrepresented sequences, but I did not trim these as is recommended for differential gene expression analysis (Liao & Shi, 2020). Mapping rate to the salmon index ranged from 74.22-92.28% for the nine samples. 

When comparing the thin vs. mature stages, 1177 genes had a significant positive log-fold change and 1248 genes had a significant negative log-fold change (adjusted p-value < 0.05) (Fig. 1). The top five most significantly differentially expressed genes were TKL2, AFR1, YPR127W, DSF1, and ADI1 (all with a negative log-fold change). All significantly expressed genes for this comparison can be found in [Supplemental Data 1.](./data/results/thin_vs_mature_shrink_results.csv)

<img width="875" height="540" alt="fig1" src="https://github.com/user-attachments/assets/b698923f-5cb4-4d35-956a-926925ee9be6" />

Figure 1. Volcano plot showing the differentially expressed genes between the thin and mature stages, where significant genes with a negative log-fold change are blue and significant genes with a positive log-fold change are red. 

When comparing the mature vs. early stages, 1561 genes had a significant positive log-fold change, and 1463 had a significant negative log-fold change (adjusted p-value < 0.05) (Fig. 2). The top five most significantly differentially expressed genes were TDH1 (negative log-fold change), FLO11 (positive), OLE1 (negative), PDC6 (negative), and HXT1 (negative). All significantly expressed genes for this comparison can be found in [Supplemental Data 2.](./data/results/mature_vs_early_shrink_results.csv)

<img width="875" height="540" alt="fig2" src="https://github.com/user-attachments/assets/f91d5172-7e9b-4db1-9bac-ff5408ace39c" />

Figure 2. Volcano plot showing the differentially expressed genes between the mature and early stages, where significant genes with a negative log-fold change are blue and significant genes with a positive log-fold change are red. 

When comparing the thin vs. early stages, 1158 genes had a significant positive log-fold change, and 1115 had a significant negative log-fold change (adjusted p-value < 0.05) (Fig. 3). The top five most significantly differentially expressed genes were PDC6 (negative), ADH7 (positive), HXT1 (negative), PCK1 (positive), and CTT1 (negative). All significantly expressed genes for this comparison can be found in [Supplemental Data 3.](./data/results/thin_vs_early_shrink_results.csv)

<img width="875" height="540" alt="fig3" src="https://github.com/user-attachments/assets/6b66bafe-2371-4b57-8c2f-d7e1834e9fc7" />

Figure 3. Volcano plot showing the differentially expressed genes between the thin and early stages, where significant genes with a negative log-fold change are blue and significant genes with a positive log-fold change are red. 

Plotting a PCA of the variance-stabilised transformed data (Fig. 4) revealed that the data clustered well by stage. Expression of the top 20 most highly expressed genes across samples and stages can be seen in Figure 5. To highlight one notable gene, EGO4 was downregulated in the early stage compared to the other stages.

<img width="875" height="540" alt="fig4" src="https://github.com/user-attachments/assets/9ad4ea74-5135-4b6c-a903-858a422dba96" />

Figure 4. PCA plot of the variance-stabilised transformed data, coloured by stage.

<img width="875" height="540" alt="fig5" src="https://github.com/user-attachments/assets/d5ea39ea-730f-4a60-822a-61ce193dfda3" />

Figure 5. Heatmap of the top 20 most highly expressed genes across samples and stages of velum formation, where warmer colours indicate upregulation and cooler colours indicate downregulation.

Gene set enrichment analysis based on GO terms indicated that for thin vs. mature stages, the GO term for the gene set with the lowest adjusted p-value was “energy reserve metabolic process”, containing 35 genes. For KEGG pathways, the pathway with the lowest adjusted p-value was “metabolic pathways” with 753 genes. There were 79 significantly enriched GO terms ([Supplemental Data 5,](./data/gsea_results/gse_t_vs_m_results.csv) Fig. 6) and 13 KEGG pathways ([Supplemental Data 6,](./data/gsea_results/t_vs_m_kegg_gsea_results.csv) Fig. 7). 

<img width="875" height="540" alt="fig6" src="https://github.com/user-attachments/assets/61ab2412-c5c1-409c-b050-f7dbecef3696" />

Figure 6. Ridgeplot showing the distribution of log-fold changes of the most enriched GO terms between thin and mature stages, where warmer colours indicate adjusted lower p-values, and cooler colours indicate higher p-values. 

<img width="875" height="540" alt="fig7" src="https://github.com/user-attachments/assets/6d32a93f-83ae-40cd-928a-c52367cbef06" />

Figure 7. Ridgeplot showing the distribution of log-fold changes of the most enriched KEGG pathways between thin and mature stages, where warmer colours indicate lower adjusted p-values, and cooler colours indicate higher p-values. 

For mature vs. early stages, the GO term with the lowest adjusted p-value was “structural constituent of ribosome” with 223 genes. The KEGG pathway with the lowest adjusted p-value was “ribosome” with 204 genes. There were 227 significantly enriched GO terms ([Supplemental Data 7,](./data/gsea_results/gse_m_vs_e_results.csv) Fig. 8) and 20 KEGG pathways ([Supplemental Data 8,](./data/gsea_results/m_vs_e_kegg_gsea_results.csv) Fig. 9). 

<img width="875" height="540" alt="fig8" src="https://github.com/user-attachments/assets/a9f90c32-4f88-4f1a-81c3-0fb111937ce3" />

Figure 8. Ridgeplot showing the distribution of log-fold changes of the most enriched GO terms between mature and early stages, where warmer colours indicate adjusted lower p-values, and cooler colours indicate higher p-values. 

<img width="875" height="540" alt="fig9" src="https://github.com/user-attachments/assets/f1edc8ca-dcd9-45b3-840a-a01de7de1a7e" />

Figure 9. Ridgeplot showing the distribution of log-fold changes of the most enriched KEGG pathways between mature and early stages, where warmer colours indicate lower adjusted p-values, and cooler colours indicate higher p-values. 

Finally, for thin vs. early stages, the GO term with the lowest adjusted p-value was “organellar ribosome” with 86 genes. The KEGG pathway with the lowest adjusted p-value was “ribosome” with 204 genes. There were 238 significantly enriched GO terms ([Supplemental Data 9,](./data/gsea_results/gse_t_vs_e_results.csv) Fig.10) and 18 KEGG pathways ([Supplemental Data 10,](./data/gsea_results/t_vs_e_kegg_gsea_results.csv) Fig. 11). 

<img width="875" height="540" alt="fig10" src="https://github.com/user-attachments/assets/7507b5f2-9ba8-4e16-9d8d-1bb7ccc64a5a" />

Figure 10. Ridgeplot showing the distribution of log-fold changes of the most enriched GO terms between thin and early stages, where warmer colours indicate adjusted lower p-values, and cooler colours indicate higher p-values. 

<img width="875" height="540" alt="fig11" src="https://github.com/user-attachments/assets/2ed1b448-32f2-41a5-9357-c0097996500f" />

Figure 11. Ridgeplot showing the distribution of log-fold changes of the most enriched KEGG pathways between thin and early stages, where warmer colours indicate lower adjusted p-values, and cooler colours indicate higher p-values. 

## Discussion

As mentioned above, I specified the -k 31 parameter incorrectly when creating the index with salmon. A smaller k-mer length may have improved mapping rate and thus affected downstream predictions, and is something to note for future analysis (Patro et al., 2017). 

There were thousands of genes significantly differentially expressed when comparing each stage of velum formation, and many functionally enriched pathways as well. For brevity, I chose to focus on one of the top five most significant genes for each comparison that was also identified by Mardanov et al. (2020), to increase confidence in my findings. 

I found that TKL2 was significantly downregulated in the thin stage compared to the mature stage (log-fold change = -4.59, adjusted p-value < 0.05) as did Mardanov et al. (2020). TKL2 is a transketolase that is part of the pentose phosphate pathway (Schaaff‐Gerstenschläger et al., 1993). Interestingly, in a study of differential gene expression between the aerial and root stages of a yeast infection biofilm, Maršíková et al. (2017) found that TKL2 was upregulated in the aerial stage, highlighting the potential of TKL2 as a target for future research on yeast biofilm formation. 

When comparing mature vs. early stages, FLO11 was significantly upregulated in mature stages, with a log-fold change of 5.48 (adjusted p-value < 0.05). FLO11 is known to be one of the most important genes involved in velum formation, as it is an adhesin involved in both cell-cell adhesion and cell-substrate adhesion (known as flocculation) (Mardanov et al., 2020; Zara et al., 2009). Finding that this gene was significantly upregulated in the mature velum stage corroborates existing evidence and indicates its critical role in velum formation.

Between the thin and early stages, PDC6 was downregulated in the thin stage (log-fold change = -4.78, adjusted p-value < 0.05). PDC6 is a pyruvate decarboxylase, which are isoenzymes involved in alcoholic fermentation (De Assis et al., 2013; Flikweert et al., 1996). Mardanov et al. (2020) also identified downregulation of PDC6 in their analyses. As the switch from alcoholic fermentation to oxidative metabolism kickstarts velum formation, this finding makes sense, as expression of genes involved in alcoholic fermentation would no longer be required when the velum has reached the thin stage. 

There were quite a few more significantly enriched GO terms in the mature and thin stages compared to the early stage, potentially indicating that more changes in cell function occur as the velum stages become more advanced. To highlight one enriched biological process out of many, “glycogen metabolic process” was significantly negatively enriched (i.e. downregulated) in the thin stage compared to the mature stage (Fig. 6). This likely reflects the fact that as the velum reaches the mature stage, the yeast has already switched over to oxidative metabolism, and there are no more sugar resources available to cells. Therefore, the metabolism of glycogen being downregulated indicates that there is no more glycogen available to metabolise, which would reflect the ideal conditions for velum formation (Legras et al., 2016). 

In conclusion, many genes were significantly differentially expressed between the different stages of velum formation of flor yeast. These genes encompass a broad range of functions, categories, and pathways, and present an extensive opportunity for further research. 

## References

Anders, S., & Huber, W. (2010). Differential expression analysis for sequence count data. Genome Biology, 11(10), R106. https://doi.org/10.1186/gb-2010-11-10-r106

Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. Nature Biotechnology, 34(5), 525–527. https://doi.org/10.1038/nbt.3519

Carlson, M. (2024). _org.Sc.sgd.db: Genome wide annotation for Yeast_. R package version 3.20.0.

Corchete, L. A., Rojas, E. A., Alonso-López, D., De Las Rivas, J., Gutiérrez, N. C., & Burguillo, F. J. (2020). Systematic comparison and assessment of RNA-seq procedures for gene expression quantitative analysis. Scientific Reports, 10(1), 19737. https://doi.org/10.1038/s41598-020-76881-x

De Assis, L. J., Zingali, R. B., Masuda, C. A., Rodrigues, S. P., & Montero-Lomelí, M. (2013). Pyruvate decarboxylase activity is regulated by the Ser/Thr protein phosphatase Sit4p in the yeast Saccharomyces cerevisiae. FEMS Yeast Research, 13(6), 518–528. https://doi.org/10.1111/1567-1364.12052

Flikweert, M. T., Van Der Zanden, L., Janssen, W. M. Th. M., Yde Steensma, H., Van Dijken, J. P., & Pronk, J. T. (1996). Pyruvate decarboxylase: An indispensable enzyme for growth of Saccharomyces cerevisiae on glucose. Yeast, 12(3), 247–257. https://doi.org/10.1002/(SICI)1097-0061(19960315)12:3%3C247::AID-YEA911%3E3.0.CO;2-I

Geistlinger, L., Csaba, G., Santarelli, M., Ramos, M., Schiffer, L., Turaga, N., Law, C., Davis, S., Carey, V., Morgan, M., Zimmer, R., & Waldron, L. (2021). Toward a gold standard for benchmarking gene set enrichment analysis. Briefings in Bioinformatics, 22(1), 545–556. https://doi.org/10.1093/bib/bbz158

Goeman, J. J., & Bühlmann, P. (2007). Analyzing gene expression data in terms of gene sets: Methodological issues. Bioinformatics, 23(8), 980–987. https://doi.org/10.1093/bioinformatics/btm051

Gu, Z., & Hübschmann, D. (2023). simplifyenrichment: A bioconductor package for clustering and visualizing functional enrichment results. Genomics, Proteomics & Bioinformatics, 21(1), 190–202. https://doi.org/10.1016/j.gpb.2022.04.008

Khalfan, M. (2019). Gene Set Enrichment Analysis and Over Representation Analysis analysis using R . https://github.com/gencorefacility/r-notebooks/tree/master

Legras, J.-L., Moreno-Garcia, J., Zara, S., Zara, G., Garcia-Martinez, T., Mauricio, J. C., Mannazzu, I., Coi, A. L., Bou Zeidan, M., Dequin, S., Moreno, J., & Budroni, M. (2016). Flor yeast: New perspectives beyond wine aging. Frontiers in Microbiology, 7. https://doi.org/10.3389/fmicb.2016.00503

Liao, Y., & Shi, W. (2020). Read trimming is not required for mapping and quantification of RNA-seq reads at the gene level. NAR Genomics and Bioinformatics, 2(3), lqaa068. https://doi.org/10.1093/nargab/lqaa068

Love, M. I., Hogenesch, J. B., & Irizarry, R. A. (2016). Modeling of RNA-seq fragment sequence bias reduces systematic errors in transcript abundance estimation. Nature Biotechnology, 34(12), 1287–1291. https://doi.org/10.1038/nbt.3682

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8

Love, M. I., Anders, S., & Huber, W. (2025). Analyzing RNA-seq data with DESeq2. https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#references

Mardanov, A. V., Eldarov, M. A., Beletsky, A. V., Tanashchuk, T. N., Kishkovskaya, S. A., & Ravin, N. V. (2020). Transcriptome profile of yeast strain used for biological wine aging revealed dynamic changes of gene expression in course of flor development. Frontiers in Microbiology, 11. https://doi.org/10.3389/fmicb.2020.00538

Mardanov, A. V., Gruzdev, E. V., Beletsky, A. V., Ivanova, E. V., Shalamitskiy, M. Yu., Tanashchuk, T. N., & Ravin, N. V. (2023). Microbial communities of flor velums and the genetic stability of flor yeasts used for a long time for the industrial production of sherry-like wines. Fermentation, 9(4), 367. https://doi.org/10.3390/fermentation9040367

Maršíková, J., Wilkinson, D., Hlaváček, O., Gilfillan, G. D., Mizeranschi, A., Hughes, T., Begany, M., Rešetárová, S., Váchová, L., & Palková, Z. (2017). Metabolic differentiation of surface and invasive cells of yeast colony biofilms revealed by gene expression profiling. BMC Genomics, 18(1), 814. https://doi.org/10.1186/s12864-017-4214-4

Pagès, H., Carlson, M., Aboyoun, P., Falcon, S., Morgan, M. (2025a). txdbmaker: Tools for making TxDb objects from genomic annotations. doi:10.18129/B9.bioc.txdbmaker, R package version 1.6.2, https://bioconductor.org/packages/txdbmaker.

Pagès, H., Carlson, M., Falcon, S., Li, N. (2025b). AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor. doi:10.18129/B9.bioc.AnnotationDbi, R package version 1.72.0, https://bioconductor.org/packages/AnnotationDbi.

Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods, 14(4), 417–419. https://doi.org/10.1038/nmeth.4197

Patro, R., Mount, S. M., & Kingsford, C. (2014). Sailfish enables alignment-free isoform quantification from RNA-seq reads using lightweight algorithms. Nature Biotechnology, 32(5), 462–464. https://doi.org/10.1038/nbt.2862

Posit team (2025). RStudio: Integrated Development Environment for R. Posit Software, PBC, Boston, MA. http://www.posit.co/.

Rapaport, F., Khanin, R., Liang, Y., Pirun, M., Krek, A., Zumbo, P., Mason, C. E., Socci, N. D., & Betel, D. (2013). Comprehensive evaluation of differential gene expression analysis methods for RNA-seq data. Genome Biology, 14(9), 3158. https://doi.org/10.1186/gb-2013-14-9-r95

Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). Limma powers differential expression analyses for rna-sequencing and microarray studies. Nucleic Acids Research, 43(7), e47. https://doi.org/10.1093/nar/gkv007

Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). Edger: A bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26(1), 139–140. https://doi.org/10.1093/bioinformatics/btp616

Rosati, D., Palmieri, M., Brunelli, G., Morrione, A., Iannelli, F., Frullanti, E., & Giordano, A. (2024a). Differential gene expression analysis pipelines and bioinformatic tools for the identification of specific biomarkers: A review. Computational and Structural Biotechnology Journal, 23, 1154–1168. https://doi.org/10.1016/j.csbj.2024.02.018

Schaaff‐Gerstenschläger, I., Mannhaupt, G., Vetter, I., Zimmermann, F. K., & Feldmann, H. (1993). tkl2 , a second transketolase gene of saccharomyces cerevisiae: Cloning, sequence and deletion analysis of the gene. European Journal of Biochemistry, 217(1), 487–492. https://doi.org/10.1111/j.1432-1033.1993.tb18268.x

Seyednasrollah, F., Laiho, A., & Elo, L. L. (2015). Comparison of software packages for detecting differential expression in RNA-seq studies. Briefings in Bioinformatics, 16(1), 59–70. https://doi.org/10.1093/bib/bbt086

Soneson, C., Love, M. I., & Robinson, M. D. (2016). Differential analyses for RNA-seq: Transcript-level estimates improve gene-level inferences (4:1521). F1000Research. https://doi.org/10.12688/f1000research.7563.2

Taylor, A., Macaulay, V. M., Miossec, M. J., Maurya, A. K., & Buffa, F. M. (2025). GeneFEAST: The pivotal, gene-centric step in functional enrichment analysis interpretation. Bioinformatics, 41(3), btaf100. https://doi.org/10.1093/bioinformatics/btaf100

Wu, T., Hu, E., Xu, S., Chen, M., Guo, P., Dai, Z., Feng, T., Zhou, L., Tang, W., Zhan, L., Fu, X., Liu, S., Bo, X., & Yu, G. (2021). Clusterprofiler 4. 0: A universal enrichment tool for interpreting omics data. The Innovation, 2(3), 100141. https://doi.org/10.1016/j.xinn.2021.100141

Yu, G., Wang, L.-G., Han, Y., & He, Q.-Y. (2012). Clusterprofiler: An r package for comparing biological themes among gene clusters. OMICS : A Journal of Integrative Biology, 16(5), 284–287. https://doi.org/10.1089/omi.2011.0118

Yu, G., Wang, L.-G., Yan, G.-R., & He, Q.-Y. (2015). DOSE: An R/Bioconductor package for disease ontology semantic and enrichment analysis. Bioinformatics, 31(4), 608–609. https://doi.org/10.1093/bioinformatics/btu684

Yu, G. (2025). Enrichplot [Computer software]. Bioconductor. https://doi.org/10.18129/B9.BIOC.ENRICHPLOT

Zara, G., Zara, S., Pinna, C., Marceddu, S., & Budroni, M. (2009). FLO11 gene length and transcriptional level affect biofilm-forming ability of wild flor strains of Saccharomyces cerevisiae. Microbiology, 155(12), 3838–3846. https://doi.org/10.1099/mic.0.028738-0
