# Trajectory and signature analyses for iPSC/iNSC reprogramming intermdeiate cells

##Title: "Analyzing intermediate populations during OSKM-mediated reprogramming"


### Preprocessing.Rmd
* Rmarkdown for loading raw dataset & quality control
* `Seurat v2.3.x` is required.
* [View & download the code](https://github.com/jeongminha90/scRNAseq/blob/main/Preprocessing.Rmd)
* Download dataset used in the analysis
  * [`raw.zip`](https://figshare.com/s/ecf794cfe2776980f4de): raw 10xGenomics output data
  * [`Aggregation`](https://github.com/jeongminha90/scRNAseq/blob/main/aggregation_csv.csv): to assign sample group names to indexes



### Extended_data_Fig.7_8.Rmd
* Rmarkdown for generation of UMAPs for whole cells and masking the signatures in Extended data figure 7 & 8
* [View & download the code](https://github.com/jeongminha90/scRNAseq/blob/main/Extended%20Data%20Fig.7%2C8.Rmd)
* Download datasets used in the analysis
 * [`miPSC_seurat.RData`](https://figshare.com/s/2d5e45d42f50dc3c6d9c): Seurat object of whole cells
 * [`Phenotyping criteria`](https://github.com/jeongminha90/scRNAseq/blob/main/Phenotyping%20Criteria.csv): lists of gene sets used in the calculation of signature scores



### Fig3_D5-D7_tranjectory_signature.Rmd
* Rmarkdown for trajectory and cell signature analyses of D5~D7 iPSC/iNSC reprogramming intermediate cells
* [View & download the code](https://github.com/jeongminha90/scRNAseq/blob/main/Fig3_D5-D7_trajectory_signature.Rmd)
* Download datasets used in the analysis
 * [`mipsc_D5-D7_Seurat3_normalized_scaled.RData`](https://figshare.com/articles/dataset/D5-D7_mipsc_normalized_scaled/13383191): Seurat object of cells in D5~D7
 * [`gene set scores`](https://figshare.com/articles/dataset/gene_set_scores_csv/13383212): for plotting the cell signatures by [`WOT 'gene set scores' function`](https://broadinstitute.github.io/wot/cli_documentation/).
