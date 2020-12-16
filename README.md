# Trajectory and signature analysis for iPSC/iNSC reprogramming intermdeiate cells

__Title: "Analyzing intermediate populations during OSKM-mediated reprogramming"__


### Preprocessing.Rmd
* Rmarkdown for loading raw dataset & quality control
* `Seurat v2.3.x` is required.
* [View & download the code](https://github.com/jeongminha90/scRNAseq/blob/main/Preprocessing.Rmd)
* Download the dataset used in the analysis
[raw 10xGenomics output data](https://figshare.com/s/ecf794cfe2776980f4de)
[Aggregation](https://github.com/jeongminha90/scRNAseq/blob/main/aggregation_csv.csv) file to assign sample group names to indexes



### Extended_data_Fig.7_8.Rmd
* Rmarkdown for Extended data figure 7 & 8
* [View & download the code](https://github.com/jeongminha90/scRNAseq/blob/main/Extended%20Data%20Fig.7%2C8.Rmd)
* [Click here to download the dataset used in the analysis](https://figshare.com/s/2d5e45d42f50dc3c6d9c)
* [Phenotyping criteria](https://github.com/jeongminha90/scRNAseq/blob/main/Phenotyping%20Criteria.csv) for calculation of signature scores



### Fig3_D5-D7_tranjectory_signature.Rmd
* Rmarkdown for trajectory and cell signature analyses of D5~D7 iPSC/iNSC reprogramming intermediate cells
* [View & download the code](https://github.com/jeongminha90/scRNAseq/blob/main/Fig3_D5-D7_trajectory_signature.Rmd)
* [Click here to download the dataset used in the analysis](https://figshare.com/articles/dataset/D5-D7_mipsc_normalized_scaled/13383191)
* Precalcualted [`gene set scores`](https://figshare.com/articles/dataset/gene_set_scores_csv/13383212) for plotting the cell signatures by [`WOT: 'gene set scores' function`](https://broadinstitute.github.io/wot/).
