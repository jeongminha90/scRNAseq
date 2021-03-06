# Trajectory and signature analyses for iPSC/iNSC reprogramming intermdeiate cells<br><br><sub>Title: "Analyzing intermediate populations during OSKM-mediated reprogramming"</sub>


### Preprocessing.Rmd
* Rmarkdown for loading raw dataset & quality control
* Download dataset used in the analysis
  * [`raw.zip`](https://figshare.com/s/ecf794cfe2776980f4de): raw 10xGenomics output data
  * [`aggregation_csv.csv`](https://github.com/jeongminha90/scRNAseq/blob/main/aggregation_csv.csv): to assign sample group names to indexes


### Extended_data_Fig7_9_whole_cells.Rmd
* Rmarkdown for generation of UMAPs for whole cells and masking the signatures in Extended data figure 7 & 8
* [`dIC_seurat3_whole.rds`](https://figshare.com/s/d083a7d6649a6e83f875): Seurat object of whole cells

### Fig3_D5-D7_slingshot.Rmd
* Rmarkdown for preprocessing and trajectory of D5~D7 iPSC/iNSC reprogramming intermediate cells

### wot_command_line.sh
* Rmarkdown for wot analysis of D5~D7 iPSC/iNSC reprogramming intermediate cells
