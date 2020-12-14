
# R code for iPSC/iNSC reprogramming  cells at D5~D7
# data are already normalized/scaled
# UMAP/clusters included

library(Seurat)
library(slingshot)
library(tidyverse)
library(patchwork)
library(splitstackshape)

setwd("dirctory/containing/mipsc.RData/")

# load seurat3 object 
load(file="mipsc_D5-D7_Seurat3_normalized_scaled.RData")

# plot UMAP
DimPlot(mipsc, reduction = "umap", 
        group.by = "conds", pt.size = 0.2, 
        label = F, repel = TRUE,
        order =c("D6_iNSC","D6_iPSC","D5","D7_iNSC","D7_iPSC")) +
  coord_fixed() +
  scale_color_manual(values = c("pink","turquoise","darkgreen","red","blue"))

DimPlot(mipsc, reduction = "umap", 
        group.by = "conds", pt.size = 0.2, 
        label = F, repel = TRUE, 
        split.by="conds", ncol = 3,
        order =c("D7_iNSC","D7_iPSC","D6_iNSC","D6_iPSC","D5")) +
  coord_fixed() +
  scale_color_manual(values = c("darkgreen","red","blue","pink","turquoise"))

# plot Clusters
library(patchwork)
p1<- DimPlot(mipsc, reduction = "umap",label = T, size=.7, label.size = 5) +
  coord_fixed() +  ggtitle("clustering: D5 - D7")

dta <- mipsc@meta.data %>% group_by(conds, seurat_clusters)
total <- mipsc@meta.data %>% count(seurat_clusters)

p2<- ggplot(dta, aes(x=seurat_clusters)) + theme_classic() +
  geom_bar(aes(fill = conds)) + #position="fill") +
  geom_text(data=total, aes(x=seurat_clusters, y=2600, label=as.factor(n)), size = 4, angle=90, hjust=1) +
  labs(x="Clusters", y = "Cell counts") +
  scale_color_manual(values = c("pink","turquoise","darkgreen","red","blue"))
p3<- ggplot(dta, aes(x=seurat_clusters)) + theme_classic() +
  geom_bar(position="fill", aes(fill=conds)) +
  geom_text(data=total, aes(x=seurat_clusters, y=1.15, label=as.factor(n)), size = 4, angle=90, hjust=1) +
  scale_y_continuous(labels = scales::percent) +
  labs(x="Clusters", y = "Fraction") +
  scale_color_manual(values = c("pink","turquoise","darkgreen","red","blue"))

p1 + (p2/p3) + plot_layout(ncol = 2, widths = c(5,5))

# Signature scores
# open gene set signature scores from WOT 'gene_set_scores'
gs <- read.csv("gene_set_scores.csv", header = T, row.names = 1)
embeddings <- Embeddings(mipsc, reduction = "umap")
df <- merge(data.frame(mipsc$conds), gs, by=0)
df2 <- merge(df, embeddings, by.x="Row.names", by.y=0)
head(df2)

data.frame(df$Row.names) %>% 
  cSplit("df.Row.names", "-") %>% 
  dplyr::count(df.Row.names_1)

theme_set(theme_classic())
p0<- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, color=Pou5f1)) + 
  ggtitle("Pou5f1") +
  geom_point(size=.2) +
  coord_fixed() +
  scale_color_gradient2(low="darkblue", mid="lightgrey", high="red", midpoint = -.5) + 
  theme(legend.title = element_blank())
p1<- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, color=Pluripotency)) + 
  ggtitle("Pluripotency") +
  geom_point(size=.2) +
  coord_fixed() +
  scale_color_gradient2(low="darkblue", mid="lightgrey", high="red") + 
  theme(legend.title = element_blank())
p2<- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, color=Neural.identity)) + 
  ggtitle("Neural.identity") +
  geom_point(size=.2) +
  coord_fixed() +
  scale_color_gradient2(low="darkblue", mid="lightgrey", high="red") + 
  theme(legend.title = element_blank())
p3<- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, color=MEF.identity)) + 
  ggtitle("MEF.identity") +
  geom_point(size=.2) +
  coord_fixed() +
  scale_color_gradient2(low="darkblue", mid="lightgrey", high="red") + 
  theme(legend.title = element_blank())
p4<- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, color=Epithelial.identity)) + 
  ggtitle("Epithelial.identity") +
  geom_point(size=.2) +
  coord_fixed() +
  scale_color_gradient2(low="darkblue", mid="lightgrey", high="red") + 
  theme(legend.title = element_blank())
p5<- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, color=Shisa8)) + 
  ggtitle("Shisa8") +
  geom_point(size=.2) +
  coord_fixed() +
  scale_color_gradient2(low="darkblue", mid="lightgrey", high="red") + 
  theme(legend.title = element_blank())
p6<- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, color=Dsp)) + 
  ggtitle("Dsp") +
  geom_point(size=.2) +
  coord_fixed() +
  scale_color_gradient2(low="darkblue", mid="lightgrey", high="red") + 
  theme(legend.title = element_blank())
p7<- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, color=Cell.cycle)) + 
  ggtitle("Cell.cycle") +
  geom_point(size=.2) +
  coord_fixed() +
  scale_color_gradient2(low="darkblue", mid="lightgrey", high="red") + 
  theme(legend.title = element_blank())
p8<- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, color=SASP)) + 
  ggtitle("SASP") +
  geom_point(size=.2) +
  coord_fixed() +
  scale_color_gradient2(low="darkblue", mid="lightgrey", high="red") + 
  theme(legend.title = element_blank())

(p1|p2|p3|p4) / (p5|p6|p7|p8)

########################################################################################
# Slingshot
library(slingshot)
library(scales)
library(viridis)
library(shape)

###################################
# Functions for color/alpha control

# # FUNCTION 1: cell_pal
#' Assign a color to each cell based on some value
#' @param cell_vars Vector indicating the value of a variable associated with cells.
#' @param pal_fun Palette function that returns a vector of hex colors, whose
#' argument is the length of such a vector.
#' @param ... Extra arguments for pal_fun.
#' @return A vector of hex colors with one entry for each cell.
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

# # FUNCTION 2: add alpha 
addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

# # Function3: colorRampPaletteAlpha()
colorRampPaletteAlpha <- function(colors, n=32, interpolate='linear') {
  # Create the color ramp normally
  cr <- colorRampPalette(colors, interpolate=interpolate)(n)
  # Find the alpha channel
  a <- col2rgb(colors, alpha=T)[4,]
  # Interpolate
  if (interpolate=='linear') {
    l <- approx(a, n=n)
  } else {
    l <- spline(a, n=n)
  }
  l$y[l$y > 255] <- 255 # Clamp if spline is > 255
  cr <- addalpha(cr, l$y/255.0)
  return(cr)
}

########################################

cl <- mipsc$seurat_clusters

# run Slingshot
sds <- slingshot(Embeddings(mipsc, "umap"), 
                 clusterLabels = cl,
                 start.clus = 0,
                 end.clus = NULL,
                 stretch = 0,
                 allow.breaks = FALSE)

par(mfrow = c(1, 2))
cell_colors <- cell_pal(cl, hue_pal())
plot(reducedDim(sds), col = addalpha(cell_colors,0.7), pch = 19, cex = .3, asp=1)
lines(sds, lwd = 2, type = 'lineages', col = 'black')
plot(reducedDim(sds), col = addalpha(cell_colors,0.7), pch = 19, cex = .3, asp=1)
lines(sds, lwd = 3, col = 'black', lty = 3)
#

