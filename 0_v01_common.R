#### Default packages and shortcuts ####

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(gridExtra)
library(reshape2) #melt()
library(stringr)
library(RColorBrewer)
library(extrafont)
library(vegan)
library(colorRamps)
library(gdata) # for NAToUnknown()
library(phyloseq)
library(DESeq2)
library(ggpubr) # ggboxplot
library(microbiome)
library(agricolae) # HSD.test
library(eulerr) # venn()
library(ggtern)

#### Helper functions ####

'%!in%' <- function(x,y)!('%in%'(x,y))

PlotOrd <- function(data2plot, components = c("PC1", "PC2"), groups = groups, treatment = treatment, variance = TRUE, connect = FALSE, path = groups, species = NULL, ...){
  ## plot ordinations as a biplot 
  # Arguments: 
  # data - a data frame containing sample and/or species scores from an ordination analysis
  # species - an optional species scores to be plotted as vectors or points
  # variance <-   print variance explained to the axis titles
  # axes <- choose axes to plot
  #   two.groups <- c("#006837", "#1A9850", "#A6D96A", "#D9EF8B", "#A50026", "#D73027", "#FDAE61", "#FEE08B") # based on brewer: RdYlGn
  p <- ggplot(data2plot, aes_string(x = components[1], y = components[2])) +
    geom_point(aes_string(colour = groups, shape = treatment), size = 5, alpha = 3/4) +
    theme(panel.background = element_rect(fill = "#F2F2F2"), 
          panel.grid.minor = element_blank()) + 
    guides(colour = guide_legend(title = ""), shape = guide_legend(title = "")) +
    if (connect == TRUE) {
      p <- p + geom_path(aes_string(group = path, colour = groups), alpha = 1/2) + guides(colour = guide_legend(title = "Gradient"))
    }
  if (!is.null(species)) {
    p + geom_path(data2plot = arrws, aes(PC1, PC2, group = species), 
                  colour = "red", arrow = arrow(length = unit(0.05, "npc"))) + 
      geom_text(data = arrws[5:8, ], aes(label = species), hjust = 0.5, vjust = -1, colour = "black")
  }
  return(p)
}

PlotOrd.single <- function(data2plot, components = c("PC1", "PC2"), groups = groups, variance = TRUE, connect = FALSE, path = groups, species = NULL, ...){
  ## plot ordinations as a biplot 
  # Arguments: 
  # data - a data frame containing sample and/or species scores from an ordination analysis
  # species - an optional species scores to be plotted as vectors or points
  # variance <-   print variance explained to the axis titles
  # axes <- choose axes to plot
  #   two.groups <- c("#006837", "#1A9850", "#A6D96A", "#D9EF8B", "#A50026", "#D73027", "#FDAE61", "#FEE08B") # based on brewer: RdYlGn
  p <- ggplot(data2plot, aes_string(x = components[1], y = components[2])) +
    geom_point(aes_string(colour = groups), size = 5, alpha = 3/4) +
    scale_shape_manual(values=c(18, 17, 15, 16)) + # 7
    theme(panel.background = element_rect(fill = "#F2F2F2"), 
          panel.grid.minor = element_blank()) + 
    guides(colour = guide_legend(title = ""), shape = guide_legend(title = "")) +
    if (connect == TRUE) {
      p <- p + geom_path(aes_string(group = path, colour = groups), alpha = 1/2) + guides(colour = guide_legend(title = "Gradient"))
    }
  if (!is.null(species)) {
    p + geom_path(data2plot = arrws, aes(PC1, PC2, group = species), 
                  colour = "red", arrow = arrow(length = unit(0.05, "npc"))) + 
      geom_text(data = arrws[5:8, ], aes(label = species), hjust = 0.5, vjust = -1, colour = "black")
  }
  return(p)
}

PlotOrd.single.bP <- function(data2plot, components = c("PC1", "PC2"), groups = groups, variance = TRUE, connect = FALSE, path = groups, species = NULL, ...){
  ## plot ordinations as a biplot 
  # Arguments: 
  # data - a data frame containing sample and/or species scores from an ordination analysis
  # species - an optional species scores to be plotted as vectors or points
  # variance <-   print variance explained to the axis titles
  # axes <- choose axes to plot
  #   two.groups <- c("#006837", "#1A9850", "#A6D96A", "#D9EF8B", "#A50026", "#D73027", "#FDAE61", "#FEE08B") # based on brewer: RdYlGn
  p <- ggplot(data2plot, aes_string(x = components[1], y = components[2])) +
    geom_point(aes_string(colour = groups), size = 5, alpha = 3/4) +
    scale_shape_manual(values=c(18, 17, 15, 16)) + # 7
    theme(panel.background = element_rect(fill = "#F2F2F2"), 
          panel.grid.minor = element_blank()) + 
    guides(colour = guide_legend(title = ""), shape = guide_legend(title = "")) +
    scale_colour_manual(values = c(colPalette03)) +
    if (connect == TRUE) {
      p <- p + geom_path(aes_string(group = path, colour = groups), alpha = 1/2) + guides(colour = guide_legend(title = "Gradient"))
    }
  if (!is.null(species)) {
    p + geom_path(data2plot = arrws, aes(PC1, PC2, group = species), 
                  colour = "red", arrow = arrow(length = unit(0.05, "npc"))) + 
      geom_text(data = arrws[5:8, ], aes(label = species), hjust = 0.5, vjust = -1, colour = "black")
  }
  return(p)
}

TaxThresh <- function(Taxonomy, tax.level = "Order", thresh = 1){
  ## clump together all taxa with frequency below threshold as "rare" and all those with bs values below threshold as "unclassified"
  # vars: Taxonomy=Taxonomy df, bs.vals=bs values df, bs.thresh=threshold of bs value to drop, tax.level=which level should be considered, thresh=threshold for "rare"
  # make all below bs threshold unclassified
  Taxonomy[, tax.level] <- factor(Taxonomy[, tax.level], levels = c(levels(Taxonomy[, tax.level]), "z.Rare")) # add "Rare" and "Unclassified" levels
  
  # count frequencies of taxonomies in samples (unique sequences only * their counts = all sequences)
  count.tax <- plyr::count(df = Taxonomy, vars = tax.level, wt_var = "Freq") # data.frame, Vars, weight
  count.tax[, 2] <- (count.tax[, 2] / sum((count.tax[, 2]))) * 100 # norm to 100%
  
  ab.tax <- droplevels(count.tax[count.tax$freq > thresh, ]) # remove all taxa below threshold (keep abundent taxa)
  
  # generate logical vector of abundant taxa
  inds <- vector(mode = "logical", length = nrow(Taxonomy)) # assign a logical vector length Taxonomy
  # cycle through taxa and generate logical vector of abundant taxa
  for (i in 1:length(levels(ab.tax[, 1]))) inds <- as.logical(inds + grepl(levels(ab.tax[,1])[i], Taxonomy[, tax.level]))
  possibleError <- tryCatch(Taxonomy[!inds, tax.level] <- "z.Rare", error = function(e) e) # replace all rare taxa with "Rare"
  ab.Taxonomy <- droplevels(Taxonomy) # remove factor levels reclassified as either "rare" or "unclassified"
  return(ab.Taxonomy)
}

Palette <- RColorBrewer::brewer.pal(n = 11, "RdYlBu")[c(1,3,10)]

# The palette with black:
barPalette <- c("#000000", "#FFCCCC", "#99FF33", "#6600FF", "#FF0099", "#CCFFFF", "#FF6600", "#003399", "#009900", "#FFFF33", "#999966", "#000033", "#CCCCCC", "#FF9933", "#00CCCC", "#990033", "#666666")

prim.col.1 <- gsub("0","1", colorRamps::primary.colors(n=26))[c(-8,-16,-22,-25)]
prim.col.4 <- gsub("0","4", colorRamps::primary.colors(n=26))[c(-8,-16,-22,-25)]
colPalette01 <- c(prim.col.1,prim.col.4)
colPalette03 <- c("#111111", "#811111", "#FF1111", "#118111", "#FF8111", "#818111", "#11FF11", "#111181", "#FFFF11", "#811181", 
                  "#FF1181", "#118181", "#FF8181", "#818181", "#81FF81", "#1111FF", "#FFFF81", "#FF11FF", "#8111FF", "#FF81FF", "#8181FF", "#81FFFF", 
                  "#444444", "#844444", "#FF4444", "#448444", "#FF8444", "#848444", "#44FF44", "#444484", "#FFFF44", "#844484",
                  "#FF4484", "#448484", "#FF8484", "#848484", "#84FF84", "#4444FF", "#FFFF84", "#FF44FF", "#8444FF", "#FF84FF", "#8484FF", "#84FFFF",
                  "#666666", "#866666", "#FF6666", "#668666", "#FF8666", "#868666", "#66FF66", "#666686", "#FFFF66", "#866686",
                  "#FF6686", "#668686", "#FF8686", "#868686", "#86FF86", "#6666FF", "#FFFF86", "#FF66FF", "#8666FF", "#FF86FF", "#8686FF", "#86FFFF",
                  "#999999", "#899999", "#FF9999", "#998999", "#FF8999", "#898999", "#99FF99", "#999989", "#FFFF99", "#899989",
                  "#FF9989", "#998989", "#FF8989", "#898989", "#89FF89", "#9999FF", "#FFFF89", "#FF99FF", "#8999FF", "#FF89FF", "#8989FF", "#89FFFF")

#### Project settings ####

RESULTS_DIR <- "Results"
DATA_DIR <-
  RESULTS_DIR %>%
  list.dirs(recursive = FALSE) %>%
  head(1) # tail(1) if the target folder is the last one
RD_DIR <- "RD"
PLOTS_DIR <- "Plots"

PROJECT_NAME <-
  DATA_DIR %>%
  basename()

cat(
  crayon::blue(crayon::bold("Common settings:")), "\n",
  crayon::blue("PROJECT_NAME: "), PROJECT_NAME, "\n",
  crayon::blue("DATA_DIR: "), DATA_DIR, "\n",
  crayon::blue("RD_DIR: "), RD_DIR, "\n",
  crayon::blue("RESULTS_DIR: "), RESULTS_DIR, "\n",
  crayon::blue("PLOTS_DIR: "), PLOTS_DIR, "\n",
  sep = ""
)


#### Default options ####

# setting the seed makes all output reproducible:
set.seed(42)


#### I/O ####

# make sure all folders exist:
dir.create(RD_DIR, showWarnings = FALSE)
dir.create(RESULTS_DIR, showWarnings = FALSE)
dir.create(PLOTS_DIR, showWarnings = FALSE)
