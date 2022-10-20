#### Default packages and shortcuts ####

# library(gdata)
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
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("phyloseq")
# BiocManager::install("DESeq2")
# biocLite('phyloseq')
library(DESeq2)
library(ggpubr) # ggboxplot
library(microbiome)
# devtools::install_github("gauravsk/ranacapa")
# library(ranacapa) #rarefaction
# library(car) # leveneTest()
library(agricolae) # HSD.test
library(eulerr) # venn()
# library(ggtern)

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
    #     scale_size_discrete(range = c(8,8)) +
    theme(panel.background = element_rect(fill = "#F2F2F2"), 
          panel.grid.minor = element_blank()) + 
    guides(colour = guide_legend(title = ""), shape = guide_legend(title = "")) +
    #     scale_colour_manual(values = NevadaRainbow) +
    #     scale_colour_brewer(palette = "Dark2") # w/o this command -> unlimited standard colours
    #      scale_colour_brewer(palette = "Paired") # testing!
    #     theme(legend.position = c(0.9, 0.76), legend.background = element_rect(fill = "white", colour = NA))
    if (connect == TRUE) {
      p <- p + geom_path(aes_string(group = path, colour = groups), alpha = 1/2) + guides(colour = guide_legend(title = "Gradient"))
    }
  if (!is.null(species)) {
    p + geom_path(data2plot = arrws, aes(PC1, PC2, group = species), 
                  colour = "red", arrow = arrow(length = unit(0.05, "npc"))) + 
      # scale_colour_gradient(limits = c(1,4) legend = F) +
      geom_text(data = arrws[5:8, ], aes(label = species), hjust = 0.5, vjust = -1, colour = "black")
  }
  return(p)
}

PlotOrd.bP <- function(data2plot, components = c("PC1", "PC2"), groups = groups, treatment = treatment, variance = TRUE, connect = FALSE, path = groups, species = NULL, ...){
  ## plot ordinations as a biplot 
  # Arguments: 
  # data - a data frame containing sample and/or species scores from an ordination analysis
  # species - an optional species scores to be plotted as vectors or points
  # variance <-   print variance explained to the axis titles
  # axes <- choose axes to plot
  #   two.groups <- c("#006837", "#1A9850", "#A6D96A", "#D9EF8B", "#A50026", "#D73027", "#FDAE61", "#FEE08B") # based on brewer: RdYlGn
  p <- ggplot(data2plot, aes_string(x = components[1], y = components[2])) +
    geom_point(aes_string(colour = groups, shape = treatment), size = 5, alpha = 3/4) +
    #     scale_size_discrete(range = c(8,8)) +
    theme(panel.background = element_rect(fill = "#F2F2F2"), 
          panel.grid.minor = element_blank()) + 
    guides(colour = guide_legend(title = ""), shape = guide_legend(title = "")) +
    scale_colour_manual(values = colPalette03) +
    #     scale_colour_brewer(palette = "Dark2") # w/o this command -> unlimited standard colours
    #      scale_colour_brewer(palette = "Paired") # testing!
    #     theme(legend.position = c(0.9, 0.76), legend.background = element_rect(fill = "white", colour = NA))
    if (connect == TRUE) {
      p <- p + geom_path(aes_string(group = path, colour = groups), alpha = 1/2) + guides(colour = guide_legend(title = "Gradient"))
    }
  if (!is.null(species)) {
    p + geom_path(data2plot = arrws, aes(PC1, PC2, group = species), 
                  colour = "red", arrow = arrow(length = unit(0.05, "npc"))) + 
      # scale_colour_gradient(limits = c(1,4) legend = F) +
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
    # scale_size_discrete(range = c(8,8)) +
    theme(panel.background = element_rect(fill = "#F2F2F2"), 
          panel.grid.minor = element_blank()) + 
    guides(colour = guide_legend(title = ""), shape = guide_legend(title = "")) +
    #  stat_ellipse(type = "t") +
    # scale_colour_manual(values = c(colorRamps::primary.colors(n=25))) +
    # scale_colour_brewer(palette = "Dark2") # w/o this command -> unlimited standard colours
    #      scale_colour_brewer(palette = "Paired") # testing!
    #     theme(legend.position = c(0.9, 0.76), legend.background = element_rect(fill = "white", colour = NA))
    if (connect == TRUE) {
      p <- p + geom_path(aes_string(group = path, colour = groups), alpha = 1/2) + guides(colour = guide_legend(title = "Gradient"))
    }
  if (!is.null(species)) {
    p + geom_path(data2plot = arrws, aes(PC1, PC2, group = species), 
                  colour = "red", arrow = arrow(length = unit(0.05, "npc"))) + 
      # scale_colour_gradient(limits = c(1,4) legend = F) +
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
    # scale_size_discrete(range = c(8,8)) +
    theme(panel.background = element_rect(fill = "#F2F2F2"), 
          panel.grid.minor = element_blank()) + 
    guides(colour = guide_legend(title = ""), shape = guide_legend(title = "")) +
    #  stat_ellipse(type = "t") +
    scale_colour_manual(values = c(colPalette03)) +
    # scale_colour_brewer(palette = "Dark2") # w/o this command -> unlimited standard colours
    #      scale_colour_brewer(palette = "Paired") # testing!
    #     theme(legend.position = c(0.9, 0.76), legend.background = element_rect(fill = "white", colour = NA))
    if (connect == TRUE) {
      p <- p + geom_path(aes_string(group = path, colour = groups), alpha = 1/2) + guides(colour = guide_legend(title = "Gradient"))
    }
  if (!is.null(species)) {
    p + geom_path(data2plot = arrws, aes(PC1, PC2, group = species), 
                  colour = "red", arrow = arrow(length = unit(0.05, "npc"))) + 
      # scale_colour_gradient(limits = c(1,4) legend = F) +
      geom_text(data = arrws[5:8, ], aes(label = species), hjust = 0.5, vjust = -1, colour = "black")
  }
  return(p)
}

PlotOrd.Fig3a <- function(data2plot, components = c("PC1", "PC2"), groups = groups, treatment = treatment, variance = TRUE, connect = FALSE, path = groups, species = NULL, ...){
  ## plot ordinations as a biplot 
  # Arguments: 
  # data - a data frame containing sample and/or species scores from an ordination analysis
  # species - an optional species scores to be plotted as vectors or points
  # variance <-   print variance explained to the axis titles
  # axes <- choose axes to plot
  #   two.groups <- c("#006837", "#1A9850", "#A6D96A", "#D9EF8B", "#A50026", "#D73027", "#FDAE61", "#FEE08B") # based on brewer: RdYlGn
  p <- ggplot(data2plot, aes_string(x = components[1], y = components[2])) +
    geom_point(aes_string(colour = groups, shape = treatment), size = 5, alpha = 3/4) +
    scale_shape_manual(values=c(17, 16)) + #17
    #     scale_size_discrete(range = c(8,8)) +
    theme(panel.background = element_rect(fill = "#F2F2F2"), 
          panel.grid.minor = element_blank()) + 
    guides(colour = guide_legend(title = ""), shape = guide_legend(title = "")) +
    #     scale_colour_manual(values = NevadaRainbow) +
    #     scale_colour_brewer(palette = "Dark2") # w/o this command -> unlimited standard colours
    #      scale_colour_brewer(palette = "Paired") # testing!
    #     theme(legend.position = c(0.9, 0.76), legend.background = element_rect(fill = "white", colour = NA))
    if (connect == TRUE) {
      p <- p + geom_path(aes_string(group = path, colour = groups), alpha = 1/2) + guides(colour = guide_legend(title = "Gradient"))
    }
  if (!is.null(species)) {
    p + geom_path(data2plot = arrws, aes(PC1, PC2, group = species), 
                  colour = "red", arrow = arrow(length = unit(0.05, "npc"))) + 
      # scale_colour_gradient(limits = c(1,4) legend = F) +
      geom_text(data = arrws[5:8, ], aes(label = species), hjust = 0.5, vjust = -1, colour = "black")
  }
  return(p)
}

# f?r Heatmap! definiert pro OTU+Taxonomie wie h?ufig in welcher einzelnen Probe
PrepOtuTax <- function(OTU.file, ...) {
  headers <- read.table(OTU.file, header = F, nrow = 1)
  headers  <- apply(headers, 1, function(x) gsub("SAMPLE\\.","",x))
  
  read.table(OTU.file, header = F, row.names = 1, skip = 1, colClasses = c("character", rep("integer", 84), rep("NULL", 12))) -> OTU.mat
  colnames(OTU.mat) <- headers[2:85]
  OTU.mat$Total <- rowSums(OTU.mat)
  
  read.table(OTU.file, header = F, row.names = 1, skip = 1, colClasses = c("character", rep("NULL", 84), rep("factor", 12)))[, seq(1, 12, 2)] -> Taxonomy
  colnames(Taxonomy) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
  Taxonomy$OTU <- rownames(Taxonomy)
  
  Taxonomy <- cbind(Taxonomy, OTU.mat) # combine data.frames (OTU/Taxonomy & Sample Abundance; merge doesnt work on such large dbs)
  Taxonomy <- melt(Taxonomy, id = c("OTU", "Total", "Domain", "Phylum", "Class", "Order", "Family", "Genus"), variable.name = 'Name', value.name = 'Freq', na.rm = TRUE) # --> jede OTU in jedem sample einen eintrag mit frequenz
  return(Taxonomy)
}

TaxThresh <- function(Taxonomy, tax.level = "Order", thresh = 1){
  ## clump together all taxa with frequency below threshold as "rare" and all those with bs values below threshold as "unclassified"
  # vars: Taxonomy=Taxonomy df, bs.vals=bs values df, bs.thresh=threshold of bs value to drop, tax.level=which level should be considered, thresh=threshold for "rare"
  # make all below bs threshold unclassified
  Taxonomy[, tax.level] <- factor(Taxonomy[, tax.level], levels = c(levels(Taxonomy[, tax.level]), "z.Rare")) # add "Rare" and "Unclassified" levels
  
  #   Taxonomy[bs.vals[, tax.level] < bs.thresh,tax.level] <- factor("Unclassified") # replace in Taxonomy all rows under tax.level whose corresponding bs.vals are below bs.thresh with "unclassified"
  
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

PlotTaxHt <- function(tax.data2plot, xvar = "Sample", yvar = "Order", colours...){
  ## plot taxa as heat map
  # parameters: data.frame;
  # xvar;
  # yvar - # set tax level as label for plot. The real taxa level is the one found in the column "Taxa"!!
  
  ggplot(tax.data2plot, aes_string(x = xvar, y = "Taxa", fill = "Rel.Ab")) +
    geom_tile() +
    scale_fill_gradientn(colours = colours, guide = "colourbar", values = c(0,.2,1), space = "Lab") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_equal() +
    #     theme_bw() +
    labs(fill = "Seq. abundance (%)") +
    #     coord_flip() +
    ylab(yvar) +
    xlab(xvar) +
    theme(axis.text.x = element_text(hjust = 1.0, angle = 45.0))
}

PlotSummarySingle <- function(data, x = "ElevationBelt", y = "Estimated", ymin = "Lower", ymax = "Upper", colour = "Metric", ...){
  ## plot output from catchall or summary single with error bars
  # parameters: data.frame;
  ggplot(data, aes_string(x = x, y = y, ymin = ymin, ymax = ymax)) +
    geom_point(aes_string(colour = colour), size = 3) +
    geom_errorbar(alpha = 1/2, width = 0.3) +
    xlab("") +
    ylab("") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1)) +
    facet_grid(Metric ~ ., scale = "free")
}

PlotSummarySingle_overlay <- function(data2plot, x = "ElevationBelt", y = "Estimate", ymin = "lerr", ymax = "herr", colour = "Metric", ...){
  ## plot output from catchall or summary single with error bars, overlay several metrics
  # parameters: data.frame;
  ggplot(data2plot, aes_string(x = x, y = y, ymin = ymin, ymax = ymax)) +
    geom_point(data = data2plot[data2plot$Metric != "Observed S", ],
               aes_string(x = x, y = y, ymin = ymin, ymax = ymax, colour = colour), position = position_dodge(0.4), size = 4) +
    geom_errorbar(data = data2plot[data2plot$Metric != "Observed S", ],
                  aes_string(x = x, y = y, ymin = ymin, ymax = ymax, colour = colour), alpha = 1/2, position = position_dodge(0.4), width = 0.3) +
    geom_bar(data = data2plot[data2plot$Metric == "Observed S", ],
             aes_string(x = x, y = y, ymin = ymin, ymax = ymax), stat = "identity", fill = "green4")  +
    geom_errorbar(data = data2plot[data2plot$Metric == "Observed S", ],
                  aes_string(x = x, y = y, ymin = ymin, ymax = ymax), colour = "black", alpha = 1/2, width = 0.3) +
    geom_errorbar(data = data2plot[data2plot$Metric == "Observed S", ],
                  aes_string(x = x, y = y, ymin = ymin, ymax = ymax, colour = colour), alpha = 1/2, width = 0.3) +
    #     guides(colour=guide_legend(c(levels(data2plot$Metric)))) +
    scale_colour_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8"), breaks = c("Observed S", "Ace", "Parametric")) +
    xlab("") +
    ylab("OTUs") +
    theme(axis.text.x = element_text(angle = 45.0, vjust = 0.9, hjust = 1))
}

Palette <- RColorBrewer::brewer.pal(n = 11, "RdYlBu")[c(1,3,10)]

# The palette with black:
barPalette <- c("#000000", "#FFCCCC", "#99FF33", "#6600FF", "#FF0099", "#CCFFFF", "#FF6600", "#003399", "#009900", "#FFFF33", "#999966", "#000033", "#CCCCCC", "#FF9933", "#00CCCC", "#990033", "#666666")
#+ scale_fill_manual(values=barPalette) #add this after the stacked barplot command! 

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
# barPalette2 <- c("#000000", "#666666", "#999966", "#FFCCCC", "#CCFFFF", "#00CCCC")
# barPalette10 <- colorRamps::primary.colors(n=10)
# barPalette15 <- colorRamps::primary.colors(n=15)
# barPalette17 <- colorRamps::primary.colors(n=17)
# barPalette20 <- colorRamps::primary.colors(n=20)
# barPalette22 <- colorRamps::primary.colors(n=22)
# barPalette25 <- colorRamps::primary.colors(n=25)
# barPalette26 <- colorRamps::primary.colors(n=26)
# barPalette29 <- colorRamps::primary.colors(n=29)
# barPalette31 <- colorRamps::primary.colors(n=31)
# barPalette35 <- colorRamps::primary.colors(n=35)
# barPalette41 <- colorRamps::primary.colors(n=41)
# barPalette53 <- colorRamps::primary.colors(n=53)
# barPalette71 <- colorRamps::primary.colors(n=71)
# barPalette80 <- colorRamps::primary.colors(n=80)
# colours <- c("#FFFFFF",  "#fd8d3c", "#feb24c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026", "#5e011d", "#330010", "#000000")

# write_tsv <- function(x, file, ...) data.table::fwrite(x = x, file = file, sep = "\t", ...)

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
  # crayon::blue("MIN_READS: "), MIN_READS, "\n",
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
