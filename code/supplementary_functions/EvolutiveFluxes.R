
#' Variabilty of antibiotic therapy response is reflected by the bacterial metabolism
#' 
#' Varying the set of parameters in the dynamical model is reflected 
#' by the evolution of metabolic fluxes as represented by PCA over simulation time

library(ggplot2)
library(scatterplot3d)
library(plot3D)
library(tidyverse)
library(patchwork)

################################################################################
################## Evolution of metabolic fluxes (PCA) #########################
################################################################################

setwd("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Results")

conditions = c("NoDrug", "Therapy")

times = seq(0, 48, 1); vecolor = 1; vecolor_tsne = 1

listcol = list(c("yellow", "red", "black"), 
               c("cyan", "blue", "black"))

viridis = c("turbo", "mako")

for (con in conditions) {
  
  TracesPath = paste("./CDiffFBA", con, "_analysis", 
                     "/EpitCellDifficileHemeSink-analysis-1.trace", sep = "")
  
  trace = read.table(TracesPath, header = F)
  
  ReactionsNames = trace[1, ]
  ReactionsNames = gsub("\\(", replacement = "_", ReactionsNames)
  ReactionsNames = gsub("\\)", replacement = "", ReactionsNames)
  
  TracesPath = paste("./CDiffFBA", con, "_analysis", 
                     "/EpitCellDifficileHemeSink-analysis-1.trace", sep = "")
  
  trace = read.table(TracesPath, header = T)
  
  colnames(trace) = ReactionsNames
  
  subflux = trace[-1, c(1, 21:length(colnames(trace)))]
  
  subflux$Time = as.character(subflux$Time)
  rownames(subflux) <- subflux$Time
  subflux = subflux[, -1]
  
  # remove zero variance columns from the dataset
  subflux = subflux[, which(apply(subflux, 2, var) != 0)]
  
  ###############################################
  ######## Principal Component Analysis #########
  ###############################################
  
  princ <- prcomp(subflux, scale = TRUE)
  # plot(princ$x[, 1], princ$x[, 2])
  
  princ.var = princ$sdev^2
  princ.var.perc = round(princ.var/sum(princ.var)*100, 1)
  pc = c()
  
  for (i in 1:length(princ.var.perc)) {
    pc[i] = paste("PC", i, sep = ".")
  }
  
  princ.var.perc = data.frame(perc = princ.var.perc, pc)
  
  pca.bar = ggplot(princ.var.perc[1:5, ]) + geom_bar(aes(order(pc, perc), perc, 
                          fill = perc), stat = 'identity') + 
    viridis::scale_fill_viridis(direction = 1, option = viridis[vecolor]) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 0, hjust = 1)) + labs(
      title = paste("Components: ", con, sep = ""),
      y = "% Variance Explained", x = "Dimentions")
  
  pca.data = data.frame(Sample = rownames(princ$x), X = princ$x[,1], Y = princ$x[,2] )
  pca.data$Time = as.double(pca.data$Sample)
  
  x = matrix(NA, length(pca.data$X), 6)
  #x value; #y value; #z value
  x[, 1] = pca.data$X; x[, 2] = (pca.data$Time); x[, 3] = pca.data$Y
  #additional value
  x[, 4] = seq(0.5, 48, 0.5)
  #find out in which interval each additional value is
  intervals = seq(0.5, 48, 0.5)
  x[, 5] = findInterval(x[, 4], intervals)
  # produce gradient of colors
  gradient = colorRampPalette(colors = listcol[[vecolor]])
  colours = gradient(length(intervals))
  x[, 6] = colours[x[, 5]]
  
  # layout(matrix(1:2, ncol = 2), 
  #        width = c(3, 1), 
  #        height = c(1, 1))
  # scatterplot3d(x[,1], x[,3], x[,2], type = 'p', grid = F, box = F,
  #               cex.symbols = 1.25, color = x[, 6], pch = 16,
  #               xlab = "Time (hour)", angle = 40, 
  #               ylab = paste("PC1 - ", princ.var.perc$perc[1], "%", sep=""), 
  #               zlab = paste("PC2 - ", princ.var.perc$perc[2], "%", sep=""),
  #               main="Evolution of metabolic fluxes represented by PCA over simulation time")
  # FactoClass::addgrids3d(as.double(x[,1]), 
  #                        as.double(x[,3]), 
  #                        as.double(x[,2]), 
  #                        grid = c("xy", "xz", "yz"))
  # plot(x = rep(1, length(pca.data$X)), y = seq_along(x[,6]), 
  #      pch = 15, cex = 1.5, 
  #      col = gradient(length(x[,6])), 
  #      ann = F, axes = F, xlim = c(1, 2))
  # axis(side = 2, at = seq(1, nrow(x), length.out = 12), 
  #      labels = seq(1, 48, 4),
  #      line = 0.15)
  
  pdf(paste("3dplot", con, ".pdf", sep = ""), width = 18, height = 5, bg = "transparent")
  layout(matrix(1:2, ncol = 2), width = c(3, 1), height = c(0.5, 0.5))
  scatter3D(as.double(x[,1]), 
            as.double(x[,3]), 
            as.double(x[,2]), 
            colvar = 1:nrow(x), 
            col = x[, 6], colkey = F,
            xlab = "Time (hour)", angle = 90, 
            zlab = paste("PC1 - ", princ.var.perc$perc[1], "%", sep=""), 
            ylab = paste("PC2 - ", princ.var.perc$perc[2], "%", sep=""),
            main = "Evolution of metabolic fluxes represented by PCA over simulation time",
            theta = 30, phi = 10,
            pch = 19, # filled circles
            alpha = 0.75, # partially transparent
            bty = "b2")
  
  plot(x = rep(1, length(pca.data$X)), y = seq_along(x[,6]),
       pch = 15, cex = 1.5,
       col = gradient(length(x[,6])),
       ann = F, axes = F, xlim = c(1, 2)) + 
    axis(side = 2, at = seq(1, nrow(x), length.out = 12),
       labels = seq(1, 48, 4),
       line = 0.15)
  
  dev.off()
  
  pca.plt = ggplot(data = pca.data, 
                   aes(x = X, y = Y, color = Time)) + 
    xlab(paste("PC1 - ", princ.var.perc$perc[1], "%", sep="")) +
    ylab(paste("PC2 - ", princ.var.perc$perc[2], "%", sep="")) +
    geom_point() + theme_bw() +
    ggtitle("Evolution of metabolic fluxes \n over simulation time (PCA)") + 
    scale_colour_gradientn(
      limits  = range(pca.data$Time),
      colours = colours[c(1, seq_along(x[, 6]), length(x[, 6]))],
      values  = c(0, scales::rescale(as.double(x[, 4]), from = range(pca.data$Time)), 1))
  
  pdf(paste("2dplot", con, ".pdf", sep = ""),
      width = 18, height = 5, bg = "transparent")
  pca.bar + pca.plt
  ggsave(paste("2dplot", con, ".pdf", sep = ""),
         width = 8, height = 4, dpi = 500)
  
  vecolor = vecolor + 1
  
  ##############################################################
  ######## t-distributed stochastic neighbor embedding #########
  ##############################################################
  
  tsnePDX = tsne::tsne(subflux)
  tsnePDX = cbind(tsnePDX, as.double(pca.data$Sample))
  colnames(tsnePDX) = c("tsne1", "tsne2", "Time")
  tsnePDX = data.frame(tsnePDX)
  
  # pdf(paste("2dplot_tsne", con, ".pdf", sep = ""), 
  #     width = 18, height = 5, bg = "transparent")
  # 
  # ggplot(data = tsnePDX, aes(x = tsne1, y = tsne2, color = Time)) +
  #   geom_point() + theme_bw() + ggtitle("tsnePlot") + 
  #   scale_fill_continuous(limits = c(0, max(tsnePDX$Time)))
  # 
  # dev.off()
  
  tsne = matrix(NA, length(tsnePDX$tsne1), 6)
  tsne[, 1] = tsnePDX$tsne1
  tsne[, 2] = tsnePDX$Time
  tsne[, 3] = tsnePDX$tsne2
  tsne[, 4] = seq(0.5, 48, 0.5)
  
  tsne[, 5] = findInterval(tsne[, 4], seq(0.5, 48, 0.5))
  gradient = colorRampPalette(colors = listcol[[vecolor_tsne]])
  colours = gradient(length(intervals))
  tsne[, 6] = colours[tsne[, 5]]
  
  pdf(paste("3dplot_tsne", con, ".pdf", sep = ""), 
      width = 18, height = 5, bg = "transparent")
  
  layout(matrix(1:2, ncol = 2), width = c(3, 1), height = c(0.5, 0.5))
  scatter3D(as.double(tsne[,1]), 
            as.double(tsne[,3]), 
            as.double(tsne[,2]), 
            colvar = 1:nrow(x), 
            col = tsne[, 6], colkey = F,
            xlab = "Time (hour)", angle = 90, 
            zlab = "t-SNE_1", ylab = "t-SNE_2",
            main = "t-SNE visualization over simulation time",
            theta = 30, phi = 10,
            pch = 19, # filled circles
            alpha = 0.75, # partially transparent
            bty = "b2")
  
  plot(x = rep(1, length(tsnePDX$tsne1)), 
       y = seq_along(tsne[,6]), pch = 15, cex = 1.5,
       col = gradient(length(x[,6])), 
       ann = F, axes = F, xlim = c(1, 2)) + 
    axis(side = 2, at = seq(1, nrow(x), length.out = 12),
         labels = seq(1, 48, 4),
         line = 0.15)
  
  dev.off()
  
  vecolor_tsne = vecolor_tsne + 1
  
}
