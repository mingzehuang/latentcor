# Script for example of mtcars
library(ggplot2)
mtcars$cyl = as.character(mtcars$cyl)
mtcars$vs = as.character(mtcars$vs)
mtcars$am = as.character(mtcars$am)
hist_mpg = ggplot(mtcars, aes(x=mpg)) + geom_histogram(aes(y=..density..), binwidth = 2, color="black", fill="white")+geom_density(alpha=.2, fill="#FF6666")
hist_cyl = ggplot(mtcars, aes(cyl)) + geom_bar(color="black", fill="white")
hist_disp = ggplot(mtcars, aes(x=disp)) + geom_histogram(aes(y=..density..), binwidth = 40, color="black", fill="white")+geom_density(alpha=.2, fill="#FF6666")
hist_hp = ggplot(mtcars, aes(x=hp)) + geom_histogram(aes(y=..density..), binwidth = 20, color="black", fill="white")+geom_density(alpha=.2, fill="#FF6666")
hist_drat = ggplot(mtcars, aes(x=drat)) + geom_histogram(aes(y=..density..), binwidth = .2, color="black", fill="white")+geom_density(alpha=.2, fill="#FF6666")
hist_wt = ggplot(mtcars, aes(x=wt)) + geom_histogram(aes(y=..density..), binwidth = .5, color="black", fill="white")+geom_density(alpha=.2, fill="#FF6666")
hist_qsec = ggplot(mtcars, aes(x=qsec)) + geom_histogram(aes(y=..density..), binwidth = 1, color="black", fill="white")+geom_density(alpha=.2, fill="#FF6666")
hist_vs = ggplot(mtcars, aes(vs)) + geom_bar(color="black", fill="white")
hist_am = ggplot(mtcars, aes(am)) + geom_bar(color="black", fill="white")
hist_gear = ggplot(mtcars, aes(gear)) + geom_bar(color="black", fill="white")
hist_carb = ggplot(mtcars, aes(x=carb)) + geom_histogram(aes(y=..density..), binwidth = 1, color="black", fill="white")+geom_density(alpha=.2, fill="#FF6666")
plot_all_1 = gridExtra::marrangeGrob(grobs = list(hist_mpg, hist_cyl, hist_disp, hist_hp, hist_drat, hist_wt), widths = rep(1, 2), top = NULL, layout_matrix = matrix(c(1:6), 3, 2, byrow = T))
ggsave("hist_mtcars_1.pdf", plot_all_1, width = 10, height = 15)
plot_all_2 = gridExtra::marrangeGrob(grobs = list(hist_qsec, hist_vs, hist_am, hist_gear, hist_carb), widths = rep(1, 2), top = NULL, layout_matrix = matrix(c(1:5, NA), 3, 2, byrow = T))
ggsave("hist_mtcars_2.pdf", plot_all_2, width = 10, height = 15)

plot_12 = gridExtra::marrangeGrob(grobs = list(hist_mpg, hist_cyl), widths = rep(1, 2), top = NULL, layout_matrix = matrix(c(1:2), 1, 2, byrow = T))
ggsave('hist_mtcars_12.pdf', plot_12, width = 10, height = 5)

rm(list=ls())
library(latentcor)
library(heatmaply)
library(ggplot2)
library(plotly)
library(graphics)
library(gridExtra)
library(png)
library(grid)
mtcars$cyl[mtcars$cyl == 4] = 0; mtcars$cyl[mtcars$cyl == 6] = 1; mtcars$cyl[mtcars$cyl == 8] = 2
mtcars$gear[mtcars$gear == 3] = 0; mtcars$gear[mtcars$gear == 4] = 1; mtcars$gear[mtcars$gear == 5] = 2
mtcars_pearson = cor(mtcars)
heatmaply(mtcars_pearson, dendrogram = "none", main = "Pearson correlation", margins = c(80,80,80,80),
                                grid_color = "white", grid_width = 0.00001)

mtcars_latentcor = latentcor(X = mtcars, types = c("con", "ter", "con", "con", "con", "con", "con", "bin", "bin", "ter", "con"), method = "original")$R
heatmaply(mtcars_latentcor, dendrogram = "none", main = "Estimated latent correlation", margins = c(80,80,80,80),
                            grid_color = "white", grid_width = 0.00001)

mtcars_diff = mtcars_latentcor - mtcars_pearson
heatmaply(mtcars_diff, dendrogram = "none", main = "Difference", col = cool_warm(50), margins = c(80,80,80,80),
                             grid_color = "white", grid_width = 0.00001)


plots <- lapply(ll <- c("~/latentcor/latentcor/supplement/mtcars_pearson.png", "~/latentcor/latentcor/supplement/mtcars_latentcor.png", "~/latentcor/latentcor/supplement/mtcars_diff.png"),function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})
ggsave("~/latentcor/latentcor/supplement/mtcars_all_heatmap.pdf", marrangeGrob(grobs=plots, nrow=1, ncol=3, top = NULL), width = 15, height = 5)
