library(heatmaply)
library(ggplot2)
library(plotly)
library(graphics)
library(gridExtra)
library(png)
library(grid)
load("~/latentcor/latentcor/supplement/all_evaluation.rda")
grid_list = list(LatentR = seq(-0.9, 0.9, by = 0.1), TruncRate = seq(0.1, 0.9, by = 0.1))
names(grid_list)=c("Latent Correlation", "Proportion of Zeros")
eval_data = evaluation_BC$meanAE_diff
rownames(eval_data) = as.character(grid_list[[1]]); colnames(eval_data) = as.character(grid_list[[2]])
heatmaply(eval_data, dendrogram = "none", main = "binary vs. continuous (hybrid)", margins = c(80,80,80,80), ylab = names(grid_list[1]), xlab = names(grid_list[2]),
          grid_color = "white", grid_width = 0.00001, label_names = c(paste0(names(grid_list[2]), ":"), paste0(names(grid_list[1]), ":"), "Mean Absolute Approximation Error:"))

eval_data = evaluation_BB$meanAE_diff
rownames(eval_data) = as.character(grid_list[[1]]); colnames(eval_data) = as.character(grid_list[[2]])
heatmaply(eval_data, dendrogram = "none", main = "binary vs. binary (hybrid)", margins = c(80,80,80,80), ylab = names(grid_list[1]), xlab = names(grid_list[2]),
          grid_color = "white", grid_width = 0.00001, label_names = c(paste0(names(grid_list[2]), ":"), paste0(names(grid_list[1]), ":"), "Mean Absolute Approximation Error:"))

eval_data = evaluation_TC$meanAE_diff
rownames(eval_data) = as.character(grid_list[[1]]); colnames(eval_data) = as.character(grid_list[[2]])
heatmaply(eval_data, dendrogram = "none", main = "truncated vs. continuous (hybrid)", margins = c(80,80,80,80), ylab = names(grid_list[1]), xlab = names(grid_list[2]),
          grid_color = "white", grid_width = 0.00001, label_names = c(paste0(names(grid_list[2]), ":"), paste0(names(grid_list[1]), ":"), "Mean Absolute Approximation Error:"))

eval_data = evaluation_TB$meanAE_diff
rownames(eval_data) = as.character(grid_list[[1]]); colnames(eval_data) = as.character(grid_list[[2]])
heatmaply(eval_data, dendrogram = "none", main = "truncated vs. binary (hybrid)", margins = c(80,80,80,80), ylab = names(grid_list[1]), xlab = names(grid_list[2]),
          grid_color = "white", grid_width = 0.00001, label_names = c(paste0(names(grid_list[2]), ":"), paste0(names(grid_list[1]), ":"), "Mean Absolute Approximation Error:"))

eval_data = evaluation_TT$meanAE_diff
rownames(eval_data) = as.character(grid_list[[1]]); colnames(eval_data) = as.character(grid_list[[2]])
heatmaply(eval_data, dendrogram = "none", main = "truncated vs. truncated (hybrid)", margins = c(80,80,80,80), ylab = names(grid_list[1]), xlab = names(grid_list[2]),
          grid_color = "white", grid_width = 0.00001, label_names = c(paste0(names(grid_list[2]), ":"), paste0(names(grid_list[1]), ":"), "Mean Absolute Approximation Error:"))

eval_data = evaluation_NC$meanAE_diff
rownames(eval_data) = as.character(grid_list[[1]]); colnames(eval_data) = as.character(grid_list[[2]])
heatmaply(eval_data, dendrogram = "none", main = "ternary vs. continuous (hybrid)", margins = c(80,80,80,80), ylab = names(grid_list[1]), xlab = names(grid_list[2]),
          grid_color = "white", grid_width = 0.00001, label_names = c(paste0(names(grid_list[2]), ":"), paste0(names(grid_list[1]), ":"), "Mean Absolute Approximation Error:"))

eval_data = evaluation_NB$meanAE_diff
rownames(eval_data) = as.character(grid_list[[1]]); colnames(eval_data) = as.character(grid_list[[2]])
heatmaply(eval_data, dendrogram = "none", main = "ternary vs. binary (hybrid)", margins = c(80,80,80,80), ylab = names(grid_list[1]), xlab = names(grid_list[2]),
          grid_color = "white", grid_width = 0.00001, label_names = c(paste0(names(grid_list[2]), ":"), paste0(names(grid_list[1]), ":"), "Mean Absolute Approximation Error:"))

eval_data = evaluation_NT$meanAE_diff
rownames(eval_data) = as.character(grid_list[[1]]); colnames(eval_data) = as.character(grid_list[[2]])
heatmaply(eval_data, dendrogram = "none", main = "truncated vs. ternary (hybrid)", margins = c(80,80,80,80), ylab = names(grid_list[1]), xlab = names(grid_list[2]),
          grid_color = "white", grid_width = 0.00001, label_names = c(paste0(names(grid_list[2]), ":"), paste0(names(grid_list[1]), ":"), "Mean Absolute Approximation Error:"))

eval_data = evaluation_NN$meanAE_diff
rownames(eval_data) = as.character(grid_list[[1]]); colnames(eval_data) = as.character(grid_list[[2]])
heatmaply(eval_data, dendrogram = "none", main = "ternary vs. ternary (hybrid)", margins = c(80,80,80,80), ylab = names(grid_list[1]), xlab = names(grid_list[2]),
          grid_color = "white", grid_width = 0.00001, label_names = c(paste0(names(grid_list[2]), ":"), paste0(names(grid_list[1]), ":"), "Mean Absolute Approximation Error:"))

plots <- lapply(ll <- c("~/latentcor/latentcor/supplement/BC_Hybrid.png", "~/latentcor/latentcor/supplement/BB_Hybrid.png", "~/latentcor/latentcor/supplement/TC_Hybrid.png"),function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})
ggsave("~/latentcor/latentcor/supplement/BC_BB_TC.pdf", marrangeGrob(grobs=plots, nrow=1, ncol=3, top = NULL), width = 15, height = 5)

plots <- lapply(ll <- c("~/latentcor/latentcor/supplement/TB_Hybrid.png", "~/latentcor/latentcor/supplement/TT_Hybrid.png", "~/latentcor/latentcor/supplement/TN_Hybrid.png"),function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})
ggsave("~/latentcor/latentcor/supplement/TB_TT_TN.pdf", marrangeGrob(grobs=plots, nrow=1, ncol=3, top = NULL), width = 15, height = 5)

plots <- lapply(ll <- c("~/latentcor/latentcor/supplement/NC_Hybrid.png", "~/latentcor/latentcor/supplement/NB_Hybrid.png", "~/latentcor/latentcor/supplement/NN_Hybrid.png"),function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})
ggsave("~/latentcor/latentcor/supplement/NC_NB_NN.pdf", marrangeGrob(grobs=plots, nrow=1, ncol=3, top = NULL), width = 15, height = 5)


eval_time = c(median(evaluation_BC$mediantime_1), median(evaluation_BC$mediantime_2),
                 median(evaluation_BB$mediantime_1), median(evaluation_BB$mediantime_2),
                 median(evaluation_TC$mediantime_1), median(evaluation_TC$mediantime_2),
                 median(evaluation_TB$mediantime_1), median(evaluation_TB$mediantime_2),
                 median(evaluation_TT$mediantime_1), median(evaluation_TT$mediantime_2),
                 median(evaluation_NT$mediantime_1), median(evaluation_NT$mediantime_2),
                 median(evaluation_NC$mediantime_1), median(evaluation_NC$mediantime_2),
                 median(evaluation_NB$mediantime_1), median(evaluation_NB$mediantime_2),
                 median(evaluation_NN$mediantime_1), median(evaluation_NN$mediantime_2)) / 10^6

comb = rep(c("BC", "BB", "TC", "TB", "TT", "TN", "NC", "NB", "NN"), each = 2)
method = rep(c("original", "approx"), 9)
time = data.frame(comb, time = eval_time, method)

timing_plot = ggplot(time, aes(x = comb, y = time, fill = method, group = method)) + geom_bar(stat="identity", position=position_dodge()) + ggtitle("Speed comparison") + xlab("Combination of data types") + ylab("Median computation time (in milliseconds)")

ggsave("~/latentcor/latentcor/supplement/computation_time.pdf", timing_plot, width = 10, height = 5)

plots <- lapply(ll <- c("~/latentcor/latentcor/supplement/speed1.png", "~/latentcor/latentcor/supplement/speed2.png"),function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})
ggsave("~/latentcor/latentcor/supplement/speed.pdf", marrangeGrob(grobs=plots, nrow=1, ncol=2, top = NULL), width = 15, height = 5)
