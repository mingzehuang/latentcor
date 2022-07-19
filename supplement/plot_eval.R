library(heatmaply)
library(ggplot2)
library(plotly)
library(graphics)
load("~/latentcor_git/latentcor/supplement/all_evaluation.rda")
grid_list = list(LatentR = seq(-0.9, 0.9, by = 0.1), TruncRate = seq(0.1, 0.9, by = 0.1))
eval_data = evaluation_BC$meanAE_diff
rownames(eval_data) = as.character(grid_list[[1]]); colnames(eval_data) = as.character(grid_list[[2]])
heatmaply(eval_data, dendrogram = "none", main = "Mean Absolute Approximation Error for Binary/Continuous", margins = c(80,80,80,80), ylab = names(grid_list[1]), xlab = names(grid_list[2]),
          grid_color = "white", grid_width = 0.00001, label_names = c(paste0(names(grid_list[2]), ":"), paste0(names(grid_list[1]), ":"), "Mean Absolute Approximation Error:"))

eval_data = evaluation_BB$meanAE_diff
rownames(eval_data) = as.character(grid_list[[1]]); colnames(eval_data) = as.character(grid_list[[2]])
heatmaply(eval_data, dendrogram = "none", main = "Mean Absolute Approximation Error for Binary/Binary", margins = c(80,80,80,80), ylab = names(grid_list[1]), xlab = names(grid_list[2]),
          grid_color = "white", grid_width = 0.00001, label_names = c(paste0(names(grid_list[2]), ":"), paste0(names(grid_list[1]), ":"), "Mean Absolute Approximation Error:"))

eval_data = evaluation_TC$meanAE_diff
rownames(eval_data) = as.character(grid_list[[1]]); colnames(eval_data) = as.character(grid_list[[2]])
heatmaply(eval_data, dendrogram = "none", main = "Mean Absolute Approximation Error for Truncated/Continuous", margins = c(80,80,80,80), ylab = names(grid_list[1]), xlab = names(grid_list[2]),
          grid_color = "white", grid_width = 0.00001, label_names = c(paste0(names(grid_list[2]), ":"), paste0(names(grid_list[1]), ":"), "Mean Absolute Approximation Error:"))

eval_data = evaluation_TB$meanAE_diff
rownames(eval_data) = as.character(grid_list[[1]]); colnames(eval_data) = as.character(grid_list[[2]])
heatmaply(eval_data, dendrogram = "none", main = "Mean Absolute Approximation Error for Truncated/Binary", margins = c(80,80,80,80), ylab = names(grid_list[1]), xlab = names(grid_list[2]),
          grid_color = "white", grid_width = 0.00001, label_names = c(paste0(names(grid_list[2]), ":"), paste0(names(grid_list[1]), ":"), "Mean Absolute Approximation Error:"))

eval_data = evaluation_TT$meanAE_diff
rownames(eval_data) = as.character(grid_list[[1]]); colnames(eval_data) = as.character(grid_list[[2]])
heatmaply(eval_data, dendrogram = "none", main = "Mean Absolute Approximation Error for Truncated/Truncated", margins = c(80,80,80,80), ylab = names(grid_list[1]), xlab = names(grid_list[2]),
          grid_color = "white", grid_width = 0.00001, label_names = c(paste0(names(grid_list[2]), ":"), paste0(names(grid_list[1]), ":"), "Mean Absolute Approximation Error:"))

eval_data = evaluation_NC$meanAE_diff
rownames(eval_data) = as.character(grid_list[[1]]); colnames(eval_data) = as.character(grid_list[[2]])
heatmaply(eval_data, dendrogram = "none", main = "Mean Absolute Approximation Error for Ternary/Continuous", margins = c(80,80,80,80), ylab = names(grid_list[1]), xlab = names(grid_list[2]),
          grid_color = "white", grid_width = 0.00001, label_names = c(paste0(names(grid_list[2]), ":"), paste0(names(grid_list[1]), ":"), "Mean Absolute Approximation Error:"))

eval_data = evaluation_NB$meanAE_diff
rownames(eval_data) = as.character(grid_list[[1]]); colnames(eval_data) = as.character(grid_list[[2]])
heatmaply(eval_data, dendrogram = "none", main = "Mean Absolute Approximation Error for Ternary/Binary", margins = c(80,80,80,80), ylab = names(grid_list[1]), xlab = names(grid_list[2]),
          grid_color = "white", grid_width = 0.00001, label_names = c(paste0(names(grid_list[2]), ":"), paste0(names(grid_list[1]), ":"), "Mean Absolute Approximation Error:"))

eval_data = evaluation_NT$meanAE_diff
rownames(eval_data) = as.character(grid_list[[1]]); colnames(eval_data) = as.character(grid_list[[2]])
heatmaply(eval_data, dendrogram = "none", main = "Mean Absolute Approximation Error for Ternary/Truncated", margins = c(80,80,80,80), ylab = names(grid_list[1]), xlab = names(grid_list[2]),
          grid_color = "white", grid_width = 0.00001, label_names = c(paste0(names(grid_list[2]), ":"), paste0(names(grid_list[1]), ":"), "Mean Absolute Approximation Error:"))

eval_data = evaluation_NN$meanAE_diff
rownames(eval_data) = as.character(grid_list[[1]]); colnames(eval_data) = as.character(grid_list[[2]])
heatmaply(eval_data, dendrogram = "none", main = "Mean Absolute Approximation Error for Ternary/Ternary", margins = c(80,80,80,80), ylab = names(grid_list[1]), xlab = names(grid_list[2]),
          grid_color = "white", grid_width = 0.00001, label_names = c(paste0(names(grid_list[2]), ":"), paste0(names(grid_list[1]), ":"), "Mean Absolute Approximation Error:"))

