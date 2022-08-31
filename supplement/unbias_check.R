library(latentcor)
library(ggplot2)
library(gridExtra)

rm(list=ls())
PlotPair = function(datapair, namepair = c("X", "Y"), title = "Plot X vs Y") {
  df = data.frame(datapair)
  colnames(df) = namepair
  print(ggplot(df, aes(x = datapair[ , 1], y = datapair[ , 2]))
        + geom_point(color = "blue") + geom_abline(intercept = 0, slope = 1, color = "red")
        +ggtitle(title) + xlab(namepair[1]) + ylab(namepair[2]) + xlim(-1, 1) + ylim(-1, 1))
}

cp1 = "expo"; cp2 = "cube"
for (tp1 in c("con", "bin", "ter", "tru")) {
  for (tp2 in c("con", "bin", "ter", "tru")) {
    rhorep = seq(-1, 1, by = .01); Rrep = matrix(NA, length(rhorep), 4)
    for (rep in 1:length(rhorep)) {
      rho = rhorep[rep]
      X = gen_data(n = 1000, rhos = rho, types = c(tp1, tp2), copulas = c(cp1, cp2))$X
      R_org = latentcor(X, types = c(tp1, tp2), method = "original")$R
      R_approx = latentcor(X, types = c(tp1, tp2), method = "approx")$R
      R_conapprox = latentcor(X, types = c(ifelse(tp1 == "ter", "con", tp1), ifelse(tp2 == "ter", "con", tp2)), method = "original")$R
      Rrep[rep, 1] = R_org[2, 1]; Rrep[rep, 2] = R_approx[2, 1]; Rrep[rep, 3] = cor(X)[2, 1]; Rrep[rep, 4] = R_conapprox[2, 1]
    }
    assign(paste("R", cp1, cp2, tp1, tp2, "org", sep = "_"), Rrep[ , 1])
    assign(paste("R", cp1, cp2, tp1, tp2, "approx", sep = "_"), Rrep[ , 2])
    assign(paste("R", cp1, cp2, tp1, tp2, "pearson", sep = "_"), Rrep[ , 3])
    assign(paste("R", cp1, cp2, tp1, tp2, "conapprox", sep = "_"), Rrep[ , 4])

    if (tp1 == "con") {type1 = "continuous"}
    else if (tp1 == "bin") {type1 = "binary"}
    else if (tp1 == "tru") {type1 = "truncated"}
    else if (tp1 == "ter") {type1 = "ternary"}
    if (tp2 == "con") {type2 = "continuous"}
    else if (tp2 == "bin") {type2 = "binary"}
    else if (tp2 == "tru") {type2 = "truncated"}
    else if (tp2 == "ter") {type2 = "ternary"}

    plot_org = PlotPair(datapair = cbind(rhorep, c(get(paste("R", cp1, cp2, tp1, tp2, "org", sep = "_")))),
             namepair = c("True latent correlation", "Estimated latent correlation (original)"),
             title = paste(type1, " vs. ", type2, " (", cp1, " vs. ", cp2, ")", sep = ""))

    plot_approx = PlotPair(datapair = cbind(rhorep, c(get(paste("R", cp1, cp2, tp1, tp2, "approx", sep = "_")))),
             namepair = c("True latent correlation", "Estimated latent correlation (approx)"),
             title = paste(type1, " vs. ", type2, " (", cp1, " vs. ", cp2, ")", sep = ""))

    plot_pearson = PlotPair(datapair = cbind(rhorep, c(get(paste("R", cp1, cp2, tp1, tp2, "pearson", sep = "_")))),
             namepair = c("True latent correlation", "Pearson correlation"),
             title = paste(type1, " vs. ", type2, " (", cp1, " vs. ", cp2, ")", sep = ""))

    plot_all = gridExtra::marrangeGrob(grobs = list(plot_pearson, plot_org, plot_approx), widths = c(10, 10, 10), top = NULL, layout_matrix = matrix(c(1, 2, 3), 1, 3))
    ggsave(paste(tp1, "vs.", tp2, cp1, cp2, ".pdf", sep = " "), plot_all, width = 15, height = 5)
    plot_pearson_org = gridExtra::marrangeGrob(grobs = list(plot_pearson, plot_org), widths = c(10, 10), top = NULL, layout_matrix = matrix(c(1, 2), 1, 2))
    ggsave(paste(tp1, "vs.", tp2, cp1, cp2, "pearson_org.pdf", sep = " "), plot_pearson_org, width = 10, height = 5)
    # pdf(file = paste(tp1, "vs.", tp2, cp1, cp2, "(conapprox).pdf", sep = " "), width = 5, height = 5)
    # PlotPair(datapair = cbind(rhorep, c(get(paste("R", cp1, cp2, tp1, tp2, "conapprox", sep = "_")))),
    #          namepair = c("True latent correlation", "Estimated latent correlation (approximated by continuous)"),
    #          title = paste(type1, " vs. ", type2, " (", cp1, " vs. ", cp2, ")", sep = ""))
    # dev.off()
  }
}


cp1 = "no"; cp2 = "no"
tp1 = "con"; tp2 = "con"; rhorep = seq(-1, 1, by = .01);  Rrep = matrix(NA, length(rhorep), 2)
for (rep in 1:length(rhorep)) {
  rho = rhorep[rep]
  X = gen_data(n = 1000, rhos = rho, types = c(tp1, tp2), copulas = c(cp1, cp2))$X
  rhorep[rep] = rho; Rrep[rep, 1] = cor(X)[2, 1]; Rrep[rep, 2] = pcaPP::cor.fk(X)[2, 1]
  }
plot_pearson = PlotPair(datapair = cbind(rhorep, Rrep[ , 1]),
                            namepair = c("True latent correlation", "Pearson correlation"),
                            title = "continuous vs. continuous")
plot_kendall = PlotPair(datapair = cbind(rhorep, Rrep[ , 2]),
                                namepair = c("True latent correlation", paste0("Kendall's ", expression(tau))),
                                title = "continuous vs. continuous")
cp1 = "expo"; cp2 = "cube"
tp1 = "con"; tp2 = "con"; rhorep = seq(-1, 1, by = .01);  Rrep = matrix(NA, length(rhorep), 2)
for (rep in 1:length(rhorep)) {
  rho = rhorep[rep]
  X = gen_data(n = 1000, rhos = rho, types = c(tp1, tp2), copulas = c(cp1, cp2))$X
  rhorep[rep] = rho; Rrep[rep, 1] = cor(X)[2, 1]; Rrep[rep, 2] = pcaPP::cor.fk(X)[2, 1]
}
plot_pearson_exp_cub = PlotPair(datapair = cbind(rhorep, Rrep[ , 1]),
                                namepair = c("True latent correlation", "Pearson correlation"),
                                title = "continuous vs. continuous (expo vs. cube)")

plot_kendall_exp_cub = PlotPair(datapair = cbind(rhorep, Rrep[ , 2]),
                                        namepair = c("True latent correlation", paste0("Kendall's ", expression(tau))),
                                        title = "continuous vs. continuous (expo vs. cube)")
plot_pearson = gridExtra::marrangeGrob(grobs = list(plot_pearson, plot_pearson_exp_cub), widths = c(10, 10), top = NULL, layout_matrix = matrix(c(1, 2), 1, 2))
    ggsave("con_con_pearson.pdf", plot_pearson, width = 10, height = 5)
plot_Kendall = gridExtra::marrangeGrob(grobs = list(plot_kendall, plot_kendall_exp_cub), widths = c(10, 10), top = NULL, layout_matrix = matrix(c(1, 2), 1, 2))
    ggsave("con_con_kendall.pdf", plot_Kendall, width = 10, height = 5)


#
# cp1 = "cube"; cp2 = "cube"
# for (tp1 in c("con")) {
#   for (tp2 in c("con", "bin", "tru", "ter", "qua", "qui", "sen", "sep", "oct", "nov", "den", "dtr")) {
#     rhorep = rep(NA, 100); Rrep = matrix(NA, 100, 3)
#     for (rep in 1:100) {
#       rho = runif(1, -1, 1)
#       X = GenData(n = 1000, rhos = rho, types = c(tp1, tp2), copulas = c(cp1, cp2))$X
#       R_org = estR(X, types = c(tp1, tp2), method = "original")$R
#       R_conapprox = estR(X, types = c("con", "con"), method = "original")$R
#       rhorep[rep] = rho; Rrep[rep, 1] = R_org[2, 1]; Rrep[rep, 2] = R_conapprox[2, 1]; Rrep[rep, 3] =  cor(X)[2, 1]
#     }
#     assign(paste("R", cp1, cp2, tp1, tp2, "org", sep = "_"), Rrep[ , 1])
#     assign(paste("R", cp1, cp2, tp1, tp2, "conapprox", sep = "_"), Rrep[ , 2])
#     assign(paste("R", cp1, cp2, tp1, tp2, "pearson", sep = "_"), Rrep[ , 3])
#     PlotPair(datapair = cbind(rhorep, c(get(paste("R", cp1, cp2, tp1, tp2, "org", sep = "_")))),
#              namepair = c("rho", paste("R", cp1, cp2, tp1, tp2, "org", sep = "_")),
#              title = paste(tp1, "vs.", tp2, "(org)", sep = " "))
#     PlotPair(datapair = cbind(rhorep, c(get(paste("R", cp1, cp2, tp1, tp2, "conapprox", sep = "_")))),
#              namepair = c("rho", paste("R", cp1, cp2, tp1, tp2, "conapprox", sep = "_")),
#              title = paste(tp1, "vs.", tp2, "(conapprox)", sep = " "))
#     PlotPair(datapair = cbind(rhorep, c(get(paste("R", cp1, cp2, tp1, tp2, "pearson", sep = "_")))),
#              namepair = c("rho", paste("R", cp1, cp2, tp1, tp2, "pearson", sep = "_")),
#              title = paste(tp1, "vs.", tp2, "(pearson)", sep = " "))
#   }
# }

pdf(file = "", width = , height = )
dev.off()

## Speed comparison
library(microbenchmark)
library(latentcor)
# p = 20
types_20 = rep(c("con", "bin", "ter", "tru"), 5)
X_20 = gen_data(types = types_20)$X
# p = 40
types_40 = rep(c("con", "bin", "ter", "tru"), 10)
X_40 = gen_data(types = types_40)$X
# p = 100
types_100 = rep(c("con", "bin", "ter", "tru"), 25)
X_100 = gen_data(types = types_100)$X
# p = 200
types_200 = rep(c("con", "bin", "ter", "tru"), 50)
X_200 = gen_data(types = types_200)$X
# p = 400
types_400 = rep(c("con", "bin", "ter", "tru"), 100)
X_400 = gen_data(types = types_400)$X
library(parallel)
library(doFuture)
plan(multisession, workers = availableCores())
timing = microbenchmark(latentcor(X = X_20, types = types_20, method = "original"), latentcor(X = X_20, types = types_20),
                        latentcor(X = X_40, types = types_40, method = "original"), latentcor(X = X_40, types = types_40),
                        latentcor(X = X_100, types = types_100, method = "original"), latentcor(X = X_100, types = types_100),
                        latentcor(X = X_200, types = types_200, method = "original"), latentcor(X = X_200, types = types_200),
                        latentcor(X = X_400, types = types_400, method = "original"), latentcor(X = X_400, types = types_400), times = 5L, unit = "s")
dim = rep(c(log10(20), log10(40), log10(100), log10(200), log10(400)), each = 2)
method = rep(c("original", "approx"), 5)
time = data.frame(dim, time = log10(c(1.5571201, 0.0063703, 6.1404127, 0.0184289, 42.2046845, 0.1080221, 168.5614688, 0.9871049, 707.5246815, 8.8118233)), method)

library(ggplot2)

timing_plot = print(ggplot(time, aes(x = dim, y = time, color = method, group = method)) + geom_line() + geom_point()
      + ggtitle("Speed comparison") + xlab("log10 of dimension") + ylab("log10 of computation time (in seconds)")
+ scale_x_continuous(breaks=c(log10(20), log10(40), log10(100), log10(200), log10(400)), labels=c("log10(20)", "log10(40)", "log10(100)", "log10(200)", "log10(400)")))
save(time, file = "timing.rda")
pdf(file = "./supplement/timing_plot.pdf", width = 21, height = 15)
timing_plot
dev.off()
