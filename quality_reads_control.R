
#from plotQualityProfile function

data_path <- '/Users/cigom/metagenomics/COI/run15/'
fl <- list.files(data_path, pattern = 'fastq.gz', full.names = TRUE)
f <- fl

n = 5e+05
FIRST <- TRUE
statdf <- data.frame(Cycle = integer(0), Mean = numeric(0), 
                     Q25 = numeric(0), Q50 = numeric(0), Q75 = numeric(0), 
                     Cum = numeric(0), file = character(0))
anndf <- data.frame(minScore = numeric(0), label = character(0), 
                    rclabel = character(0), rc = numeric(0), file = character(0))


for (f in fl[!is.na(fl)]) {
  library(ShortRead)
  srqa <- qa(f, n = n)
  df <- srqa[["perCycle"]]$quality
  rc <- sum(srqa[["readCounts"]]$read)
  if (rc >= n) {
    rclabel <- paste("Reads >= ", n)
  } else {
    rclabel <- paste("Reads: ", rc)
  }
  means <- rowsum(df$Score * df$Count, df$Cycle)/rowsum(df$Count, 
                                                        df$Cycle)
  get_quant <- function(xx, yy, q) {
    xx[which(cumsum(yy)/sum(yy) >= q)][[1]]
  }
  q25s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                   foo$Count, 0.25), simplify = TRUE)
  q50s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                   foo$Count, 0.5), simplify = TRUE)
  q75s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                   foo$Count, 0.75), simplify = TRUE)
  cums <- by(df, df$Cycle, function(foo) sum(foo$Count), 
             simplify = TRUE)
  if (!all(sapply(list(names(q25s), names(q50s), names(q75s), 
                       names(cums)), identical, rownames(means)))) {
    stop("Calculated quantiles/means weren't compatible.")
  }
  if (FIRST) {
    plotdf <- cbind(df, file = basename(f))
    FIRST <- FALSE
  } else {
    plotdf <- rbind(plotdf, cbind(df, file = basename(f)))
  }
  statdf <- rbind(statdf, data.frame(Cycle = as.integer(rownames(means)), 
                                     Mean = means, Q25 = as.vector(q25s), Q50 = as.vector(q50s), 
                                     Q75 = as.vector(q75s), Cum = 10 * as.vector(cums)/rc, 
                                     file = basename(f)))
  anndf <- rbind(anndf, data.frame(minScore = min(df$Score), 
                                   label = basename(f), rclabel = rclabel, rc = rc, 
                                   file = basename(f)))
}

anndf$minScore <- min(anndf$minScore)
plotdf.summary <- aggregate(Count ~ Cycle + Score, plotdf, sum)
plotdf.summary$label <- paste(nrow(anndf), "files (aggregated)")

means <- rowsum(plotdf.summary$Score * plotdf.summary$Count, 
                plotdf.summary$Cycle)/rowsum(plotdf.summary$Count, 
                                             plotdf.summary$Cycle)

q25s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, 
                                                                         foo$Count, 0.25), simplify = TRUE)
q50s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, 
                                                                         foo$Count, 0.5), simplify = TRUE)
q75s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, 
                                                                         foo$Count, 0.75), simplify = TRUE)
cums <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) sum(foo$Count), 
           simplify = TRUE)

statdf.summary <- data.frame(Cycle = as.integer(rownames(means)), 
                             Mean = means, Q25 = as.vector(q25s), 
                             Q50 = as.vector(q50s), 
                             Q75 = as.vector(q75s), 
                             Cum = 10 * as.vector(cums)/rc)

# 2c7fb8 # blue
# de2d26 # red
# FC8D62 # orange 
# 66C2A5 # green
p <- ggplot(data = plotdf.summary, aes(x = Cycle, y = Score)) + 
  geom_tile(aes(fill = Count)) + scale_fill_gradient(low = "#F5F5F5", 
                                                     high = "black") +
  geom_line(data = statdf.summary, aes(y = Mean), color = "#2c7fb8") +  
  geom_line(data = statdf.summary, aes(y = Q25), color = "#de2d26", size = 0.7, linetype = "dashed") + 
  geom_line(data = statdf.summary, aes(y = Q50), color = "#de2d26", size = 0.7) + 
  geom_line(data = statdf.summary, aes(y = Q75), color = "#de2d26", size = 0.7, linetype = "dashed") + 
  ylab("Quality Score") + xlab("Cycle") + 
  annotate("text", x = 0, y = 0, label = sprintf("Total reads: %d", sum(anndf$rc)), color = "red", hjust = 0) +   
  theme_classic(base_size = 16) + 
  theme(panel.grid = element_blank()) + guides(fill = FALSE) + 
  facet_wrap(~label) + ylim(c(0, NA))

p + geom_vline(xintercept=220, linetype="dashed", 
           color = "#636363", size=0.4)


if (length(unique(statdf$Cum)) > 1) {
  p <- p + geom_line(data = statdf, aes(y = Cum), color = "red", 
                     size = 0.25, linetype = "solid") + scale_y_continuous(sec.axis = sec_axis(~. * 
                                                                                                 10, breaks = c(0, 100), labels = c("0%", "100%"))) + 
    theme(axis.text.y.right = element_text(color = "red"), 
          axis.title.y.right = element_text(color = "red"))
}
# plot ----
p
# or
srqa <- qa(f, n = n)
dim(df <- srqa[["perCycle"]]$quality)
rc <- sum(srqa[["readCounts"]]$read)
dim(fqs <- srqa[['frequentSequences']])
fqs$size <- nchar(as.character(fqs$sequence))

# and
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("Rqc")

library(Rqc)

# https://bioconductor.org/packages/release/bioc/vignettes/Rqc/inst/doc/Rqc.html

files <- sort(list.files(dir, "fastq.gz", full.names=TRUE))
qa <- rqcQA(files, pair = c(1,1), workers = 1, group = rep('Mock',length(files)))

rqcReadQualityBoxPlot(qa)
rqcReadQualityPlot(qa)
plot <- rqcCycleAverageQualityPlot(qa); plot + theme_bw() # good
rqcReadFrequencyPlot(qa) 
rqcFileHeatmap(qa[[2]]) # good
rqcReadWidthPlot(qa)
rqcCycleGCPlot(qa)
rqcCycleQualityPlot(qa) #good 

#rqcCycleAverageQualityPcaPlot(qa)
rqcCycleQualityBoxPlot(qa) #good
rqcCycleBaseCallsPlot(qa)


# Qual_vs_MaxEE_plot
data_path <- '/Users/cigom/metagenomics/Franccesco/callahan_multirun/fastqc_clipped/'
#data_path <- '/Users/cigom/metagenomics/COI/run012/quality/no_group/'
fileqc <- list.files(data_path, pattern = '*zip', full.names = TRUE)

# fastq.r1 <- fileqc[1]
# fastq.r2 <- fileqc[1]
#source('https://raw.githubusercontent.com/RJEGR/infovis/master/qualMaxEEplot.R')

#qualMaxEEplot(fastq.r1 = fastq.r1, fastq.r2 = fastq.r2)

run <- sapply(strsplit(basename(fileqc), "[-]"), `[`, 1)

fastqr <- vector("list", length(unique(run)))
names(fastqr) <- unique(run)


set.seed(202007)
for(runs in unique(run)) {
  which_samples <- which(run %in% runs)
  fastqr[[runs]] <- fileqc[sample(which_samples, 10)]
}


#fastqr <- fileqc[sample(which_samples,1)]
tag <- function(filename) {sapply(strsplit(basename(filename), "[_]"), `[`, 1)}
basenameF <- unlist(lapply(fastqr, tag))
files <- fileqc[which(unlist(lapply(fileqc, tag)) %in% basenameF)]

# check datavis ----
raw_path <- '/Users/cigom/metagenomics/Franccesco/callahan_multirun/raw/'
fq_raw <- list.files(raw_path, pattern = '*gz', full.names = TRUE)
fq_raw <- fq_raw[tag(fq_raw) %in%  tag(files)]
group <- sapply(strsplit(basename(fq_raw), "[-]"), `[`, 1)
pn <- length(fq_raw) /2
pair <- rep(1:pn, each = 2)
qa <- rqcQA(fq_raw, pair = pair, workers = 1, group = group)

ggsave(rqcReadQualityBoxPlot(qa), path = raw_path, 
       filename = 'ReadQualityBoxPlo.png', width = 11)

ReadQualityPlot <- rqcReadQualityCalc(qa) %>%
  mutate(group = sapply(strsplit(as.character(filename), "[-]"), `[`, 1)) %>%
  ggplot(aes_string(y = "quantile", x = "value", 
                    colour = "group")) + 
  geom_point(alpha = 3/5, shape = 4) + 
  geom_smooth(aes(group = group, fill = group), method = "auto", se = T) +
  guides(fill = 'none') +
  labs(y = "% of Reads Exceeding Quality (q)", 
                     x = "Quality (q)", colour = "Run",
       subtitle = 'per-sample Quality fitted to y ~ x  polynomial model')

ggsave(ReadQualityPlot, path = raw_path, 
                filename = 'ReadQualityPlot.png', width = 11)

len <- max(as.integer(rqcCycleAverageQualityCalc(qa)$cycle))

CycleAverageQuality <- rqcCycleAverageQualityCalc(qa) %>%
  ggplot(aes_string(x = "cycle", y = "quality",colour = "group")) + 
  geom_point(alpha = 1/5) + 
  # geom_line(aes_string(group = "filename")) + 
  stat_smooth(aes(group = group, fill = group), se = T) +
  labs(x = "Cycle", y = "Average Quality", colour = "Run",
       subtitle = 'per-sample Quality fitted to y ~ x  polynomial model') + 
  scale_y_continuous(limits = c(34,40)) +
  guides(fill = 'none') +
  scale_x_discrete(breaks = seq(from = 1, to = len, by = len%/%20))

ggsave(CycleAverageQuality, path = raw_path, 
       filename = 'CycleAverageQuality.png', width = 11)


df <- rqcCycleQualityCalc(qa)
len <- max(as.integer(df$cycle))

CycleQualityPlot <- df %>%
  mutate(group = sapply(strsplit(as.character(filename), "[-]"), `[`, 1)) %>%
  ggplot(aes_string(x = "cycle", y = "percentiles", fill = "value")) + 
  geom_tile() + viridis::scale_fill_viridis() + 
    facet_grid(group ~.) + 
  labs(x = "Cycle", y = "%", fill = "Quality") + 
    scale_x_discrete(breaks = seq(from = 1, to = len, by = len%/%10)) +
  guides(fill = guide_colorbar(barheight = unit(4, "in"),
                               ticks.colour = "black", 
                               frame.colour = "black",
                               label.theme = element_text(size = 12)))
ggsave(CycleQualityPlot, path = raw_path, 
       filename = 'CycleQualityPlot.png', width = 11)


# continue ----
# fastqc --nogroup yourFastqFile_R1.fastq yourFastqFile_R2.fastq
fastqStats <- function(fileqc) {
  fastq <- data.frame(qc_read(fileqc, modules = "Per base sequence quality", verbose = TRUE)$per_base_sequence_quality)
  fastq.tmp <- rbind(data.frame(R=fastq$Base, 
                                Q=fastq$Mean, S=c("Mean"), 
                                E=10^(-fastq$Mean/10), 
                                A=Reduce('+', 10^(-fastq$Mean/10), accumulate = TRUE)),
                     data.frame(R=fastq$Base, 
                                Q=fastq$Median, 
                                S=c("Median"), 
                                E=10^(-fastq$Median/10), 
                                A=Reduce('+', 10^(-fastq$Median/10), accumulate = TRUE)),
                     data.frame(R=fastq$Base, 
                                Q=fastq$Lower.Quartile, 
                                S=c("Lower.Quartile"), 
                                E=10^(-fastq$Lower.Quartile/10), 
                                A=Reduce('+', 10^(-fastq$Lower.Quartile/10), accumulate = TRUE)),
                     data.frame(R=fastq$Base, 
                                Q=fastq$Upper.Quartile, 
                                S=c("Upper.Quartile"), 
                                E=10^(-fastq$Upper.Quartile/10), 
                                A=Reduce('+', 10^(-fastq$Upper.Quartile/10), accumulate = TRUE)),
                     data.frame(R=fastq$Base, 
                                Q=fastq$X10th.Percentile, 
                                S=c("X10th.Percentile"), 
                                E=10^(-fastq$X10th.Percentile/10), 
                                A=Reduce('+', 10^(-fastq$X10th.Percentile/10), accumulate = TRUE)),
                     data.frame(R=fastq$Base, 
                                Q=fastq$X90th.Percentile, 
                                S=c("X90th.Percentile"), 
                                E=10^(-fastq$X90th.Percentile/10), 
                                A=Reduce('+', 10^(-fastq$X90th.Percentile/10), accumulate = TRUE)))
  
  fastq.tmp <- data.frame(fastq.tmp, filename = basename(fileqc))
  return(fastq.tmp)
}
fastTemp <- function(fileqc) {
  fastq.tmp <- data.frame(qc_read(fileqc, modules = "Per base sequence quality", verbose = TRUE)$per_base_sequence_quality)
  
  fastq.tmp <- data.frame(fastq.tmp, filename = basename(fileqc))
  return(fastq.tmp)
}
  
fqdata <- do.call(rbind, lapply(files, fastqStats))
fqTmp <- do.call(rbind, lapply(files, fastTemp))


fqdata %>% 
  mutate(group = sapply(strsplit(as.character(filename), "[-]"), `[`, 1)) %>%
  mutate(filename = sapply(strsplit(as.character(filename), "[_]"), `[`, 1)) %>%
  ggplot(aes(color = S, group = filename)) + 
  geom_point(aes(x = R, y = Q), size=1, alpha = 3/5) + 
  labs(x="Reads position", y="Reads Quality", color = 'Stats') +
  # scale_color_brewer(palette = 'Dark2') +
  theme_bw() + facet_grid(group ~ .)

fqdata %>% 
  as_tibble() %>%
  filter(S %in% 'Mean' & !is.na(S)) %>%
  mutate(group = sapply(strsplit(as.character(filename), "[-]"), `[`, 1)) %>%
  mutate(filename = sapply(strsplit(as.character(filename), "[_]"), `[`, 1)) %>%
  group_by(group) %>%
  summarise(mean = mean(log10(A)), 
            Min = min(log10(A)), 
            Max = max(log10(A)), 
            sd = sd(log10(A)), 
            var = var(log10(A)))

# Error
fqdata %>% 
  mutate(group = sapply(strsplit(as.character(filename), "[-]"), `[`, 1)) %>%
  mutate(filename = sapply(strsplit(as.character(filename), "[_]"), `[`, 1)) %>%
  filter(S %in% 'Mean' & !is.na(S)) %>%
  ggplot() +
  geom_point(aes(x = R, y = log10(A), color = group), alpha = 1/5) +
  stat_smooth(aes(group = group, x = R, y = log10(A), color = group), se = F) +
  #coord_cartesian(ylim = c(-3, -1.5)) +
  labs(x="Reads position", y="EE = sum(10^(-Q/10)) log10",
       subtitle = 'per-sample Quality fitted to y ~ x  polynomial model')
#  or
stats <- c('Mean', 'Lower.Quartile', 'Upper.Quartile')
errplot <- fqdata %>% 
  mutate(group = sapply(strsplit(as.character(filename), "[-]"), `[`, 1)) %>%
  mutate(filename = sapply(strsplit(as.character(filename), "[_]"), `[`, 1)) %>%
  filter(S %in% stats & !is.na(S)) %>%
  ggplot(aes(x=R, y = Q)) +
  stat_summary_2d(aes(z = log10(A))) +
  geom_line(aes(y = Q, color = S), size = 0.4, linetype = "dashed") + #
  scale_color_manual(values = c("#e31a1c", "#1f78b4", "#ff7f00")) +
  viridis::scale_fill_viridis() +
  theme_classic() +
  theme(legend.position = "right") +
  facet_grid(group ~ .) +
  labs(x="Reads position", y="Reads Quality", color = 'Stats', fill = 'Expected\nError') +
  guides(fill = guide_colorbar(barheight = unit(4, "in"),
                               ticks.colour = "black", 
                               frame.colour = "black",
                               label.theme = element_text(size = 12)))

ggsave(errplot, path = raw_path, 
       filename = 'qualMaxEEPlot.png', width = 11, height = 7)

# grid search method for optimization in r
# https://machinelearningmastery.com/tuning-machine-learning-models-using-the-caret-r-package/
# issue in 
# https://github.com/benjjneb/dada2/issues/236
# http://topepo.github.io/caret/index.html

library(caret)
set.seed(7)
data(iris)

featurePlot(x = iris[, 1:4], 
            y = iris$Species, 
            plot = "pairs",
            ## Add a key at the top
            auto.key = list(columns = 3))

# scatterplot matrix with ellipses
featurePlot(x = iris[, 1:4], 
            y = iris$Species, 
            plot = "ellipse",
            ## Add a key at the top
            auto.key = list(columns = 3))

# Overlayed Density Plots
featurePlot(x = iris[, 1:4], 
            y = iris$Species,
            plot = "density", 
            ## Pass in options to xyplot() to 
            ## make it prettier
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")), 
            adjust = 1.5, 
            pch = "|", 
            layout = c(4, 1), 
            auto.key = list(columns = 3))

# ....

featurePlot(x = statdf[, 1:6], 
            y = statdf$file, 
            plot = "ellipse",
            ## Add a key at the top
            auto.key = list(columns = 1))

featurePlot(x = statdf[, 2:5], 
            y = statdf$file,
            plot = "density", 
            ## Pass in options to xyplot() to 
            ## make it prettier
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")), 
            adjust = 1.5, 
            pch = "|", 
            layout = c(4, 1), 
            auto.key = list(columns = 2))

# Boxplot
featurePlot(x = statdf[, 2:5], 
            y = statdf$file, 
            plot = "box", 
            ## Pass in options to bwplot() 
            scales = list(y = list(relation="free"),
                          x = list(rot = 90)),  
            layout = c(4,1 ), 
            auto.key = list(columns = 2))

theme1 <- trellis.par.get()
theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
theme1$plot.symbol$pch = 16
theme1$plot.line$col = rgb(1, 0, 0, .7)
theme1$plot.line$lwd <- 2
trellis.par.set(theme1)

featurePlot(x = statdf[, 2:5], 
            y = statdf$Mean, 
            plot = "scatter",
            type = c("p", "smooth"),
            span = 0.5,
            layout = c(4, 1))

# Creating Dummy variables
# The function dummyVars can be used to generate a 
# complete (less than full rank parameterized) #
# set of dummy variables from one or more factors
library(caret)
str(model.matrix(Genus ~., data = alluv_in))
head(dummies <- dummyVars(Genus ~ ., data = alluv_in)) # instead of model.matrix use it
str(predict(dummies, newdata = alluv_in))

# Note there is no intercept and each factor 
# has a dummy variable for each level, 
# so this parameterization may not be useful for some model functions, such as lm

# Model training and parameter tuning ----

# The train function can be used to :
# - evaluate, using resampling, the effect of model tuning parameters on performance
# - choose the “optimal” model across these parameters
# - estimate model performance from a training set

# method: manb, Model Averaged Naive Bayes Classifier

library(caret)
set.seed(998)
# statdf[, 2:5]
inTraining <- createDataPartition(statdf$Cycle, p = .75, list = FALSE)

training <- statdf[ inTraining,]
testing  <- statdf[-inTraining,]

# ---

fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)

gbmFit1 <- train(Cycle ~ ., data = training, 
                 method = "gbm", 
                 trControl = fitControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = FALSE)
gbmFit1
plot(gbmFit1)
# ---

gbmGrid <-  expand.grid(interaction.depth = c(1, 5, 9), 
                       n.trees = (1:30)*50, 
                       shrinkage = 0.1,
                       n.minobsinnode = 20)

nrow(gbmGrid)

set.seed(825)

gbmFit2 <- train(Cycle ~ ., data = training, 
                 method = "gbm", 
                 trControl = fitControl, 
                 verbose = FALSE, 
                 ## Now specify the exact models 
                 ## to evaluate:
                 tuneGrid = gbmGrid)

gbmFit2
# and visualize
trellis.par.set(caretTheme())
plot(gbmFit2) 
plot(gbmFit2, metric = "MAESD")

plot(gbmFit2, metric = "RsquaredSD", plotType = "level",
     scales = list(x = list(rot = 90)))

ggplot(gbmFit2)

# and ROC
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary,
                           search = "random")
# Error: wrong model type for regression
rda_fit <- train(Cycle ~ ., data = training, 
                 method = "rda",
                 metric = "ROC",
                 tuneLength = 30,
                 trControl = fitControl)
rda_fit
ggplot(rda_fit) + theme(legend.position = "top")

