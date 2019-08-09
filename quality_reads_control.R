
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
data_path <- '/Users/cigom/metagenomics/COI/run15/'
fileqc <- list.files(data_path, pattern = 'zip', full.names = TRUE)
fastq.r1 <- fileqc[1]
fastq.r2 <- fileqc[2]
source(qualMaxEEplot.R)
qualMaxEEplot(fastq.r1 = fastq.r1, fastq.r2 = fastq.r2)

# fastqc --nogroup yourFastqFile_R1.fastq yourFastqFile_R2.fastq
# From the fastQC output file:
# sed -n '/>>Per\sbase\ssequence\squality/,/>>END_MODULE/p' fastqc_data.txt  | sed '1d;$d' > fastq.qual.csv







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

