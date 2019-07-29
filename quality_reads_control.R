
#from plotQualityProfile function

fl <- '/Users/cigom/metagenomics/COI/run15/015-Mock27-COI-Zoo_S56_L001_R1_001.fastq.gz'
#fl <- list.files(data_path, pattern = 'fastq.gz', full.names = TRUE)
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
