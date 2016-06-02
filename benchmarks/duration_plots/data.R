#!/usr/bin/env Rscript

# Packages
check_package <- function(x) {
  if (!require(x,character.only = TRUE)) {
    install.packages(x,dep=TRUE)
  }
}

check_package("Hmisc")
check_package("gtools")
check_package("plyr")
check_package("ggplot2")
check_package("getopt")

# Option Parsing
spec = matrix(c(
    'directory', 'd', 1, "character",
    'help', 'h', 0, "logical"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

if (!is.null(opt$help)) {
cat(getopt(spec, usage=TRUE));
q(status=1);
}

if (!is.null(opt$directory)) {
setwd(opt$directory)
}
infiles <- mixedsort(list.files(pattern = '*.out'))

evaluate <- function(file) {
  data <- read.csv(file, header = TRUE, sep=";", dec = ".", comment.char='#')
  row_sum_diff <- c(rowSums(data[,grep("^diff_eos_mfe_", colnames(data))]))
  row_sum_prob <- c(rowSums(data[,grep("^prob_", colnames(data))]))
  result <-subset(data, select=c("num_mutations", "score", "mode"))
  result$number_of_structures <- as.factor(data$number_of_structures)
  result$diff <- row_sum_diff
  result$prob <- row_sum_prob
  summary(result)
  return(result)
}

all_infiles <- lapply(infiles, 
                      function(file) { 
                        tryCatch(evaluate(file), error = function(e) { stop(paste("Error in file ", file, ":", e$message)) }) 
                      }
)
all_infiles <- do.call(rbind,all_infiles)

head(all_infiles)

cdata <- ddply(all_infiles, c("num_mutations", "mode", "number_of_structures"), summarise,
               N = length(score),
               mean_score = mean(score),
               mean_diff = mean(diff),
               mean_prob = mean(prob),
               sd   = sd(score),
               se   = sd / sqrt(N),
               min = min(score),
               max = max(score)
)

max_score <- max(cdata$mean_score)
min_score <- min(cdata$mean_score)

cdata$rel_score <- with(cdata, (mean_score) / (max_score))
cdata$rel_sd <- with(cdata, (sd) / (max_score))
cdata$rel_min <- with(cdata, (min) / (max_score))
cdata$rel_max <- with(cdata, (max) / (max_score))

head(cdata)

svg(paste(opt$directory, ".score.svg", sep=""), width=6, height=4)
#ggplot(cdata, aes(num_mutations,rel_score, colour = mode, shape=mode)) + 
#stat_smooth() +
#geom_point() + 
qplot(num_mutations, rel_score, colour = mode, shape=mode, data = cdata) +
#geom_errorbar(data=cdata, mapping=aes(x=num_mutations, ymin = rel_min, ymax = rel_max), width=0.1) +
scale_x_log10(limits = c(0.2, 600000), breaks = c(1,10,100, 1000, 10000, 100000)) +
#coord_trans(x="log10") +
scale_y_continuous(limits = c(0, 1)) +
geom_line() +
ylab("relative mean Score") +
xlab("Number of Sampled Sequences") +
labs(colour="Sample Mode", shape="Sample Mode")
dev.off()

svg(paste(opt$directory, ".diff.svg", sep=""), width=6, height=4)
qplot(num_mutations, mean_diff, colour = mode, shape=mode, data = cdata) +
scale_x_log10(limits = c(0.2, 600000), breaks = c(1,10,100, 1000, 10000, 100000)) +
#coord_trans(x="log10") +
#scale_y_continuous(limits = c(0, 1)) +
geom_line() +
ylab("mean Diff EOS MFE") +
xlab("Number of Sampled Sequences") +
labs(colour="Sample Mode", shape="Sample Mode")
dev.off()

svg(paste(opt$directory, ".prob.svg", sep=""), width=6, height=4)
qplot(num_mutations, mean_prob, colour = mode, shape=mode, data = cdata) +
scale_x_log10(limits = c(0.2, 600000), breaks = c(1,10,100, 1000, 10000, 100000)) +
#coord_trans(x="log10") +
scale_y_continuous(limits = c(0, 1)) +
geom_line() +
ylab("mean Probability in Ensemble") +
xlab("Number of Sampled Sequences") +
labs(colour="Sample Mode", shape="Sample Mode")
dev.off()

