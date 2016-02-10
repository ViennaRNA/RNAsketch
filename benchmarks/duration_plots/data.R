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
               se   = sd / sqrt(N)
)
head(cdata)

pdf("score.pdf")
qplot(num_mutations, mean_score, colour = number_of_structures, shape=mode, data = cdata) +
scale_x_log10() + geom_line() +
ylab("mean Score") +
xlab("Number of Mutations") +
labs(colour="Number of Structures", shape="Sample Mode")
dev.off()

pdf("diff.pdf")
qplot(num_mutations, mean_diff, colour = number_of_structures, shape=mode, data = cdata) +
scale_x_log10() + geom_line() +
ylab("mean Diff EOS MFE") +
xlab("Number of Mutations") +
labs(colour="Number of Structures", shape="Sample Mode")
dev.off()

pdf("prob.pdf")
qplot(num_mutations, mean_prob, colour = number_of_structures, shape=mode, data = cdata) +
scale_x_log10() + geom_line() +
ylab("mean Probability in Ensemble") +
xlab("Number of Mutations") +
labs(colour="Number of Structures", shape="Sample Mode")
dev.off()

