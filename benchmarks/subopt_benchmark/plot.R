# Packages
check_package <- function(x) {
  if (!require(x,character.only = TRUE)) {
      install.packages(x,dep=TRUE)
  }
}

check_package("ggplot2")
check_package("getopt")
check_package("plyr")
check_package("reshape2")
check_package("corrplot")

# Option Parsing
spec = matrix(c(
    'file','f', 1, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

# if the input filename is missing, ask for it
if ( is.null(opt$file)) {
cat(getopt(spec, usage=TRUE));
q(status=1);
}

# Script
names<-cbind("length", "structures", "graph_construction", "num_cc", "max_special_ratio", "mean_special_ratio", "nos", "construction_time", "sample_time");
data <- read.csv(opt$file, header=FALSE, sep = ";", dec = ".", comment.char='#', col.names=names);

times <- c(grep("time$", names(data), value=TRUE))
data_sub <- subset(data, graph_construction!=0)
data_melt <- melt(data_sub, id.vars="mean_special_ratio", measure.vars=times)

number_ticks <- function(n) {function(limits) pretty(limits, n)}

# plot correlation plot
pdf(paste(opt$file, "_corrplot.pdf", sep=""))
cc <- cor(data_sub[, !names(data_sub)=="graph_construction"])
corrplot(cc, method = "color")
dev.off()

# plot times vs special ratio
pdf(paste(opt$file, "_time.pdf", sep=""))
p <- ggplot(data_melt, aes(x=mean_special_ratio, y=value, color=variable, shape=variable, group=variable))
p + geom_point() + 
    scale_y_continuous(limits = c(-0, 20), breaks=number_ticks(5)) + 
    geom_smooth(method=lm, fullrange=TRUE) + 
    scale_colour_discrete(name  ="Duration", labels=c("Construction", "Sampling")) + 
    scale_shape_discrete(name  ="Duration", labels=c("Construction", "Sampling")) + 
    xlab("mean Ratio Specials/Vertices") +
    ylab("Time [s]") + 
    ggtitle("Time vs Complexity")

dev.off()
# plot done!
