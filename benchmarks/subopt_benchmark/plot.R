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



data_failed <- subset(data, graph_construction==0)
data_failed <- data_failed[, (names(data_failed) %in% c("length", "structures"))]
failed_count <- ddply(data_failed, .(structures, length), count)

# plot failed xy with correlation curve
pdf(paste(opt$file, "_failed.pdf", sep=""))
p <- ggplot(failed_count, aes(x=structures, y=length))
p + geom_point(aes(size=freq, colour=freq)) + 
    guides(color=guide_legend(), size = guide_legend()) + 
    geom_smooth(method=lm, fullrange=TRUE) + 
    xlab("Number of Structures") +
    ylab("Design Length") + 
    ggtitle("Graph Construction Timeouts (15s)")
dev.off()

times <- c(grep("time$", names(data), value=TRUE))
data_sub <- subset(data, graph_construction!=0)

cor.mtest <- function(mat, conf.level = 0.95) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    diag(lowCI.mat) <- diag(uppCI.mat) <- 1
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
            lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
            uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
        }
    }
    return(list(p.mat, lowCI.mat, uppCI.mat))
}

# plot correlation plot
pdf(paste(opt$file, "_corrplot.pdf", sep=""))
cc <- cor(data_sub[, !names(data_sub)=="graph_construction"])
res1 <- cor.mtest(cc, 0.95)
corrplot(cc, method = "color", order = "hclust", diag= FALSE, p.mat = res1[[1]], insig = "p-value")
dev.off()


number_ticks <- function(n) {function(limits) pretty(limits, n)}

data_melt <- melt(data_sub, id.vars="mean_special_ratio", measure.vars=times)
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
