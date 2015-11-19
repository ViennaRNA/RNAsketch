# Packages
check_package <- function(x) {
  if (!require(x,character.only = TRUE)) {
      install.packages(x,dep=TRUE)
  }
}

check_package("ggplot2")
check_package("getopt")
check_package("plyr")

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
names<-cbind("N", "count", "length", "n1", "n2", "d1", "d2", "sequence");
data <- read.csv(opt$file, header=FALSE, sep = ";", dec = ".", comment.char='#', col.names=names);
data <- data[, !names(data)=="sequence"]
head(data)
data[, "delta1"] <- apply(data[, c("d1", "d2")], 1, min)
data[, "delta2"] <- apply(data[, c("d1", "d2")], 1, max)
head(data)

min_data <- as.list(apply(data, 2, min))
min_data
min_delta2 <- min_data$delta2




