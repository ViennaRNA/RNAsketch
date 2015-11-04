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
names<-cbind("x", "score", "mfe", "dE1", "dE2", "dE3");
data <- read.csv(opt$file, header=FALSE, sep = ";", dec = ".", comment.char='#', col.names=names);
# calculate mean, sd, se, mfe
cdata1 <- ddply(data, c("x"), summarise,
               struct = '1',
               score = mean(score),
               mfe = mean(mfe),
               dEmean = mean(dE1),
               dEsd   = sd(dE1),
               dEse   = dEsd / sqrt(length(dE1))
)
cdata2 <- ddply(data, c("x"), summarise,
               struct = '2',
               score = mean(score),
               mfe = mean(mfe),
               dEmean = mean(dE2),
               dEsd   = sd(dE2),
               dEse   = dEsd / sqrt(length(dE2))
)
cdata3 <- ddply(data, c("x"), summarise,
               struct = '3',
               score = mean(score),
               mfe = mean(mfe),
               dEmean = mean(dE3),
               dEsd   = sd(dE3),
               dEse   = dEsd / sqrt(length(dE3))
)
# merge frames
cdata<-rbind(cdata1, cdata2, cdata3)
# add the coordinates for the error bars
cdata <- ddply(cdata,.(x),transform,ystart = cumsum(dEmean),yend = cumsum(dEmean) + dEse)

# Standard deviation of the mean as error bar
pdf(paste(opt$file, ".pdf", sep=""))
p <- ggplot(cdata, aes(x=x, y=dEmean, fill=struct)) + geom_bar(stat="identity")

p + geom_segment(aes(xend=x,y=ystart,yend=yend), size = 0.1) + 
    geom_point(aes(x=x,y=yend), shape = "-", show_guide = FALSE, size = 1.5) +
    geom_point(aes(x=x,y=.02, colour=mfe), shape = 15, show_guide = FALSE, size = 3)
dev.off()
# plot done!
