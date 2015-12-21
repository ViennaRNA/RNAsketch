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
    'file','f', 1, "character",
    'x','x', 1, "character",
    'label','l', 1, "character",
    'title','t', 1, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

# if the input filename is missing, ask for it
if ( is.null(opt$file)) {
cat(getopt(spec, usage=TRUE));
q(status=1);
}

if ( is.null(opt$label)) {
opt$label = 'x';
}
if ( is.null(opt$title)) {
opt$title = '';
}
# Script
# jump;exit;mode;score;num_mutations;seq_length;sequence;graph_construction;num_cc;max_specials;max_component_vertices;
# max_special_ratio;mean_special_ratio;nos;construction_time;sample_time;
# mfe_reached_0;mfe_reached_1;mfe_reached_2;diff_eos_mfe_0;diff_eos_mfe_1;diff_eos_mfe_2;prob_0;prob_1;prob_2;
# names<-cbind("opt$x", "score", "mfe_reached", "diff_eos_mfe_1", "diff_eos_mfe_2", "diff_eos_mfe_3", "prob_1", "prob_2", "prob_3");
data <- read.csv(opt$file, header=TRUE, sep = ";", dec = ".", comment.char='#');
# rename exit to x for plotting
data <- eval(parse(text = paste0("rename(data, c(", opt$x, "='x'))")))
#summary(data)

all_mfe_reached <- data[,grep("^mfe_reached_", colnames(data))] # find all entries concerning mfe_reached_N
number_structures <- ncol(all_mfe_reached)
number_structures
cdata <- data.frame()
for (i in 1:number_structures) {
    colNameMFE<-paste("mfe_reached", toString(i-1), sep="_")
    print(colNameMFE)
    colNameDiff<- paste("diff_eos_mfe", toString(i-1), sep="_")
    colNameProb<- paste("prob", toString(i-1), sep="_")
    
    eval(parse(text = paste0(
    "currentcdata <- ddply(data, c('x'),\
    summarise, struct = \"", toString(i), "\", score = mean(score),\
    mfe = mean(", colNameMFE, ")*100, dEmean = mean(", colNameDiff, "),\
    dEsd   = sd(", colNameDiff, "), dEse   = dEsd / sqrt(length(", colNameDiff, ")),\
    pmean  = mean(", colNameProb, "))")))
    
    #currentcdata
    cdata <- rbind(cdata, currentcdata)
}
#cdata
# add the coordinates for the error bars
cdata <- ddply(cdata, c("x"),transform,ystart = cumsum(dEmean),yend = cumsum(dEmean) + dEse)
cdata
number_ticks <- function(n) {function(limits) pretty(limits, n)}

# Standard deviation of the mean as error bar
pdf(paste(opt$file, ".pdf", sep=""))
p <- ggplot(cdata, aes(x=x, y=dEmean, fill=struct)) + geom_bar(stat="identity") +
    scale_y_continuous(limits = c(-0.2, NA), breaks=number_ticks(5)) +
    expand_limits(y=5.5) +
    ylab(expression(paste("mean ", delta, "E [kcal]"))) +
    xlab(opt$label) +
    labs(fill="Target Structure") +
    ggtitle(opt$title)

p + geom_segment(aes(xend=x,y=ystart,yend=yend), size = 0.1) +
    geom_point(aes(x=x,y=yend), shape = "-", show_guide = FALSE, size = 1.5) +
    geom_point(aes(x=x,y=-0.1, colour=mfe), shape = 15, show_guide = FALSE, size = 3) +
    labs(colour='MFE reached [%]') #+
    #geom_point(aes(x=x,y=-0.2, colour=pmean), shape = 15, show_guide = FALSE, size = 3)
dev.off()
# plot done!
