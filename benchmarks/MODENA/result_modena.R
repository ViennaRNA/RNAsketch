# Packages
check_package <- function(x) {
  if (!require(x,character.only = TRUE)) {
    install.packages(x,dep=TRUE)
  }
}

check_package("Hmisc")

infiles <- list.files(pattern = '*.out')

evaluate <- function(file) {
  data <- read.csv(file, header = TRUE, sep=";", dec = ".", comment.char='#')
  data <- cbind(data,path=file) # add column with filename
  
  sequence_length <- data[1,3] #sequence length: first entry of column 3, since all entries are the same
  sum_n1 <- sum(data$n1)
  sum_n2 <- sum(data$n2)
  
  # find min d2 and min d1, if there are more than one min d2
  min_d2 <-  data[data[,7] == min(data[,7]),] # rows where column 7(= d2) is minimum
  min_d1 <- min(min_d2[,6]) # since min_d2 contains all entries where d2 is min, it is sufficient to find one minimum (does not matter if there are more than one)

  d1_min <- min_d1 # not necessary
  d2_min <- min_d2[1,7] # all entries in this column mimima, therefore take first entry
 
  # mean, median, standard deviation and standard error of d1 and d2, not used 
  mean_d1 <- mean(data$d1)
  median_d1 <- median(data$d1)
  sd_d1 <- sd(data$d1)
  se_d1 <- sd_d1/sqrt(length(data$d1))
  mean_d2 <- mean(data$d2)
  median_d2 <- median(data$d2)
  sd_d2 <- sd(data$d2)
  se_d2 <- sd_d2/sqrt(length(data$d2))
  
  # result values
  result <-data.frame(file, sequence_length, sum_n1, sum_n2, d1_min, d2_min, mean_d1, median_d1, mean_d2, median_d2)
 
  return(result)
}

all_infiles <- lapply(infiles, evaluate)
all_infiles <- do.call(rbind,all_infiles)

colnames(all_infiles) <- c("RNA", "l", "n1", "n2", "d1", "d2",  "mean_d1", "median_d1", "mean_d2", "median_d2") # RNA = file

# calculate mean and median from all_infiles d1 and d2  
mean_all_d1 <- mean(all_infiles$d1)
median_all_d1 <- median(all_infiles$d1)
mean_all_d2 <- mean(all_infiles$d2)
median_all_d2 <- median(all_infiles$d2)

all_infiles <- rbind(all_infiles, data.frame(RNA="mean", l = NA, n1 = NA, n2= NA, d1= mean_all_d1, d2 = mean_all_d2, mean_d1= NA, median_d1= NA, mean_d2= NA, median_d2= NA))
all_infiles <- rbind(all_infiles, data.frame(RNA="median", l = NA, n1 = NA , n2= NA, d1= median_all_d1, d2 = median_all_d2, mean_d1= NA, median_d1= NA, mean_d2= NA, median_d2= NA))

# generate latex table
latex.default(all_infiles, cdec=c(0,0,0,0,2,2,2,2), na.blank=TRUE, booktabs=FALSE, table.env=FALSE, center="none", file="", title="")
