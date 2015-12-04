# Packages
check_package <- function(x) {
  if (!require(x,character.only = TRUE)) {
    install.packages(x,dep=TRUE)
  }
}

check_package("Hmisc")

infiles <- list.files(pattern = '*.out')

evaluate <- function(file) {
  if (file.info(file)$size == 0) # check if file ist empty
	{ print(file) }

  data <- read.csv(file, header = TRUE, sep=";", dec = ".", comment.char='#')
  data <- cbind(data,path=file) # add column with filename
  
  # sequence_length <- data[1,3] #sequence length: first entry of column 3, since all entries are the same
  sequence_length <- max(data$seq_length)
  mean_num_mutations <- mean(data$num_mutations)
  
  all_deltas <- data[,grep("^diff_eos_mfe_", colnames(data))] # find all entries concerning diff_eos
 
  # sort all deltas
#  sorted_deltas <- t(apply(all_deltas, 1, sort)) # all delta entries sorted row-wise for calculating mean and median
#  colnames(sorted_deltas) <- cbind(colnames(all_deltas))
#  rownames(sorted_deltas) <- rbind(rownames(all_deltas))
  
  d1 <- apply(all_deltas, 1, min) # find min delta of each row
  d2 <- apply(all_deltas, 1, max) # find max delta of each row
   
  # number of designed sequences with mfe
  n1 <- d1[d1 == 0]  
  sum_n1 <- length(n1) # length of vector, because n1 contains only entries where mfe is reached
  
  n2 <- d2[d2 == 0]
  sum_n2 <- length(n2)
 
  # find min d2 and min d1, if there are more than one min d2
  deltas <- data.frame(d1, d2) 
  min_d2 <-  deltas[deltas$d2 == min(deltas$d2),] # rows where column d2 is minimum
  min_d1 <- min(min_d2$d1) # since min_d2 contains all entries where d2 is min, it is sufficient to find one minimum (does not matter if there are more than one)
 
  d1_min <- min_d1 # not necessary
  d2_min <- min_d2[1,2] # all entries in this column mimima, therefore take first entry

  # mean, median, standard deviation and standard error of d1 and d2 
  mean_d1 <- mean(d1)
  median_d1 <- median(d1)
  sd_d1 <- sd(d1)
  se_d1 <- sd_d1/sqrt(length(d1))
  mean_d2 <- mean(d2)
  median_d2 <- median(d2)
  sd_d2 <- sd(d2)
  se_d2 <- sd_d2/sqrt(length(d2))
 
  # result values
  result <-data.frame(file, sequence_length, sum_n1, sum_n2, d1_min, d2_min, mean_d1, median_d1, mean_d2, median_d2, mean_num_mutations)

  return(result)
}

all_infiles <- lapply(infiles, evaluate)
all_infiles <- do.call(rbind,all_infiles)

colnames(all_infiles) <- c("RNA", "l", "n1", "n2", "d1", "d2",  "mean_d1", "median_d1", "mean_d2", "median_d2", "mean_num_mutations") # RNA = file

# calculate mean and median from all_infiles d1 and d2  
mean_all_d1 <- mean(all_infiles$d1)
median_all_d1 <- median(all_infiles$d1)
mean_all_d2 <- mean(all_infiles$d2)
median_all_d2 <- median(all_infiles$d2)

all_infiles <- rbind(all_infiles, data.frame(RNA="mean", l = NA, n1 = NA, n2= NA, d1= mean_all_d1, d2 = mean_all_d2, mean_d1= NA, median_d1= NA, mean_d2= NA, median_d2= NA, mean_num_mutations = NA))
all_infiles <- rbind(all_infiles, data.frame(RNA="median", l = NA, n1 = NA , n2= NA, d1= median_all_d1, d2 = median_all_d2, mean_d1= NA, median_d1= NA, mean_d2= NA, median_d2= NA, mean_num_mutations = NA))

# generate latex table
latex.default(all_infiles, cdec=c(0,0,0,0,2,2,2,2,2,2,2), na.blank=TRUE, booktabs=FALSE, table.env=FALSE, center="none", file="", title="")

