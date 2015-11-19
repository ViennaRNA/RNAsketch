# Packages
check_package <- function(x) {
  if (!require(x,character.only = TRUE)) {
      install.packages(x,dep=TRUE)
  }
}

check_package("ggplot2")
# check_package("getopt")
check_package("plyr")
check_package("xtable")

infiles <- dir(pattern='\\.out$')

change.files <- function(file){
  data <- read.csv(file, header=TRUE, sep=";", dec = ".", comment.char='#')
  
}


lapply(infiles , change.files)



n1 <- data$n1
n2 <- data$n2

d1 <- data$d1
d2 <- data$d2

length_of_sequence <- mean(data$length.of.sequence)

## length of sequence weil gleich, besser designname oder so...
cdata1 <- ddply(data, c("length.of.sequence"), summarise,
               struct = '1',
               sum_n1 = sum(n1),
               sum_n2 = sum(n2),
                       
               mean_d1 = mean(d1),
               median_d1 = median(d1),
               sd_d1 = sd(d1),
               se_d1 = sd_d1/sqrt(length(d1)),
    
               mean_d2 = mean(d2),
               median_d2 = median(d2),
               sd_d2 = sd(d2),
               se_d2 = sd_d2/sqrt(length(d2))

)


sum_n1 <- sum(n1)
sum_n2 <- sum(n2)

# find min d2 und d1
mins_d2 <- which( d2 == min(d2)) #mins_d2 vector mit Positionen wo min values
if ((length(mins_d2) >= 2)){
    # smallest_d1 <- 135
    for (i in 1: length(mins_d2)){
        current <- d1[mins_d2[i]]
        if (current < smallest_d1){
            smallest_d1 <- current 
            index <- mins_d2[i]
        }
    }
    d1 <- d1[index]
    d2 <- d2[index]
}

result <- c(length_of_sequence,d1, d2, sum_n1, sum_n2)
print(result)

   
    





