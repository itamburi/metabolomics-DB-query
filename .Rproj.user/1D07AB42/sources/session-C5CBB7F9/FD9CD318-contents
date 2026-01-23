library(tidyverse)




mycpds = read.csv("mycpds.csv")

set.seed(123)  # for reproducibility
mycpds_short <- mycpds[sample(nrow(mycpds), 50), ]
row.names(mycpds_short) <- NULL
mycpds_short$X = NULL
write.csv(mycpds_short, "mycpds_short.csv")
