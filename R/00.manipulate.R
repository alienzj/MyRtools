# transform data into package data

# library(devtools)
# library(roxygen2)
# library(dplyr)


# save package to DESCRIPTION
# use_package(package)

# common data1
# DATA <- read.csv(system.file("extdata", "fruits.csv", package="patternplot"))
# use_data(DATA)
# sessionInfo()

# common data2
# DATA2 <- mtcars
# use_data(DATA2)

# volcano_data
# read.table("inst/results.txt", header = T) -> DATA3
# use_data(DATA3)

# statistics data
# phen <- read.csv("inst/phenotype.csv") %>%
#       mutate(BMI=rnorm(150, mean = 25, sd=6),
#       Group=c(rep("ETB", 60), rep("ETP", 90)))
# prof <- read.table("inst/test.profile")
# use_data(phen, prof)
