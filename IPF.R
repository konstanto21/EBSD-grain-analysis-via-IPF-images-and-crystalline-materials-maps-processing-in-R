################################################################################################################
# load libraries and create cluster of CPU cores for parallel computing
################################################################################################################

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

library (plyr); library(dplyr); library(tibble); library(future); plan(multisession)
library(sfsmisc); library(pracma); library(caret); library(ggplot2); library(magick); library(plotly)
library(gapminder); library(mosaic); library(COUNT); library(stats4); library(bbmle); library(tidyr)
library(doSNOW); library(readr); library(readxl); library(cluster)
#cl <- makeCluster(4, type="SOCK") # 4 - number of cores
#registerDoSNOW(cl) # Register Backend Cores for Parallel Computing
library("factoextra"); library(dbscan); library(fpc); library("cluster")
library(tidyverse); library(tibble); library(Rtsne); library("PerformanceAnalytics")
library(stringr); library(reshape2); library(RColorBrewer) # https://plotly.com/r/contour-plots/
library(parallel);library(e1071); library(doParallel); library(mclust); library(plotGMM)
library(imager); library(dplyr); library(ggplot2); library(tidyr); library(ggvoronoi); library(jpeg) #library(magick) library(tidyverse) library(magrittr)

# Get and Set the number of available cores
num_cores <- detectCores(); print(num_cores); num_cores <- 20

# Create a cluster object using the specified number of cores
cl <- makeCluster(num_cores)

# Register the cluster as the parallel backend
registerDoParallel(cl)

# Stop the cluster
#stopCluster(cl)

################################################################################################################ 
# *************** Data Preparation - Phase #1 ***************** #
################################################################################################################

mypath = ".../largemap" # select the path to the image file to be read

setwd(mypath)

# Read IPF image
img <- readJPEG("IPF.jpg") # IPF
imgDm <- dim(img)

# Assign RGB channels to data frame
img_wide <- data.frame(
  x = rep(1:imgDm[2], each = imgDm[1]),
  y = rep(imgDm[1]:1, imgDm[2]),
  R = as.vector(img[,,1]),
  G = as.vector(img[,,2]),
  B = as.vector(img[,,3])
)

# do the rescalling
# Coordinates fix: IPF 59.5 x 45.0 um & round digits to facilitate merging - make sure to know the exact dimensions of the image in real world scale - image reads the coordinates starting from 1, which should be starting from 0 for scaling reasons
img_wide$x <- img_wide$x - 1
img_wide$y <- img_wide$y - 1

img_wide$x <- img_wide$x*59.5/max(img_wide$x)
img_wide$y <- img_wide$y*45.0/max(img_wide$y)

img_IPF <- img_wide

# identify the exact same coordinate locations to combine multi-modal outputs

# correct coordinates - if needed
#img_wide <- filter(img_wide, x > 2.499)
#img_wide <- filter(img_wide, between(y, 12.529, (max(img_wide$y) - 10.099)))

# fix to start from (0,0) - if needed
#summary(img_wide)
#img_wide$x <- img_wide$x - min(img_wide$x)
#img_wide$y <- img_wide$y - min(img_wide$y)
#summary(img_wide)

# Round digits
img_wide$x <- round(img_wide$x, digits = 1)
img_wide$y <- round(img_wide$y, digits = 1)

img_wide <- img_wide %>%
  mutate(color = rgb(R, G, B))

head(img_wide)

# cluster separately: colours -> crystallographic planes

# ggplot theme to be used
plotTheme <- function() {
  theme(
    panel.background = element_rect(size = 3, colour = "black", fill = "white"),
    axis.ticks = element_line(size = 2),
    panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "gray90", linetype = "dashed"),
    axis.title.x = element_text(size = rel(1.2), face = "bold"),
    axis.title.y = element_text(size = rel(1.2), face = "bold"),
    plot.title = element_text(size = 20, face = "bold", vjust = 1.5)
  )
}

# see if we need to initialise to reproduce the same colours every time we cluster
# Combined map
set.seed(1000)
kClusters <- 50 # make sure to select the number of clusters which retains the level of detail that exists in the original image
kMeans <- kmeans(img_wide[, c("R", "G", "B")], centers = kClusters)
kColours <- rgb(kMeans$centers[kMeans$cluster,])
img_wide <- cbind(img_wide, kColours) # img_wide <- img_wide[,1:6]
unique(kColours) # length(unique(img_wide$color)) # img_wide_save <- img_wide

ggplot(data = img_wide, aes(x = x, y = y)) + 
  geom_point(colour = img_wide$kColours) +
  labs(title = paste("k-Means Clustering of", kClusters, "Colours")) +
  xlab("x") + ylab("y") + plotTheme()

# try to convert kColours to RGB values
df_test <- col2rgb(img_wide$kColours)
df_test <- t(as.data.frame(df_test))
df_test <- as.data.frame(df_test)
df_test$red <- df_test$red/255
df_test$green <- df_test$green/255
df_test$blue <- df_test$blue/255
head(df_test)
img_wide <- cbind(img_wide, df_test) # img_wide <- img_wide[,1:7]
head(img_wide)

col_df <- unique(img_wide$kColours)

merged_list <- lapply(col_df, function(col_val) filter(img_wide, kColours == col_val)) # check if colours are unique

# Define the custom plot theme function
#plotTheme <- function() {
#  theme_minimal() +
#    theme(axis.text.x = element_text(angle = 45, hjust = 1))}

# Use lapply to create a plot for each data frame in the merged_list
plot_list <- lapply(merged_list, function(df) {
  ggplot(data = df, aes(x = x, y = y)) + 
    geom_point(aes(colour = df$kColours)) + # aes(colour = kColours)
    xlab("x") + ylab("y") + 
    plotTheme()
})

################################################################################################################
################      ISOLATE ALL GRAINS      #########
################################################################################################################
# Create an empty list to store the processed dataframes
processed_dataframes <- list()

# Function to compare RGB values with specific colors - more generic colors (maybe pale colors will be used)
# Recommendation: choose colors for kColours as centers and compare with kmeans "kColours" since the original image contains a lot of shades, which compromises the annotation accuracy
get_color_category <- function(R, G, B) { # the kmeans colors are pale: this introduces a complication that favors pink pred
  # Define standard RGB values for each color category
  colors_data <- data.frame(
    color = c("green", "blue", "red", "yellow", "pink", "cyan", "white","black"),
    R = c(0.4666667, 0.3529412, 0.6705882, 0.7921569, 0.8470588, 0.4901961, 1, 0), # R = c(0.9, 0.9, 1, 1, 1, 0.9, 1, 0),
    G = c(0.8313725, 0.3490196, 0.2784314, 0.8196078, 0.4313725, 0.8313725, 1, 0), # G = c(1, 0.9, 0.9, 1, 0.7529, 1, 1, 0),
    B = c(0.4392157, 0.7843137, 0.2470588, 0.4000000, 0.6352941, 0.8156863, 1, 0) # B = c(0.9, 0.9, 0.9, 0, 0.7961, 1, 1, 0)
    #R = c(0.85, 0.85, 1, 1, 1, 0.85, 1, 0), # R = c(0.85, 0.85, 1, 1, 1, 0.85, 1, 0),
    #G = c(1, 0.85, 0.85, 1, 0.85, 1, 1, 0), # G = c(1, 0.85, 0.85, 1, 0.85, 1, 1, 0),
    #B = c(0.85, 0.85, 0.85, 0, 0.85, 1, 1, 0) # B = c(0.85, 0.85, 0.85, 0, 0.85, 1, 1, 0)
  )
  
  # Function to find the closest color category for each row (based on Euclidean distance)
  find_closest_color <- function(r, g, b) {
    abs_diff <- abs(colors_data$R - r) + abs(colors_data$G - g) + abs(colors_data$B - b)
    closest_color <- colors_data$color[which.min(abs_diff)]
    return(closest_color)
  }
  
  # Apply the function to each row of R, G, and B
  color_categories <- mapply(find_closest_color, R, G, B)
  
  return(color_categories)
}

#col2rgb("#77D470")/255 col2rgb("#5A59C8")/255 col2rgb("#AB473F")/255 col2rgb("#CAD166")/255 col2rgb("#D86EA2")/255 col2rgb("#7DD4D0")/255

# Loop through each dataframe in the plot_list and perform the operations
for (i in seq_along(plot_list)) {
  # Get the current dataframe from the list
  data <- as.data.frame(plot_list[[i]]$data)
  
  # Scaling
  scaled <- scale(data[, c(1, 2)])
  
  # k-means++ clustering implementation 
  
  # Choose an appropriate value of K (number of clusters) - here since it is not possible to automate the identification of grains, the number of clusters can be extremely large so that the excess of clusters is naturally populated with zero data points
  k_value <- 80
  
  # Use K-means++ initialization and run the algorithm in parallel with multiple starts (e.g., 100 starts)
  kmeans_model <- kmeans(scaled, centers = k_value, nstart = 100, algorithm = "Lloyd", trace = FALSE, iter.max = 300)
  
  # Get the cluster labels
  cluster_labels <- kmeans_model$cluster # unique(cluster_labels)
  
  data <- as_tibble(data)
  data$PDA <- 0
  data$PDA <- cluster_labels
  data$PDA <- as.numeric(data$PDA)
  data$PDA <- data$PDA + 80*(i-1)+1 # we have to change the annotation number in clusters: add a condition of $PDA <- $PDA + i*max(d$PDA) +1) or smthng like that!
  
  # Assuming you have the 'data' dataframe with columns 'R', 'G', 'B', and 'colors'
  # Perform the get_plane_annotation function and Apply the annotation and create a new column 'plane'
  data$plane <- with(data, get_color_category(red, green, blue)) # unique(data$plane) data$plane <- with(data, get_plane_annotation(R, G, B)) # unique(data$plane)
  
    # Filter out rows where the value appears less than 25 times
  value_counts <- table(data$PDA)
  threshold <- 5 # exclude clusters with very small population of data
  data <- data[data$PDA %in% names(value_counts[value_counts >= threshold]), ]

  # Save the processed dataframe in the list
  processed_dataframes[[i]] <- data
}

# Combine the processed dataframes into a single dataframe
final_dataframe <- do.call(rbind, processed_dataframes) #final_dataframe <- lapply(final_dataframe, as.numeric); final_dataframe <- as.data.frame(final_dataframe)
#final_dataframe$plane <- with(final_dataframe, get_plane_annotation(R, G, B))
# this I used for quick eval: final_dataframe$plane <- with(final_dataframe, get_color_category(red, green, blue)) # Apply the function to your data
unique(final_dataframe$plane) 

final_dataframe$plane <- ifelse(final_dataframe$plane == 'red', '(001)', 
                                ifelse(final_dataframe$plane ==  'blue', '(111)',
                                       ifelse(final_dataframe$plane == 'green', '(101)', 
                                              ifelse(final_dataframe$plane == 'cyan', '(101)U(111)', 
                                                     ifelse(final_dataframe$plane == 'yellow', '(001)U(101)', 
                                                            ifelse(final_dataframe$plane == 'pink', '(001)U(111)', 
                                                                   ifelse(final_dataframe$plane == 'white', 'White', 
                                                                          ifelse(final_dataframe$plane == 'black', 'Black', 'NA'
                                                                                 ))))))))

head(final_dataframe) # final_dataframe$PDA <- final_dataframe$PDA - 100
nrow(filter(final_dataframe, plane == "")) # count(filter(final_dataframe, plane == ""))
nrow(filter(final_dataframe, plane == "NA"))
str(final_dataframe)

nrow(filter(final_dataframe, plane == "(001)U(111)")) 
nrow(filter(final_dataframe, plane == "(001)"))
nrow(filter(final_dataframe, plane == "(111)"))
nrow(filter(final_dataframe, plane == "White"))
nrow(filter(final_dataframe, plane == "Black"))
nrow(filter(final_dataframe, plane == "(001)U(101)"))
nrow(filter(final_dataframe, plane == "(101)"))
nrow(filter(final_dataframe, plane == "(101)U(111)"))

# Now 'final_dataframe' contains the results of applying the operations to each specific dataframe in the list
write.csv(x = final_dataframe, file = 'IPF_grains_n_planes.csv')
 
################################################################################################################
# test how it looks like a grain
################################################################################################################
test <- filter(final_dataframe, PDA == 2423) 

ggplot(data = test, aes(x = x, y = y)) + 
  geom_point(aes(colour = kColours)) + # aes(colour = kColours)
  xlab("x") + ylab("y") + plotTheme()

################################################################################################################
# test how it looks all the specific grains allocated to a specific crystallographic plane
################################################################################################################
test2 <- filter(final_dataframe, plane == "(001)U(101)")
head(test2)

ggplot(data = test2, aes(x = x, y = y)) + 
  geom_point(colour = test2$kColours, size=1.5) +
  labs(title = paste("EBSD pattern")) +
  xlab("x") +  ylab("y") #+  plotTheme()

unique(test2$kColours) # "#B14A72" "#A17061" "#3F413F" "#AB473F" "#CA5A54"
#[1] "#C48B54" "#D6A670" "#CAD166" "#AFB851" "#B6D490" "#CAD4BB"
