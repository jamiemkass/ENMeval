########################################################################
### ENMeval Vignette

library(ENMeval)


########################################################################
########################################################################
########################################################################
### PREFACE ABOUT ALL ISSUES IN GETTING DATA READY...
### Get the data ready...
### This code comes directly from Hijmans and Elith (2015), Chapter 8, pg 43. 
files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), pattern='grd', full.names=TRUE)
predictors <- stack(files)
file <- paste(system.file(package="dismo"), "/ex/bradypus.csv", sep="")
bradypus <- read.table(file, header=TRUE, sep=',')
bradypus <- bradypus[,-1]
set.seed(0)
ext <- extent(-90, -32, -33, 23)
backg <- randomPoints(predictors, n=1000, ext=ext, extf = 1.25)
########################################################################
########################################################################
########################################################################



##### TABLE OF CONTENTS #####

######################################
######################################
###### DEMONSTRATION OF OPTIONS ######
######################################
######################################

######################################################
### 1) Partitioning methods ##########################
######################################################
occs <- as.data.frame(bradypus)
bgpts <- as.data.frame(backgr)
env <- predictors

# A run of ENMevaluate begins by using one of six methods to partition occurrence localities into testing and training bins (folds) for k-fold cross-validation (Fielding and Bell 1997; Peterson et al. 2011).
# Generally, the data partitioning step is done within the main 'ENMevaluate' function call.  In this section, we illustrate the different options here.
# The first, three partitioning methods are variations of what Radosavljevic and Anderson (2014) referred to as 'masked geographically structured' data partitioning.
# Basically, these methods partition both occurrence records and background points into evaluation bins based on some spatial rules.
# The methods are designed to reduce spatial-autocorrelation between points included in testing and training bins, which can overinflate model performance.

#######################
###### 1A) Block ######
#######################
# First, the 'block' method partitions data according to the latitude and longitude lines that divide the occurrence localities into four bins of (insofar as possible) equal numbers.
# Both occurrence and background localities are assigned to each of the four bins based on their position with respect to these lines. 

blocks <- get.block(occs, bgpts)

# The resulting object is a list of two vectors that supply the bin designation for each occurrence and background point.
str(blocks)

plot(!is.na(env[[1]]))
points(occs, pch=21, bg=blocks$occ.grp)

# ENMeval also includes two variants of a 'checkerboard' approach to partition occurrence localities. 
# These generate checkerboard grids across the study extent that partition the localities into bins. 
# In contrast to the block method, the checkerboard methods subdivide geographic space equally but do not ensure a balanced number of occurrence localities in each bin.
# For these methods, the user needs to provide a raster layer on which to base the underlying checkerboard pattern.
# Here we simply use the predictor variable RasterStack.
# Additionally, the user needs to define an 'aggregation.factor'.  
# This value tells the number of grids cells to aggregate when making the underlying checkerboard pattern.

###############################
###### 1B) Checkerboard1 ######
###############################
check1 <- get.checkerboard1(occs, env, bgpts, aggregation.factor=5)
plot(env[[1]])
points(occs, pch=21, bg=check1$occ.grp)

# The partitioning method is more clearly illustrated by looking at the background points:
points(bgpts, pch=21, bg=check1$bg.grp)

# We can change the aggregation factor to better illustrate how this partitioning method works:
check1.large <- get.checkerboard1(occs, env, bgpts, aggregation.factor=30)
plot(env[[1]])
points(bgpts, pch=21, bg=check1.large$bg.grp)
points(occs, pch=21, bg=check1.large$occ.grp, col='white', cex=1.5)

###############################
###### 1C) Checkerboard2 ######
###############################
# The second checkerboard method is partitions the data into k=4 bins.  
# This is done by aggregating the input raster at 2 scales.
# Points are assigned to bins with respect to where they fall in checkerboards of both scales.

check2 <- get.checkerboard2(occs, env, bgpts, aggregation.factor=c(5,5))
plot(env[[1]])
points(bgpts, pch=21, bg=check2$bg.grp)
points(occs, pch=21, bg=check2$occ.grp, col='white', cex=1.5)



##############################
###### 1D) User-defined ######
##############################
# Finally, for maximum flexibiliby, users can define a priori partitions.
# This provides a flexible way to conduct spatially independent cross validation with background masking.
# For example, perhaps we would like to partition points based on a k-means clustering routine.

ngrps <- 10
kmeans <- kmeans(occs, ngrps)
occ.grp <- kmeans$cluster

plot(env[[1]])
points(occs, pch=21, bg=occ.grp)

# When using the user-defined partitioning method, we need to supply ENMevaluate with group identifiers for both occurrence points AND background points.
# If we want to use all background points for each group, we can set the background to zero.

bg.grp <- rep(0, nrow(bgpts))
points(bgpts, pch=16, bg=bg.grp)

# Alternatively, we may think of various ways to partition background data. 
# This depends on the goals of the study but we might, for example, find it reasonable to partition background by clustering around the centroids of the occurrence clusters.

centers <- kmeans$center
d <- pointDistance(bgpts, centers, lonlat=T)
bg.grp <- apply(d, 1, function(x) which(x == min(x)))

plot(env[[1]])
points(bgpts, pch=21, bg=bg.grp)


###############################
###### 1E) k-1 Jackknife ######
###############################

# The last two partitioning methods do not partition the background points into different groups.
# Note that neither of these methods accounts for spatial autocorrelation between testing and training localities, which can inflate evaluation metrics, at least for data sets that result from biased sampling (Veloz 2009; Hijmans 2012; Wenger and Olden 2012).

# Primarily when working with small data sets (e.g. < ca. 25 localities), users may choose a special case of k-fold cross-validation where the number of bins (k) is equal to the number of occurrence localities (n) in the data set (Pearson et al. 2007; Shcheglovitova and Anderson 2013).
# This is referred to as the k-1 jackknife.

jack <- get.jackknife(occs, bgpts)
plot(env[[1]])
points(occs, pch=21, bg=jack$occ.grp)  # note that colors are repeated here...


###############################
###### 1F) Random k-fold ######
###############################
# Finally, the 'random k-fold' method partitions occurrence localities randomly into a userspecified number of (k) bins.
# This method is equivalent to the 'cross-validate' partitioning scheme available in the current version of the MAXENT software.

# For instance, let's partition the data into five evaluation bins:
random <- get.randomkfold(occs, bgpts, k=6)

plot(env[[1]])
points(occs, pch=21, bg=random$occ.grp)


################################################
###### CONCLUSION OF PARTITIONING SECTION ######
################################################
# Choosing among the data partitioning methods depends on the research objectives and the characteristics of the study system.
# We refer to XX for additional considerations on appropriate partitioning for evaluation.





####################################
####################################
###### 2) RUNNING ENMevaluate ######
####################################
####################################

# Once you decide which method of data partitioning you would like to use, you are ready to start building models.
# We move on to the main function in ENMeval: ENMevaluate.

# Initial considerations:
# Unless you supply the function with background points (which is recommended in many cases), you need to define how many background points should be used with the 'n.bg' argument.
# If any of your predictor variables are categorical values, you need to define which layer(s) in the 'categoricals' argument.


# The two main parameters to define when calling ENMevalaute are (1) the range of regularization multiplier values and (2) the combinations of feature class to consider.
# The regularization multiplier (RM) determines the penalty for adding parameters to the model.
# Higher RM values impose a stronger penalty on model complexity and thus result in simpler (flatter) model predictions.
# The feature classes determine the potential shape of the response curves.
# A model that is only allowed to include linear feature classes will be simpler than a model that is allowed to include all possible feature classes.
# Much more description of these parameters is available in (CITE OTHER RESOURCES).
# For the purposes of this vignette, we demonstrate simply how to adjust these parameters.
# The following section deals with comparing the outputs of each model.

# ENMevaluate builds a separate model for each unique combination of RM values and feature class combinations.
# For example, the following call will build and evaluate 2 models.
# One with RM=1 and one with RM=2, both allowing only linear features.

res1 <- ENMevaluate(occs, env, bgpts, method='checkerboard2', RMvalues=c(1,2), fc=c('L'))
res1

# We may, however, want to compare a wider range of models that can use a wider variety of feature classes:
res2 <- ENMevaluate(occs, env, bgpts, method='checkerboard2', RMvalues=c(1,2), fc=c('L','LQ','LQP'))
res2

# When buiding many models, the command may take a long time to run.
# Of course this depends on the size of your dataset and the computer you are using.
# When working on big projects, running the command in parallel can be faster:
res2.par <- ENMevaluate(occs, env, bgpts, method='checkerboard2', RMvalues=c(1,2), fc=c('L','LQ','LQP'), parallel=T)
res2.par

# Another way to save time at this stage is to turn off the option that generates model predictions across the full study extent (rasterPreds).
# Note that these are needed for calculating AICc values so those are returned as NaN when rasterPreds=F.
res3 <- ENMevaluate(occs, env, bgpts, method='checkerboard2', RMvalues=c(1,2), fc=c('L','LQ','LQP'), rasterPreds=F)
res3

# No predictions generated:
res3@predictions

# No AICc values calculated
res3@results$aicc


# We can also calculate one of two niche overlap statistics directly here using the 'niche.overlap' argument.
# Note that you can also calculate this value at a later stage using the separate 'calc.niche.overlap' function.

overlap <- calc.niche.overlap(res2@predictions, )
overlap

# The 'bin.output' argument determines if separate evaluation statistics for each testing bin are included in the results file.


#################################
#################################
###### 2) EXPLORING OUTPUT ######
#################################
#################################

#################################
### PLOTTING RESULTS
# Plot complex versus smooth models to see effect
# Talk about evaluation metrics (table from paper)
# Options for selecting optimal models
# Plotting response curves (Use response(ENMeval_results@models[[1]]))


#################################
### DOWNSTREAM
# Extracting model results from object (various threshold)
# Use model object to make a new prediction if you want a logistic prediction
# Make a projection to a new extent
# Do MESS map (Use mess() is dismo)



