#### Starts with line 80 of Vignette
results.table<-models@results #Isolate the results table 
results.table #Look at the tuning results. Here, you will find information about omission rates, AUC and AICc scores
######################################################################################################################
## You can access elements within the ENMevaluate object using the @ symbol. To get a look at your options, use
# simply run the object.
models 
# The first four lines tell us about the parameters used to during model tuning.
#######################################################################################################################
# To see a dataframe of evaluation statistics, run:
models@results
# Here, you can select optimal models using various criteria. For example, by lowest AICc score:
best.mod.AICc<- results.table[results.table$delta.AICc==0,]
#######################################################################################################################
# To access a RasterStack of the model outputs, run:
models@predictions
# Lets plot the optimal model according to AICc. Notice that this model is was the 10th one that we ran (row 10), which mean that
# it will be the 10th item in our RasterStack
plot(models@predictions[[10]])
#######################################################################################################################
# We can also access the list of model objexts. You can subset this list using
# double brackets: e.g. results@models[[2]]. This will also allow you to access various 
# elements of the model. This model object can also be used for predicting the model 
# into other time periods or geographic areas. The html file no longer exists for these, but
# other important information can be accessed.
# Lets take a look at the AICc optimal model object (item 10 again).
models@models[[10]]
models@models[[10]]@lambdas # lambda file to find which variables were used
models@models[[10]]@results # results shows the Maxent model statistics
#######################################################################################################################
models@occ.pts# This simply shows you the data.frame of occurrence records that you used.
#######################################################################################################################
models@occ.grp# This shows you which group each of you occurrence records was sorted into based on your grouping method. 
# To find which method you used: results@partition.method
#######################################################################################################################
models@bg.pts# This shows you the coordinates of the background points used.
#######################################################################################################################
models@bg.grp# This should you a vector of which group each background point was sorted into.
#######################################################################################################################
# If you used the argument to calculate raster overlap (overlap = T), you can access a matrix of all pairwise Schoener's 
# *D* comparisons (by default) between all predictions.
models@overlap
# You can also calculate niche overlap post using the ENMevaluate function. This can be done using the 
# calc.niche.overlap() function. Here, you can specify Moran's *I* or Schoener's *D*. 
calc.niche.overlap(models@predictions)
# You can also calculate this using a subset of rasters using ()



