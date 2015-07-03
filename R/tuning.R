#############################################
#########	TUNING FUNCTION	#############
#############################################
# THIS FUNCTION DOES SPATIALLY-INDEPENDENT EVALUATIONS
# INPUT ARGUMENTS COME FROM WRAPPER FUNCTION

tuning <- function(occ, env, bg.coords, occ.grp, bg.grp, method, maxent.args, args.lab, categoricals, aggregation.factor, kfolds, bin.output, clamp) {

noccs <- nrow(occ)

# MAKE EVALUTION GROUPS
if (method == 'checkerboard1') group.data <- get.checkerboard1(occ, env, bg.coords, aggregation.factor)
if (method == 'checkerboard2') group.data <- get.checkerboard2(occ, env, bg.coords, aggregation.factor)
if (method == 'block') group.data <- get.block(occ, bg.coords)
if (method == 'jackknife') group.data <- get.jackknife(occ, bg.coords)
if (method == 'randomkfold') group.data <- get.randomkfold(occ, bg.coords, kfolds)
if (method == 'user') group.data <- get.user(occ.grp, bg.grp)

nk <- length(unique(group.data$occ.grp))
if(nk==1) stop("Partitioning method gave only 1 bin, cannot run evaluations.")
if(method %in% c('checkerboard2','block') & nk < 4) {
	warning(paste("Partitioning method gave only", nk, "bins.", sep=" "), call.=FALSE, immediate.=TRUE)
	}

# GET ENV VALUES AT OCC AND BG POINTS 
pres <-  as.data.frame(extract(env, occ))
bg <- as.data.frame(extract(env, bg.coords))

# DEFINE CATEGORICAL VARIABLES
if (!is.null(categoricals)) {
  for (i in 1:length(categoricals)) {
    pres[, categoricals[i]] <- as.factor(pres[, categoricals[i]])
    bg[, categoricals[i]] <- as.factor(bg[, categoricals[i]])    
  }
}

if (length(maxent.args) > 1) pb <- txtProgressBar(0, length(maxent.args), style=3)

full.mod <- list()
AUC.TEST <- data.frame()
AUC.DIFF <- data.frame()
OR10 <- data.frame()
ORmin <- data.frame()
predictive.maps <- stack()
nparam <- vector()
full.AUC <- vector()

for (a in 1:length(maxent.args)) {
	if (length(maxent.args) > 1) setTxtProgressBar(pb, a)

	# RUN THE MODEL WITH ALL POINTS
	x <- rbind(pres, bg)
	p <- c(rep(1, nrow(pres)), rep(0, nrow(bg)))
	tmpfolder <- tempfile()
	full.mod[a] <- maxent(x, p, args=maxent.args[[a]], factors=categoricals, path=tmpfolder)
	pred.args <- c("outputformat=raw", ifelse(clamp==TRUE, "doclamp=true", "doclamp=false")) 
	predictive.maps <- stack(predictive.maps, predict(full.mod[[a]], env, args=pred.args))
	full.AUC[a] <- full.mod[[a]]@results[5] # get training AUC from full model

	### RUN KFOLD EVALUATION
	for (k in 1:nk) {
		bin <- sort(unique(group.data$occ.grp))[k]
		train.val <- pres[group.data$occ.grp != bin,]
		test.val <- pres[group.data$occ.grp == bin,]
		bg.val <- bg[group.data$bg.grp != bin,]

	# RUN THE MODEL
	x <- rbind(train.val, bg.val)
	p <- c(rep(1, nrow(train.val)), rep(0, nrow(bg.val)))
	mod <- maxent(x, p, args=maxent.args[[a]], factors=categoricals, path=tmpfolder)

	# CALCULATE AUC.TEST AND AUC.DIFF
	AUC.TEST[a,k] <- evaluate(test.val, bg, mod)@auc
	AUC.DIFF[a,k] <- max(0, evaluate(train.val, bg, mod)@auc - AUC.TEST[a,k])

	# GET RAW PREDICTED SUITABILITY VALUES FOR TRAIN AND TEST
	p.train <- predict(mod, train.val, args=pred.args)
	p.test <- predict(mod, test.val, args=pred.args)

	# OMISSION RATE (OR) BASED ON "10% training omission" THRESHOLD
	if (nrow(train.val) < 10) {
		n90 <- floor(nrow(train.val) * 0.9) # if under 10, round down
	} else {
		n90 <- ceiling(nrow(train.val) * 0.9) # if over 10, round up
	}
	train.thr.10 <- rev(sort(p.train))[n90]
	OR10[a,k] <- mean(p.test < train.thr.10)

	# OMISSION RATE (OR) BASED ON "minimum training presence" THRESHOLD
	train.thr.min <- min(p.train)
	ORmin[a,k] <- mean(p.test < train.thr.min)
	}
	unlink(tmpfolder, recursive=TRUE)
}

# CALCULATE MEAN AND VARIANCE (CORRECTED FOR AUC) OF OPTIMALITY CRITERIA
bins <- sort(unique(group.data$occ.grp))

names(AUC.DIFF) <- paste("AUC.DIFF_bin", bins, sep=".")
Mean.AUC.DIFF <- rowMeans(AUC.DIFF)
Var.AUC.DIFF <- corrected.var(AUC.DIFF, nk)

names(AUC.TEST) <- paste("AUC_bin", bins, sep=".")
Mean.AUC <- rowMeans(AUC.TEST)
Var.AUC <- corrected.var(AUC.TEST, nk)

names(OR10) <- paste("OR10_bin", bins, sep=".")
Mean.OR10 <- rowMeans(OR10)
Var.OR10 <- corrected.var(OR10, nk)

names(ORmin) <- paste("ORmin_bin", bins, sep=".")
Mean.ORmin <- rowMeans(ORmin)
Var.ORmin <- corrected.var(ORmin, nk)

# DO AICc CALCULATIONS AND ADD TO "res" TABLE
for (i in 1:length(full.mod)) nparam[i] <- get.params(full.mod[[i]])
aicc <- calc.aicc(nparam, occ, predictive.maps)

# COMBINE INTO A SINGLE RESULTS TABLE "res"
features <- args.lab[[1]]
rm <- args.lab[[2]]
settings <- paste(args.lab[[1]], args.lab[[2]], sep='_')

res <- data.frame(settings, features, rm, full.AUC, Mean.AUC, Var.AUC, Mean.AUC.DIFF, Var.AUC.DIFF, Mean.OR10, Var.OR10, Mean.ORmin, Var.ORmin, aicc)

# GET EVALUATION STATS FOR EACH BIN SEPARATELY IF DESIRED
if (bin.output == TRUE) {
	res <- as.data.frame(cbind(res, AUC.TEST, AUC.DIFF, OR10, ORmin))
	}

names(predictive.maps) <- settings

results <- ENMevaluation(results=res, 
						predictions=predictive.maps, 
						partition.method=method,
						occ.pts=occ, 
						occ.grp=group.data[[1]], 
						bg.pts=bg.coords, 
						bg.grp=group.data[[2]])

return(results)
}
