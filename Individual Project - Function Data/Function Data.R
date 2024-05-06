#' ---
#' title: "Project Data Version: No 3"
#' author: "Khoirunnisa_Aulia"
#' 
#' output:
#'    html_document:
#'      toc: true
#'      toc_float: true
#'      highlight: zenburn
#' ---

#' Install the libraries 
library(ggplot2)
library(fda)

#############################################################################################################################

#' # Question a
#' Provide one composite visualisation of the raw data as connected lines over the year for every single observation, using the matplot function. Make sure the x-axis reveals the actual time scale.
data <- read.csv("data_version3.csv")
days <- as.numeric(sub("Day", "", colnames(data)))
transposed_data <- t(data)
matplot(days, transposed_data, type = 'l', lty = 1, xlab = "Day of the Year", ylab = "Enhance Vegetation Index (EVI)", main = "Composite Visualization of EVI over the Year", col = rainbow(nrow(data)))

#############################################################################################################################

#' # Question b 
#' Use a saturated B-spline basis, spanning over the whole year to fit the data using a standard roughness penalty and choose the penalty para mater using cross validation. Comment on the smoothness of the curves and provide graphs of GCV and the final smoothed data.

#' ## Define the order 
norder <- 4 #' Cubic B-splines (third-degree polynomials) are most commonly used because they provide a good balance between flexibility and smoothness.

#' ## Define the range of the observations for the whole year
days_range <- range(days)

#' ## Define the basis
nbasis <- length(days) + 4

#' ## Create the B-spline basis
bbasis <- create.bspline.basis(rangeval=days_range, nbasis=nbasis, norder=norder)

#' ## Set up the roughness penalty
init_lambda <- 1e6 #' The initial value of lambda  was set purposefully large so that we can distinguish a fit from the data.
curv.initialfdPar <- fdPar(bbasis, Lfdobj=int2Lfd(2), lambda=init_lambda) #' int2Lfd was used to specify a second-order penalty.

#' ## Fit the model using the fdSmooth function with the roughness penalty
EVI.initialSmooth <- smooth.basis(days, transposed_data, curv.initialfdPar)

#' ## Extract the fitted smooth curve:
names(EVI.initialSmooth)
fitted_values <- eval.fd(days, EVI.initialSmooth$fd)

#' ## Plot the raw data
matplot(days, transposed_data, type='l', col=rgb(0.5,0.5,0.5,0.5), xlab = "Day of the Year", ylab = "EVI", main = "EVI with Intial Fitted Smooth Curves")

#' ## Add the fitted smooth curve on top
matlines(days, fitted_values, col='red', lwd=2) #' The roughness penalty object (curv.fdPar) with an initial lambda, which is an arbitrary starting pointto give a baseline smoothing of the data, shows that it not the optimal smoothing because the lambda value has not been fine-tuned.

#' ## Perform Cross-Validation to Find Optimal Lambda
lambda_values <- 10^seq(-10, 0, length.out=100) #' lambda values to look over
gcv_values <- rep(NA, length(lambda_values)) #' define a variable to store gcv_values
for (i in seq_along(lambda_values)) {
  curv.fdPar <- fdPar(bbasis, Lfdobj=int2Lfd(2), lambda=lambda_values[i])
  smooth_fd_temp <- smooth.basis(days, transposed_data, curv.fdPar)
  gcv_values[i] <- smooth_fd_temp$gcv
}
optimal_lambda <- lambda_values[which.min(gcv_values)]

#' ## Fit the model using the optimal lambda
curv.finalfdPar <- fdPar(bbasis, Lfdobj=int2Lfd(2), lambda=optimal_lambda)
EVI.finalSmooth <- smooth.basis(days, transposed_data, curv.finalfdPar)

#' ## Plot the CGV Scores
plot(lambda_values, gcv_values, type='b', log='x', xlab='Log(Lambda)', ylab='GCV Score', main='GCV Scores for Different Lambda Values')

#' ## Plot the final smoothed data
fitted_values_optimal <- eval.fd(days, EVI.finalSmooth$fd)
matplot(days, transposed_data, type='l', col=rgb(0.5, 0.5, 0.5, 0.5), xlab="Day of the Year", ylab="EVI", main="EVI with Optimal Fitted Smooth Curves")
matlines(days, fitted_values_optimal, col='red', lwd=2)

#############################################################################################################################

#' # Question c
#'Adjust your code to use a harmonic acceleration penalty with a period of 1 year. Choose an appropriate penalty para mater using cross validation. Comment on the smoothness of the curves and provide graphs of GCV and the final smoothed data.

#' ## Define nbasis for a cubic line
nbasis <- length(unique(days)) + norder - 2 

#' ## Create the saturated B-spline basis
bbasis <- create.bspline.basis(rangeval=days_range, nbasis=nbasis, norder=norder)

#' ## Define the harmonic acceleration operator
harm_accel_Lfd <- vec2Lfd(c(0, (2*pi/365)^2), rangeval=days_range)

#' ## Perform cross-validation to find the optimal lambda for the harmonic acceleration penalty
lambda_values <- 10^seq(-10, 0, length.out=100)
gcv_values <- rep(NA, length(lambda_values))
for (i in seq_along(lambda_values)) {
  fdPar_object <- fdPar(bbasis, Lfdobj=harm_accel_Lfd, lambda=lambda_values[i])
  smooth_fd_temp <- smooth.basis(argvals=days, y=transposed_data, fdPar=fdPar_object)
  gcv_values[i] <- smooth_fd_temp$gcv
}
optimal_lambda <- lambda_values[which.min(gcv_values)]

#' ## Fit the final model with the optimal lambda
optimal_fdPar <- fdPar(bbasis, Lfdobj=harm_accel_Lfd, lambda=optimal_lambda)
EVISmoothFinal.harmonic <- smooth.basis(argvals=days, y=transposed_data, fdPar=optimal_fdPar)

#' ## Plot the CGV Scores
plot(lambda_values, gcv_values, type='b', log='x', xlab='Log(Lambda)', ylab='GCV Score', main='GCV Scores with Harmonic Acceleration Penalty')

#' ## Plot the final smoothed data
fitted_values_harmonic <- eval.fd(days, EVISmoothFinal.harmonic$fd)
matplot(days, transposed_data, type='l', col=rgb(0.5, 0.5, 0.5, 0.5), xlab="Day of the Year", ylab="EVI", main="EVI with Harmonic Optimal Fitted Smooth Curves")
matlines(days, fitted_values_harmonic, col='red', lwd=2)

#############################################################################################################################

#' # Question d
#' Based on the appropriate harmonic acceleration penalty fit calculate and plot the graphs of the first and second and derivatives of the curves.

#' ## Calculate the first derivative and store it
first_deriv <- deriv.fd(EVISmoothFinal.harmonic$fd)

#' ## Calculate the second derivative and store it
second_deriv <- deriv.fd(EVISmoothFinal.harmonic$fd, 2)

#' ## Plot the first derivative
plot(first_deriv, xlab="Day of the Year", ylab="First Derivative (Rate of Change)", main="First Derivative of EVI")

#' ## Plot the second derivative
plot(second_deriv, xlab="Day of the Year", ylab="Second Derivative (Acceleration)", main="Second Derivative of EVI")

#############################################################################################################################

#' # Question e
#' Conduct a unpenalized principal components analysis of these data. How many com- ponents do you need to recover 80% of the variation? Do the components appear satisfactory?
#' ## Conduct unpenalized PCA on the functional data
PCA_results <- pca.fd(EVISmoothFinal.harmonic$fd, nharm=3) #' Extract the first 3 principal components.

#' ## The cumulative percentage of variance explained by the principal components
cumulative_variance <- cumsum(PCA_results$varprop)
print(cumulative_variance)

#' ## Find the number of components needed to explain at least 80% of the variance
num_components <- which(cumulative_variance >= 0.8)[1]
cat("Number of components needed to explain at least 80% of variance:", num_components, "\n")

#' ## Plot the principal component functions (Each component)
par(mfrow=c(2, 2))
plot(PCA_results$harmonics[,1], main="First Principal Component")
plot(PCA_results$harmonics[,2], main="Second Principal Component")
plot(PCA_results$harmonics[,3], main="Third Principal Component")

#' ## Plot the principal component functions (full component)
plot(PCA_results$harmonics, xlab="Days", ylab="Function Value", main="Principal Components of the Unsmoothed Data")

#' ## Plot the scree plot
plot(PCA_results$varprop, xlab="Principal Component", ylab="Proportion of Variance Explained", type="b", main="Scree Plot")

#############################################################################################################################

#' # Question f
#' Penalized/Smoothed principal components analysis of the data and choice of the number of PCA’s
#' The smoothed functional data from (b), where the fitted the data with a saturated B-spline basis and the selected smoothing parameter via cross-validation, has been generated as a smoothed representation of the data. 
#' ## Conduct the PCA on the smoothed functional data
EVI.finalSmooth
PCA_results_penalized <- pca.fd(EVI.finalSmooth$fd, nharm=3)

#' ## Assess the variance explained by the PCA components
cum_var_penalized <- cumsum(PCA_results_penalized$varprop)
print(cum_var_penalized)

#' ## Find the number of components needed to explain a certain percentage of the variance
num_components_penalized <- which(cum_var_penalized >= 0.8)[1]
cat("Number of components needed to explain at least 80% of variance:", num_components, "\n")

#' ## Plot the Smoothed PCA’s with the optimal penalty
plot(PCA_results_penalized$harmonics, xlab="Days", ylab="Function Value", main="Principal Components of the Smoothed Data")

# Assuming 'lambda_values' and 'gcv_values' are already computed
plot(lambda_values, gcv_values, type='b', log='x', xlab='Log(Lambda)', ylab='GCV Score', main='GCV Scores for Different Lambda Values')

#############################################################################################################################

#' # Question h

# ## The Unpenalized PCA:
if(ncol(PCA_results$scores) >= 3) {
  cor_matrix_unpenalized <- cor(PCA_results$scores[, 1:3])
  print(cor_matrix_unpenalized)
}

#' ## The Penalized PCA:
if(ncol(PCA_results_penalized$scores) >= 3) {
  cor_matrix_penalized <- cor(PCA_results_penalized$scores[, 1:3])
  print(cor_matrix_penalized)
}

#' ## The codes to conduct a prod function on the Unpenalized PCA
#' Retrieve the basis information
bspline_PCA_result <- PCA_results$harmonics$basis
#' Create fd objects for each principal component
fd1 <- fd(PCA_results$harmonics$coefs[, "PC1"], bspline_PCA_result)
fd2 <- fd(PCA_results$harmonics$coefs[, "PC2"], bspline_PCA_result)
fd3 <- fd(PCA_results$harmonics$coefs[, "PC3"], bspline_PCA_result)
#' Calculate the inner products
inprod_unpenalized_12 <- inprod(fd1, fd2)
inprod_unpenalized_13 <- inprod(fd1, fd3)
inprod_unpenalized_23 <- inprod(fd2, fd3)
#' Print the results
print(inprod_unpenalized_12)
print(inprod_unpenalized_13)
print(inprod_unpenalized_23)

#' ## The codes to conduct a prod function on the Unpenalized PCA
#' Retrieve the basis information for penalized PCA results
bspline_PCA_penalized_result = PCA_results_penalized$harmonics$basis
#' Create fd objects for each principal component
fd1_penalized = fd(PCA_results_penalized$harmonics$coefs[, "PC1"], bspline_PCA_penalized_result)
fd2_penalized = fd(PCA_results_penalized$harmonics$coefs[, "PC2"], bspline_PCA_penalized_result)
fd3_penalized = fd(PCA_results_penalized$harmonics$coefs[, "PC3"], bspline_PCA_penalized_result)
#' Calculate the inner products
inprod_penalized_12 = inprod(fd1_penalized, fd2_penalized)
inprod_penalized_13 = inprod(fd1_penalized, fd3_penalized)
inprod_penalized_23 = inprod(fd2_penalized, fd3_penalized)
#' Print the results
print(inprod_penalized_12)
print(inprod_penalized_13)
print(inprod_penalized_23)

