# Looking at lambda selection in huge

install.packages("Rcpp")
install.packages("RcppParallel")
# Install from .tar.gz file
install.packages(pkgs = "FastGGM.tar.gz", repos = NULL, type = "source")
# Or install from GitHub
library(devtools)
install_github("wt2015-github/FastGGM")

library(huge)

data_mat <- as.matrix(q2_data)

data_mat <- as.matrix(counts_z_transposed[Her2_cases,sig_deg_2k])

ug_glasso <- ug_glasso_her2

ug_glasso_her2 <- huge(x = data_mat, method = "glasso", nlambda = 3)
#ug_glasso <- huge(x = data_mat, method = "glasso", lambda = c(0.2, 0.3, 0.4, 0.5, 0.6), nlambda = 10)
plot(ug_glasso, align = TRUE)

ug_glasso$lambda
ug_glasso$sparsity


# Try the function using the covariance matrix as input
data_cov_mat <- cov(data_mat)
ug_glasso_cov <- huge(x = data_mat, method = "glasso", lambda = 0.5)

# Number of zero coeff along graph path icov
ug_glasso$df

#log liklihood
ug_glasso$loglik

# precision matrix
k <- ug_glasso$icov

# missing edges are when the precision matrix has entries that are zero, but it's symmetric so only need the top triangle
# to get the missing edges and not double count, we take the number of zeroes in the precision matrix and divide by 2

k_zeros <- sum(k[[2]] == 0)
missing_edges <- k_zeros/2


# Now we will explore what happens as we vary lambda

# define function to count missing edges

missing_edges_function <- function(i) {
  k <- ug_glasso$icov[i][[1]]
  k_zeros <- sum(k == 0)
  missing <- k_zeros/2
  missing
}


#model selection using ebic
out.select = huge.select(ug_glasso,criterion = "ebic")
plot(out.select)

ug_glasso2 <- huge(x = data_mat, method = "glasso", lambda = c(0.2, 0.3, 0.4, 0.5), nlambda = 10)
out.select2 = huge.select(ug_glasso2,criterion = "ebic")

ug_glasso3 <- huge(x = data_mat, method = "glasso", lambda = c(0.05, 0.1, 0.15, 0.2, 0.25), nlambda = 10)
out.select3 = huge.select(ug_glasso3,criterion = "ebic")

ug_glasso4 <- huge(x = data_mat, method = "glasso", lambda = c(0.0001, 0.001, 0.01, 0.05, 0.1), nlambda = 10)
out.select4 = huge.select(ug_glasso4,criterion = "ebic")

out.select4ric = huge.select(ug_glasso4,criterion = "ric")
out.select3ric = huge.select(ug_glasso3,criterion = "ric")
out.select2ric = huge.select(ug_glasso2,criterion = "ric")
out.select1ric = huge.select(ug_glasso,criterion = "ric")

out.select1stars = huge.select(ug_glasso,criterion = "stars")
out.select2stars = huge.select(ug_glasso2,criterion = "stars")
out.select3stars = huge.select(ug_glasso3,criterion = "stars")
out.select4stars = huge.select(ug_glasso4,criterion = "stars")

ug_glasso_select_ebic <- huge.select(ug_glasso,criterion = "ebic")
ug_glasso_select_stars = huge.select(ug_glasso,criterion = "stars")
