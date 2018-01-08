#' A rmvDAG Function 
#'
#' This function allows you to generate Gaussian data from a DAG
#' @param n The sample size  
#' @param amat The adjacency matrix of a DAG
#' @param sigma The error variance of each node  
#' @return Gaussian data with the given sample size
#' @keywords rmvDAG
#' @export
#' @examples
#' amat=matrix(c(0,1,0,0),2,2)
#' rmvDAG(50,amat)
rmvDAG <-function(n,amat,sigma=NULL){
    stopifnot(is.numeric(n), is.matrix(amat)) 
    p=dim(amat)[1]
    if(is.null(sigma)) sigma <- rep(1,p)
    errMat <- matrix(rnorm(n * p),nrow = n)
     X <- matrix(0, n, p)
    
    X[, 1] <- errMat[, 1] * sqrt(sigma[1])
    for (j in 2:p) {
        ij <- 1:(j - 1)
        X[, j] <- X[, ij, drop = FALSE] %*% amat[j,ij] + errMat[, j] * sqrt(sigma[j])
    }
    X
}