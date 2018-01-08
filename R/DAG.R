#' A DAG Function 
#'
#' This function allows you to learn the DAG structure from Gaussian data
#' @param X The data matrix  
#' @param K tuning parameter
#' @param tau tuning parameter
#' @param B A p by p matrix indicating nonzero elements as initial values
#' @param A0 A p by p matrix as initial values for A
#' @param opts.tol Tolerance for convergence
#' @param equal.sigma  Whether equal error variance is assumed. 
#' @return Estimated adjacency matrix and Estimated error variances for each node
#' @keywords DAG
#' @export
#' @useDynLib GDAG
#' @examples
#' m=50
#' amat <- matrix(0, m, m)
#' amat[2:m,1]<- 1 #(rbinom((m-1),1,0.5)-0.5)*runif((m-1),0.5,2)
#' sigma <- seq(1,0.5,length.out=m)
#' X <- rmvDAG(n,amat,sigma)
#' out<-DAG(X,10,0.01)

DAG <- function(X,K,tau,B=NULL,A0=NULL,opts.tol=1e-5,equal.sigma=FALSE ){
	# B is the initial nonzero pattern
	# A0 is the inital value for A
	m = ncol(X)
	n = nrow(X)
	#treatment for B. need more work
	if(is.null(B)){
		B=matrix(0,m,m)
		nonzero = 0
	} else{
		if(all.equal(dim(B),c(m,m)))   #need more conditions
			nonzero = sum(B==1)
			else stop("Invalid input B!")
	}
	#checking errors
	#treatment for A0
	if(is.null(A0)){
		A0=matrix(0,m,m)
	}else{
		if(!all.equal(dim(A0),c(m,m))) {
			stop("Invalid input B!")
		}			
	}
	sigma=rep(1,m)
	out=.C("DAG",as.double(X), A=as.double(A0),as.integer(m), as.integer(n), as.double(K),
		   as.double(tau), as.integer(B),as.integer(nonzero), sigma=as.double(sigma),as.double(opts.tol),
		   obj=as.double(0),as.integer(equal.sigma))
	
	Bout <- matrix(out[[7]],m,m)
	nonzero = sum(Bout==1)
	obj1 <- out$obj
    if(nonzero<=15 && nonzero >=2){
	    out2 = .C("DAG",as.double(X), A=as.double(A0),as.integer(m), as.integer(n), as.double(K),
		   as.double(tau), as.integer(t(Bout)),as.integer(nonzero), sigma=as.double(sigma),as.double(opts.tol),
		   obj=as.double(0),as.integer(equal.sigma))
		if(out2$obj < obj1) out=out2
	}
	
	Aout <- matrix(out$A,m,m)
	sigma <-out$sigma
	return(list(A=Aout,sigma = sigma))				
}