#Require function findCenter in file funzioni.R
#Require kernlab package
#DB is the dataset where rows are feature to be clustered and nCenters is the 
#number of clusters we want in output
#you can specify kernel function, see ?specc

#' A MVDA Function
#'
#' This function execute spectral clustering  on single view patient prorotypes
#' 
#' @param DB is your matrix dataset
#' @param nCenters is the desidered number of cluster
#' @param kernel_f is the kernel function used in computing the affinity matrix. Default value: rbfdot
#' @keywords spectral-clustering;
#' @return a list containing three field: specc.res is the spectral clustering results. clustering is the vector with clustering assignment. center is the matrix with center prototypes.
#' @export

spectral_preprocessing<-function(DB=NULL,nCenters = NULL,kernel_f ="rbfdot"){
  require(kernlab)
  if(is.null(DB)){
    stop("You must insert a DB!\n")
  }else{
    if(is.null(nCenters)){
      stop("You must insert the number of centers!\n")
    }
    else{
      cat("Execute spectral clustering...\n")
      specc(x=DB,centers=nCenters,kernel=kernel_f)-> DB.specc
      clusters <- DB.specc@.Data
      centers <- DB.specc@centers
      toRet <- list(specc.result = DB.specc, clustering = clusters, centers = centers);
      cat("End...\n")
      return(toRet);
    }
  }
}



