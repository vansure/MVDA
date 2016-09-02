#' A MVDA Function
#'
#' This function execute spectral clustering  on single view patient prorotypes
#' 
#' @param nCenters is the number of cluster we want to obtain
#' @param prototype is the matrix of prototype we want to cluster
#' @param kernel_function is the kernel function used in computing the affinity matrix. Default value: rbfdot
#' @keywords spectral-clustering;
#' @return a list containing three field: specc.res is the spectral clustering results. clustering is the vector with clustering assignment. center is the matrix with center prototypes.
#' @export

spectral_sv <- function(nCenters,prototype,kernel_function="rbfdot"){
  require(kernlab)
  if(is.null(prototype) || is.null(nCenters)){
    stop("You must insert prototypes and the number of desidered clusters!\n")
  }else{
    N_PAT <- dim(prototype)[2]
    cat("Executing spectral clustering...\n")
    specc(x=t(as.matrix(prototype)),centers=nCenters,kernel=kernel_function)->specc_res
    cluster.m <- specc_res@.Data
    center <- findCenter(DB=prototype,clust_vector=cluster.m)
    
    toRet <- list(specc.res=specc_res,clustering = cluster.m,center=center)
    return(toRet)
  }
}

