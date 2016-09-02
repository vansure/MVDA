#' A MVDA Function
#'
#' This function implements the matrix factorization integration process
#' @param A is the matrix resulting from the concatenation of the membership matrix for each view
#' @param k the desidered number of cluster
#' @param iter_mat the maximum allowed iteration
#' @param eps is the precision parameter
#' @param nView is the number of views
#' @param pazienti is your matrix dataset
#' @param info is the factor of class labels
#' @return a list containing the matrix factorization results, the clustering assignment, the confusion matrix, the error in confusion matrix, the matrix T of views influence, the matrix of clusters centroids, the labels for each cluster and the matrix of fisher test for each cluster
#' @export
#'


MF <- function(k,A,eps,iter_max,nView,V1.lc1,V1.lc2,
                                 V1.lc3,V1.lc4,V1.lc5,pazienti, info){
  matrix_factorization(X=t(A),k1=k,eps=eps,iter_max=iter_max)->mf
  cluster <- apply(mf$H,2,which.max)

  cm <- confusion_matrix(etichette=info$x,clustering=cluster,pazienti_geni=t(pazienti),nCluster=k)
  class_label(cm) -> ClassLabels
  error <- CMsup_error(cm,dim(A)[1])

  find_T(P=mf$P,k=k, nView=nView,V1.lc1=V1.lc1,V2.lc2=V1.lc2,V3.lc3=V1.lc3,V4.lc4 = V1.lc4,V5.lc5=V1.lc5) -> T_mat
  center <- findCenter(DB=t(pazienti),clust_vector=cluster)
  fisher_test(cluster,ClassLabels,cm,info$x) -> fisher_test_res

  toRet <- list(mf_res = mf,mv_cluster=cluster,confMat = cm, error = error, T_mat = T_mat, center=center,labels = ClassLabels, fisher=fisher_test_res)
  return(toRet)
}
