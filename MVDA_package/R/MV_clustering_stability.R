#' A MVDA Function
#'
#' This function evaluate clustering stability by performing multi-view integration perturbing the dataset via leave-one-out and comparing
#' NMI between different clustering results.
#' 
#' @param method is the method you can use for the integration. Values can be matrix_factorization (default) of gli.
#' @param center is the matrix of centroid
#' @param class_cluster is the vector of class labels assigned to each cluster
#' @param X is the staked matrix of membeship matrix of each view.
#' @return a list containing the list of clustering obtained, the NxN matrix of NMI betwee each couple of clustering, the mean and the standard deviation of this matrix.
#' @export

mvclst <- function(method="matrix_factorization",center,class_cluster,X){
    if(method =="matrix_factorization"){
        X <- as.matrix(X)
        nRow <- dim(X)[1]
        cat("Removing one patient at a time...\n")
        nCol <- dim(X)[2]
        pb <- txtProgressBar(min=1,max=nCol,style=3)
        cluster_list <- list()
        nCluster <-dim(X)[1]
                
        for(i in 1:nCol){
            matrix_factorization(X=X[,-i],k1=nCluster,eps=1e-6,iter_max=10)->mf
            clustVector <- apply(mf$H,2,which.max)
            if(i > 1){
                pre <- clustVector[1:(i-1)]
            }else{
                pre <- NULL
            }
                    
            if(i < nCol){
                post <- clustVector[i:nCol]
            }else{
                post <- NULL
            }
            cluster_list[[i]]<- c(pre, NA, post)
            setTxtProgressBar(pb=pb,value=i)
        }
        close(pb)  
                
        nmi_matrix <- matrix(0, nCol, nCol)
        for(i in 1:nCol){
            for(j in 1:nCol){
                clus1 <- cluster_list[[i]]
                clus2 <- cluster_list[[j]]
                    nas <- union(which(is.na(clus1)), which(is.na(clus2)))
                    mutinfo <- mutinformation(clus1[-nas], clus2[-nas], method='emp')
                    entr1 <- entropy(clus1[-nas], method='emp')
                    entr2 <- entropy(clus2[-nas], method='emp')
                    nmi_matrix[i,j] <- mutinfo / mean(entr1, entr2)
            }
        } 
        return(list(cluster_list = cluster_list,nmi_matrix = nmi_matrix, sd=sd(nmi_matrix),mean=mean(nmi_matrix)));
    } #matrix factorization
    else{ #gli
        A <- as.matrix(X)
        nRow <- dim(A)[1]
        pb <- txtProgressBar(min=1,max=nRow,style=3)
        cluster_list <- list()
        k <-dim(A)[2]
        eps <- 0.1
        alfa <- 1
                
        for(i in 1:nRow){
            GLI(A=A[-i,],k=k,alfa=alfa,eps=eps)->GLI.res
            B <- GLI.res$B
            clustVector <- unlist(apply(B,1,which.max))
            if(i > 1){
               pre <- clustVector[1:(i-1)]
            }else{
                pre <- NULL
            }
                    
            if(i < nRow){
                post <- clustVector[i:nRow]
            }else{
                post <- NULL
            }
            cluster_list[[i]]<- c(pre, NA, post)
            setTxtProgressBar(pb=pb,value=i)
        }
        close(pb)            
                
        nmi_matrix <- matrix(0, nRow, nRow)
        for(i in 1:nRow){
            for(j in 1:nRow){
                clus1 <- cluster_list[[i]]
                clus2 <- cluster_list[[j]]
                        
                    nas <- union(which(is.na(clus1)), which(is.na(clus2)))
                    mutinfo <- mutinformation(clus1[-nas], clus2[-nas], method='emp')
                    entr1 <- entropy(clus1[-nas], method='emp')
                    entr2 <- entropy(clus2[-nas], method='emp')
                    nmi_matrix[i,j] <- mutinfo / mean(entr1, entr2)
            }
        }
        return(list(cluster_list = cluster_list,nmi_matrix = nmi_matrix, sd=sd(nmi_matrix),mean=mean(nmi_matrix)));  
    }#gli
}#end function

    
    