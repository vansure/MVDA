setwd("~/Dropbox/paper tesi/bmc_article/Script_paper/MVDA/")
library("devtools")
remove.packages(pkgs = "MVDA")
remove.packages(pkgs="CorKohonen")
install.packages("~/Dropbox/paper tesi/bmc_article/Script_paper/CorKohonen/",type="source",repos=NULL)
install.packages("~/Dropbox/paper tesi/bmc_article/Script_paper/MVDA/",type="source",repos=NULL)

library(MVDA)
data(tcga_breast)

# 1 step miRNASeq Preprocessing
miRNA_km <- kmeans_preprocessing(DB = miRNA, nCenters = 10)
miRNA_km_summary <- clustering_summary(DB=miRNA, cluster=miRNA_km$clustering)

RNA_pamk <- pamk_preprocessing(DB = RNA,nCenters = 50)
RNA_pamk_summary <- clustering_summary(DB=RNA, cluster=RNA_pamk$clustering)

# 2 step prototype ranking
sda_km_miRNA <- sda_ranking(prototype = miRNA_km$centers,
							info = info$x,ranking.score = "avg")
plot(sda_km_miRNA)
sda_RNA <- sda_ranking(prototype = RNA_pamk$centers,
					   info = info$x,ranking.score = "avg")
plot(sda_RNA)

#3 step single view clustering
prototypes_miRNA <- miRNA_km$centers[rownames(sda_km_miRNA),]
nCenter <- 4
miRNA_pat_km <- kmeans_sv(nCenters = nCenter,prototype = prototypes_miRNA)
cm_km <- confusion_matrix(etichette = info$x,
	clustering = miRNA_pat_km$clustering,
	pazienti_geni = t(prototypes_miRNA),
	nCluster = nCenter)
miRNAerror <- CMsup_error(CM_sup = cm_km,nPat = dim(prototypes_miRNA)[2])
prototypes_RNA <- RNA_pamk$centers[rownames(sda_RNA),]
nCenter <- 4
RNA_pat_km <- kmeans_sv(nCenters = nCenter,prototype = prototypes_RNA)
cm_km <- confusion_matrix(etichette = info$x,
						  clustering = RNA_pat_km$clustering,
						  pazienti_geni = t(prototypes_RNA),
						  nCluster = nCenter)
RNAerror <- CMsup_error(CM_sup = cm_km,nPat = dim(prototypes_RNA)[2])

#4 step integration
nrows <- dim(info)[1]
ncols <- length(table(miRNA_pat_km$clustering))
rows_id <- rownames(info)
cols_id <- paste("Clustering",names(table(miRNA_pat_km$clustering)),sep=" ")
X1 <- supervised_matrix(etichette = miRNA_pat_km$clustering,nRow = nrows,
						nCol = ncols,row_id = rows_id,col_id  = cols_id)
X2 <- supervised_matrix(etichette = RNA_pat_km$clustering,nRow = nrows,
						nCol = ncols,row_id = rows_id,col_id  = cols_id)
SUP <- supervised_matrix(etichette = as.numeric(info$x),nRow = nrows,
						nCol = length(table(info$x)),
						row_id = rows_id,col_id  = cols_id)
X <- cbind(X1,X2,SUP)
K <- dim(X)[2]
patients <- t(cbind(t(prototypes_RNA),t(prototypes_miRNA)))
mf_res <- MF(A = X,k = K,eps = 0.01,iter_max = 100,nView = 2,V1.lc1 = 4,
			V1.lc2 = 4,V1.lc3 = 4,V1.lc4 = 0,V1.lc5 = 0,
			info = info,nclust_s = 4,pazienti = patients)
GLI_res <- general_late_integration(A = X,k = K,alfa = 1,
			eps = 0.001,pazienti = patients,info = info)



