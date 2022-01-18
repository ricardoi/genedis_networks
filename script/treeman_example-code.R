library(treeman)
#install.packages("remotes")
remotes::install_github("DomBennett/MoreTreeTools")
library(MoreTreeTools)

ntips <- 100
tree <- rtree(ntips)
dpsi <- 10 #power difference

#Generate Data
psi_1 <- dpsi
psi_2 <- dpsi
focal <- round(ntips*.05)
mean.incid <- ntips*.05

c1 <- genCommData(tree=tree, psi=psi_1, mean.incid= mean.incid,
                  mean.abun= mean.incid, nsites=1, focal=focal)
c2 <- genCommData(tree=tree, psi=psi_2, mean.incid= mean.incid,
                  mean.abun= mean.incid, nsites=1, focal=focal)

cmatrix <- rbind(c1,c2)
cmatrix[cmatrix > 0] <- 1
commplot(cmatrix, tree, groups=c(1,2), no.margin=FALSE)
mtext(text = paste0("psi = ", dpsi))

#Permutation test with treeman
tree_tm <- as(tree, 'TreeMan')
c1_ids <- colnames(c1)[c1[1, ] > 0]
c2_ids <- colnames(c2)[c2[1, ] > 0]

obs_ovrlp <- calcOvrlp(tree_tm, c1_ids, c2_ids)
iterations <- 99
null <- rep(NA, iterations)

  for(i in 1:iterations) {
          cat('....[', i, ']\n', sep= "")
          null_tips <- sample(tree_tm['tips'], length(c1_ids))
          null[i] <- calcOvrlp(tree_tm, c1_ids, null_tips)
        }
p_value <- sum(obs_ovrlp >= null)/iterations
hist(null, main="", xlab="", ylab="")
 abline(v=obs_ovrlp, col="red")
 mtext(paste0("P-value: ", signif(p_value, 3)))
 cat("P-value", signif(p_value, 3), "\n", sep="")
