### R code from vignette source 'RankAggreg.Rnw'

###################################################
### code chunk number 1: RankAggreg.Rnw:355-356
###################################################
    library(RankAggreg)


###################################################
### code chunk number 2: RankAggreg.Rnw:358-360 (eval = FALSE)
###################################################
##     help(package="RankAggreg")
##     vignette("RankAggreg")


###################################################
### code chunk number 3: RankAggreg.Rnw:381-384
###################################################
options(width=55)
library(xtable)
options(warn = -1)


###################################################
### code chunk number 4: RankAggreg.Rnw:387-398
###################################################
    options(warn = -1)
    library(clValid)
    library(mclust)
    library(kohonen)
    data(mouse)
    express <- mouse[1:100,c("M1","M2","M3","NC1","NC2","NC3")]
    rownames(express) <- mouse$ID[1:100]
    set.seed(100)
    result <- clValid(express, 5, 
        clMethods=c("hierarchical","fanny","model", "kmeans","sota","pam","clara",
        "agnes", "diana"), validation=c("internal","stability"))


###################################################
### code chunk number 5: RankAggreg.Rnw:407-435
###################################################
    getRanks_Weights <- function(clObj){
    temp <- clObj@measures[,1,]
    ranks <- matrix(0, nrow(temp), ncol(temp))
    colnames(ranks) <- 1:ncol(ranks)
    rownames(ranks) <- rownames(temp)
    weights <- ranks

    algs <- colnames(temp)
    measures <- rownames(temp)

    for(i in 1:nrow(temp)){
        if(measures[i] %in% c("Dunn", "Silhouette"))
            incr <- TRUE
        else
            incr <- FALSE   
        sorted <- sort(temp[i,], ind=TRUE, decreasing=incr)
        ranks[i,] <- algs[sorted$ix]
        weights[i,] <- sorted$x
    }
    list(ranks=ranks, weights=weights)
    }

    res <- getRanks_Weights(result)
    shortNames <- sort(c("ST", "FN", "KM", "PM", "HR", "AG", "CL", "DI", "MO"))
    longNames <- sort(unique(as.vector(res$ranks)))
    for(i in 1:length(shortNames))
        res$ranks[res$ranks==longNames[i]] <- shortNames[i]
    


###################################################
### code chunk number 6: RankAggreg.Rnw:440-441
###################################################
    print(res$ranks, quote=FALSE)


###################################################
### code chunk number 7: RankAggreg.Rnw:444-445
###################################################
     print(res$ranks, quote=FALSE)


###################################################
### code chunk number 8: RankAggreg.Rnw:456-457
###################################################
load("aggr_res.Rdata")


###################################################
### code chunk number 9: RankAggreg.Rnw:460-461 (eval = FALSE)
###################################################
## BruteAggreg(res$ranks, 9, res$weights, "Spearman")


###################################################
### code chunk number 10: RankAggreg.Rnw:469-470 (eval = FALSE)
###################################################
## (CEWS <- RankAggreg(res$ranks, 9, res$weights, seed=123))


###################################################
### code chunk number 11: RankAggreg.Rnw:473-474
###################################################
CEWS


###################################################
### code chunk number 12: figure1
###################################################
plot(CEWS)


###################################################
### code chunk number 13: fig1
###################################################
plot(CEWS)


###################################################
### code chunk number 14: RankAggreg.Rnw:503-504 (eval = FALSE)
###################################################
## (CEWK <- RankAggreg(res$ranks, 9, res$weights, "CE", "Kendall", seed=123, verbose=FALSE))


###################################################
### code chunk number 15: RankAggreg.Rnw:507-508
###################################################
CEWK


###################################################
### code chunk number 16: RankAggreg.Rnw:518-519 (eval = FALSE)
###################################################
## (GAWS <- RankAggreg(res$ranks, 9, res$weights, "GA", "Spearman",seed=123, verbose=FALSE))


###################################################
### code chunk number 17: RankAggreg.Rnw:521-522
###################################################
GAWS


###################################################
### code chunk number 18: RankAggreg.Rnw:525-526 (eval = FALSE)
###################################################
## (GAWK <- RankAggreg(res$ranks, 9, res$weights, "GA", "Kendall",seed=123, verbose=FALSE))


###################################################
### code chunk number 19: RankAggreg.Rnw:528-529
###################################################
GAWK


###################################################
### code chunk number 20: figure2
###################################################
plot(GAWS)


###################################################
### code chunk number 21: fig2
###################################################
plot(GAWS)


###################################################
### code chunk number 22: RankAggreg.Rnw:572-573
###################################################
data(geneLists)


###################################################
### code chunk number 23: RankAggreg.Rnw:576-577
###################################################
xtable(t(geneLists), caption="Top-25 upregulated genes from 5 prostate microarray experiments.")


###################################################
### code chunk number 24: RankAggreg.Rnw:587-588 (eval = FALSE)
###################################################
## top25CE <- RankAggreg(geneLists, 25, seed=100, rho=.01)


###################################################
### code chunk number 25: RankAggreg.Rnw:591-592
###################################################
top25CE


###################################################
### code chunk number 26: RankAggreg.Rnw:604-605 (eval = FALSE)
###################################################
## top25CEw <- RankAggreg(geneLists, 25, seed=100, importance=c(1,2,1,1,2), rho=.01)


###################################################
### code chunk number 27: RankAggreg.Rnw:608-609
###################################################
top25CEw


###################################################
### code chunk number 28: RankAggreg.Rnw:625-626 (eval = FALSE)
###################################################
## top25GA <- RankAggreg(geneLists, 25, seed=100, method="GA", maxIter=3000, convIn=50)


###################################################
### code chunk number 29: RankAggreg.Rnw:629-630
###################################################
top25GA


