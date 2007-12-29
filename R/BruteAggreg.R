`BruteAggreg` <-
function(x, k, weighted=FALSE, index.weights=NULL, distance=c("Spearman", "Kendall")){
    distance <- match.arg(distance, c("Spearman", "Kendall"))
    x <- x[,1:k]
    distinct <- apply(x, 1, function(y) ifelse(length(unique(y)) < k, 1, 0))
    if (sum(distinct) >= 1)
        stop("Elements of Each Row Must Be Unique")
    if (nrow(x)<2)
        stop("X must have more than 1 row") 
    
    if(weighted){
        index.weights <- index.weights[,1:k]
        #standardize weights:
        index.weights <- t(apply(index.weights,1,function(z) (z-min(z))/(max(z)-min(z))))
        if(index.weights[1,k]!=0)
            index.weights <- 1-index.weights
        if(dim(x)[1] != dim(index.weights)[1] || dim(x)[2] != dim(index.weights)[2])
        stop("Dimensions of x and weight matrices have to be the same")
    }

       
    comp.list <- unique(sort(as.vector(x)))
    n <- length(comp.list)

    if (k > n)
        stop("k must be smaller or equal to n") 

    library(gtools)
    perms <- permutations(n,k,comp.list)
    if(distance=="Spearman")
        f.y <- spearman(x, perms, index.weights, weighted)
    else
        f.y <- kendall(x, perms, index.weights, weighted)
    list(top.list=perms[which.min(f.y),],fn.score=min(f.y), distance=distance)
}

