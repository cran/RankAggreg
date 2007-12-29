`kendall` <-
function(x, y, weights=NULL, weighted=TRUE)
{
    # Inputs:   x - lists to be combined
    #           y - candidate lists
    #           weights - weight matrix if weights to be used
    #           weighted - boolean if weights to be used
    # Outputs:  f.y - distances between each of candidate list and 
    #               x measured by Spearman formula
    
    k <- ncol(y)
    N <- nrow(x)

    # Kendall distances: regular and weighted
    kendall.dist <- function(x, y)
    {
        K=0
        n <- length(x)
        for (i in 1:(n-1))
            for (j in i:n)
                if((x[i] > x[j] & y[i] < y[j]) | (x[i] < x[j] & y[i] > y[j]))
                    K=K+1
        K
    }
    
    kendall.dist.w <- function(x, y, weight, k)
    {   
        K=0
        n <- length(x)
        for (i in 1:(n-1))
            for (j in (i+1):n)
                if((x[i] > x[j] & y[i] < y[j]) | (x[i] < x[j] & y[i] > y[j]))
                    K=K+abs(weight[ifelse(x[i]>k,k,x[i])]-weight[ifelse(y[j]>k,k,y[j])])
        K
    }


    if(weighted){   
        f.y <- apply(y, 1, function(z){
        sum <- 0
        for(i in 1:N){
            ul <- unique(c(z,x[i,]))
            rank.z <- match(ul,z,nomatch=k+1)
            rank.x <- match(ul,x[i,],nomatch=k+1)
            sum <- sum + kendall.dist.w(rank.x, rank.z, weights[i,], k)
        }
        sum
        })}
    else{
        f.y <- apply(y, 1, function(z){
            sum(apply(x,1,function(q){
                ul <- unique(c(z,q))
                rank.z <- match(ul,z,nomatch=k+1)
                rank.q <- match(ul,q,nomatch=k+1)
                kendall.dist(rank.z,rank.q)
            }))
        })
    }       
        
}

