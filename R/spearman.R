`spearman` <-
function(x, y, weights, weighted=TRUE)
{   
    # Inputs:   x - lists to be combined
    #           y - candidate lists
    #           weights - weight matrix if weights to be used
    #           weighted - boolean if weights to be used
    # Outputs:  f.y - distances between each of candidate list and 
    #               x measured by Spearman formula
    
    k <- ncol(y)
    N <- nrow(x)

    ### Spearman footrule functions (regular and weighted)
    subtr <- function(x,y)
        sum(abs(x-y))

    subtrw <- function(x,y,weight,k){
        # puts a weight of 0 for all x > k (weight[k]=0 due to normalizaton)
        sum(abs(weight[ifelse(x>k,k,x)]-weight[ifelse(y>k,k,y)])*abs(x-y))
        }
    ####
    
    if(weighted){   
        f.y <- apply(y, 1, function(z){
        sum <- 0
        for(i in 1:N){
            ul <- unique(c(z,x[i,]))
            rank.z <- match(ul,z,nomatch=k+1)
            rank.x <- match(ul,x[i,],nomatch=k+1)
            sum <- sum + subtrw(rank.x, rank.z, weights[i,], k)
        }
        sum
        })}
    else{
        f.y <- apply(y, 1, function(z){
            sum(apply(x,1,function(q){
                ul <- unique(c(z,q))
                rank.z <- match(ul,z,nomatch=k+1)
                rank.q <- match(ul,q,nomatch=k+1)
                subtr(rank.z,rank.q)
            }))
        })
    }       
        
}

