`spearman` <-
function(x, y, importance, weights)
{   
    # Inputs:   x - lists to be combined
    #           y - candidate lists
    #           importance - the weight factors indicating the importance of ordered lists
    #           weights - weight matrix if weights to be used
    
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
    
    if(!is.null(weights)){   
        f.y <- apply(y, 1, function(z){
        sum <- 0
        for(i in 1:N){
            ul <- unique(c(z,x[i,]))
            rank.z <- match(ul,z,nomatch=k+1)
            rank.x <- match(ul,x[i,],nomatch=k+1)
            sum <- sum + importance[i]*subtrw(rank.x, rank.z, weights[i,], k)
        }
        sum
        })}
    else{
        f.y <- apply(y, 1, function(z){
            sum(importance*apply(x,1,function(q){
                ul <- unique(c(z,q))
                rank.z <- match(ul,z,nomatch=k+1)
                rank.q <- match(ul,q,nomatch=k+1)
                subtr(rank.z,rank.q)
            }))
        })
    }       
        
}

