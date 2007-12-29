`RankAggreg` <-
function(x, k, index.weights=NULL, use.weights=FALSE, method=c("CrossEntropy", "GeneticAlgorithm"), distance=c("Spearman", "Kendall"), 
            rho=.01, weight=.5, N=5*k*length(unique(sort(as.vector(x)))), error = .001, maxIter = 100,
            popSize=100, CP=1, MP=.001, informative=FALSE, v1=NULL, verbose=TRUE)
{
   
    method <- match.arg(method, c("CrossEntropy", "GeneticAlgorithm"))
    distance <- match.arg(distance, c("Spearman", "Kendall"))
    x <- x[,1:k]
    distinct <- apply(x, 1, function(y) ifelse(length(unique(y)) < k, 1, 0))
    if (sum(distinct) >= 1)
        stop("Elements of Each Row Must Be Unique")
    if (nrow(x)<2)
        stop("X must have more than 1 row") 

    compr.list <- unique(sort(as.vector(x)))
    n <- length(compr.list)
    
    if(method=="CrossEntropy"){
        comp.list <- 1:n
        x <- t(apply(x,1, function(xx) match(xx,compr.list)))}

    if(use.weights){
        index.weights <- index.weights[,1:k]
        #standardize weights:
        index.weights <- t(apply(index.weights,1,function(z) (z-min(z))/(max(z)-min(z))))
        if(index.weights[1,k]!=0)
            index.weights <- 1-index.weights
        if(dim(x)[1] != dim(index.weights)[1] || dim(x)[2] != dim(index.weights)[2])
        stop("Dimensions of x and weight matrices have to be the same")
    }

    if (k > n)
        stop("k must be smaller or equal to n") 

   
    if(method=="CrossEntropy"){
        v <- array(dim=c(n,k,maxIter))
        if(informative)
            v[,,1] <- as.matrix(v1)
        else
            v[,,1] <- 1/(n)
        y <- vector("numeric")
        diff <- vector("numeric")

        Nhat <- round(rho*N)
        if (Nhat < 5)
            stop("rho is too small")

        t <- 1
        repeat
        {
            cands <- mcmcProc(v[,,t], N, comp.list=comp.list)
            
            if(verbose)
                if(distance=="Spearman")
                    cat(" Calculating Spearman scores... \n")
                else
                    cat(" Calculating Kendall scores... \n")
                    
            if(use.weights)
                if(distance=="Spearman")
                    f.y <- spearman(x, cands, index.weights, weighted=TRUE)
                else
                    f.y <- kendall(x, cands, index.weights, weighted=TRUE)
            else
                if(distance=="Spearman")
                    f.y <- spearman(x, cands, weighted=FALSE)
                else
                    f.y <- kendall(x, cands, weighted=FALSE)
                    
            if(verbose)
                cat("Iteration", t, "Done! ")

            fy <- sort(f.y, ind=TRUE)
            y[t] <- fy$x[Nhat]
            good.cand <- cands[f.y <= y[t],]
            
            
            v[,,t+1] <- upd.prob(good.cand, v[,,t], weight, comp.list)
            
            best.cand <- compr.list[cands[fy$ix[1],]]
            rm(cands) # clean up
            
            y.l <- paste(best.cand, sep="", collapse=",")

            diff[t] <- sum(abs(v[,,t+1]-v[,,t]))/(n*k)
            
            if(verbose)     
                cat("\n",    c("Optimal value: ", y[t],
                        "\n Optimal List:  ",y.l, 
                       " \n Prob. Diff:    ",  diff[t]), 
                      ifelse(diff[t] < error,"<",">"), error, "\n\n")   

            if(diff[t] < error)          
                break 

            t <- t + 1
            if (t > 200){
            cat("Did not converge after 200 iterations. Please increase sample size N")
            break}
        }
    } else{
        #generate initial population randomly
        cands <- matrix(0, popSize, k)
        for(i in 1:popSize)
            cands[i,] <- sample(compr.list, k)
            
        #calculate obj. fn
        if(use.weights)
            if(distance=="Spearman")
                f.y <- spearman(x, cands, index.weights, weighted=TRUE)
            else
                f.y <- kendall(x, cands, index.weights, weighted=TRUE)
        else
            if(distance=="Spearman")
                f.y <- spearman(x, cands, weighted=FALSE)
            else
                f.y <- kendall(x, cands, weighted=FALSE)
        
        best.cand <- cands[which.min(f.y),]
        bestevery <- min(f.y)    
        
        conv=FALSE
        t <- 1
        iter <- 0
        while(!conv){
            #selection probability
            minf <- min(f.y)
            p.y <- (max(f.y)+1-f.y)/sum((max(f.y)+1-f.y))
            cpy <- cumsum(p.y)
            
            #select cands for the next generation
            ind <- runif(popSize)
            ind2 <- rep(0, popSize)
            for(i in 1:popSize)
                ind2[i] <- sum(ind[i] > cpy)+1
            cands <- cands[ind2,]

            # cross-over
            pairstocross <- floor(popSize*CP/2)
            samp <- sample(1:popSize, pairstocross*2)
            pointsofcross <- sample(2:k, pairstocross, replace=TRUE)
            for(i in 1:pairstocross){ 
                for(j in pointsofcross[i]:k){
                # this loop performs partially matched crossover (PMX) described in Section 10.5 of 
                # Data Mining: Concepts, Models, Methods, and Algorithms by Mehmed Kantardzic (2003)

                    t1 <- cands[samp[i],j]
                    t2 <- cands[samp[i+pairstocross],j]
                    
                    if(!is.na(t3 <- match(t2, cands[samp[i],])))
                        cands[samp[i], t3] <- t1
                    if(!is.na(t3 <- match(t1, cands[samp[i+pairstocross],])))
                        cands[samp[i+pairstocross], t3] <- t2
                                                
                    cands[samp[i], j] <- t2
                    cands[samp[i+pairstocross], j] <- t1
                }
            }              
                          
            # random mutations with probability MP
            mutations <- round(popSize*k*MP)
            rows <- sample(1:popSize, mutations, replace=TRUE)
            cols <- sample(1:k, mutations, replace=TRUE)
            for(i in 1:mutations)
                if(length(temp <- setdiff(compr.list, cands[rows[i],]))!=0)
                    cands[rows[i], cols[i]] <- sample(temp,1)
            
            
            #calculate obj. fn
            if(use.weights)
                if(distance=="Spearman")
                    f.y <- spearman(x, cands, index.weights, weighted=TRUE)
                else
                    f.y <- kendall(x, cands, index.weights, weighted=TRUE)
            else
                if(distance=="Spearman")
                    f.y <- spearman(x, cands, weighted=FALSE)
                else
                    f.y <- kendall(x, cands, weighted=FALSE)
  
            y.l <- paste(cands[which.min(f.y),], sep="", collapse=",")
            if(verbose)     
                cat("\n", "Iteration", t, ": ",  c("Optimal value: ", min(f.y),
                        "\n Optimal List:  ", y.l, "\n"))
            
            if(minf == min(f.y))
                iter <- iter+1
            else
                iter <- 1
            
            if(iter == 5 || t >= maxIter)
                conv=TRUE
            
            if(min(f.y) < bestevery){
                best.cand <- cands[which.min(f.y),]
                bestevery <- min(f.y)
            }
            t <- t+1
        }
    }
    list(top.list=best.cand, optimal.value=ifelse(method=="CrossEntropy", fy$x[1], bestevery), 
    sample.size = ifelse(method=="CrossEntropy", N, popSize), num.iter=t, method=method, distance=distance)
}

