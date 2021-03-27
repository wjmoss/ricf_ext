# INPUTS ---
# L, O: the 0-1 valued adjacency matrices indicating the structure of the graph G
# S: the p-by-p covariance matrix of observations from the true model
# Linit; Oinit: the initial values of the edge parameters for running the algorithm
# sigconv: Boolean (True = look for conv. in Sigma, False = look for conv. in L/O )
# tol: gives the max average entrywise absolute deviation in Sigma permited for convergence
# maxiter: gives the max number of iterations before divergence is accepted
# out: String: options: None/False, Final, All/True
# maxkap: maximum condition number accepted before error thrown
# B: (optional instead of L -- here B = t(L))
# iterreduction: Boolean (True = find stronly connected components and avoid iteration
#                         on nodes with no spouses and not in an scc,
#                         False = iteration on all nodes)
# acc: String, the method of acceleration
# K: length of window for Anderson accerelation
#
# OUTPUTS ---
# Sigmahat: the fitted value for Sigma resulting from the algorithm
# Bhat: the fitted value for the directed edge matrix (note: equal to I-t(Lambdahat))
# Omegahat, Lambdahat: the MLE for the parameter values resulting from the algorithm
# iterations: the number of iterations run by the algorithm before convergence or divergence accepted
# converged: TRUE or FALSE - based on whether or not the algorithm converged before maxiter reached
ricf_R_acc <- function(L = NULL, O, S, Linit = NULL, Oinit = NULL, sigconv=TRUE, tol=10^(-6),
                   maxiter=5000, out="none", maxkap = 1e13, B = NULL,
                   iterreduction=TRUE, acc="Anderson", K=5){
  if (is.null(L)) {
    if (!is.matrix(B))
      stop("B must be a matrix!")
    if (nrow(B) != ncol(B))
      stop("B must be a square matrix!")
    L = t(B)
  }
  if (!is.matrix(L))
    stop("L must be a matrix!")
  if (nrow(L) != ncol(L) || nrow(O) != nrow(L))
    stop("L must be a square matrix of the same size as O!")
  p <- nrow(L)
  
  # Initialize the directed edge parameters via OLS
  initL <- function(L, S) {
    Linit <- L
    parents <- apply(L, 2, function(x) which(x > 0))
    len <- apply(L, 2, sum)
    for (j in 1:p) {
      if (len[j] > 0) {
        for (k in 1:len[j]) {
          Linit[parents[[j]][k], j] <- S[parents[[j]][k], j] / S[parents[[j]][k], parents[[j]][k]]
        }
      }
    }
    return(Linit)
  }
  
  # Initialize the bidirected edge parameters at random
  # Initialize the bidirected edge parameters at random
  initO <- function(O) {
    ## O2 <- diag(diag(S))  ## for BAPs
    ## R <- diag(p)
    ## if(p>1)
    ##   for(i in 1:(p-1)) { for(j in (i+1):p) { R[j, i] <- R[i, j] <- rnorm(1) }}
    ## O2 <- R*O
    ## diag(O2) <- rep(0, p)
    ## diag(O2) <- rowSums(abs(O2)) + abs(rnorm(p))
    R <- diag(diag(S))
    s <- min(diag(R))
    if(p>1)
      for(i in 1:(p-1)) { for(j in (i+1):p) { R[j, i] <- R[i, j] <- runif(1,0,0.5)*s }}
    O2 <- R*(O!=0)
    diag(O2) <- diag(S)
    lll <- min(eigen(O2)$value)
    if(lll<0){
      diag(O2) <- diag(O2) + abs(lll) + 0.1
    }
    return(O2)
  }
  
  if (is.null(Linit)) { Linit <- initL(L, S) }
  if (is.null(Oinit)) { Oinit <- initO(O) }
  if (!is.matrix(Linit) || !is.matrix(Oinit))
    stop("O, S, Linit, and Oinit need to be matrices!")
  if (nrow(Linit) != ncol(Linit) || nrow(Oinit) != ncol(Oinit))
    stop("Linit and Oinit must be square matrices!")
  if (length(unique(c(nrow(Linit), nrow(Oinit), nrow(L), ncol(S)))) > 1)
    stop("One of the input matrices has the wrong dimension!")
  if (sum(abs(L * O)) != 0) {
    warning("The graph either contains a bow or a self-loop.",call. = FALSE, immediate. = TRUE)
  }
  if (maxiter <= 0 || maxiter %% 1 != 0)
    stop("A positive integer is needed for the max number of iterations!")
  if (tol <= 0)
    stop("A positive tolerance is needed for convergence to be possible!")
  if (!is.logical(sigconv))
    stop("sigconv needs to take on a logical value!")
  out = tolower(out)
  if (!is.character(out))
    stop("Output needs to be a string: none/false, final, or all/true!")
  if (!(out == "true" || out == "all" || out == "final" || out == "false" || out == "none")) {
    stop("Output variable needs to be: none/false, final, or all/true!")
  }
  
  Det <- function(Lcur) {
    return(det(diag(nrow(Lcur)) - Lcur))
  }
  
  ## Tarjan's algorithm, to find all strongly connected components
  tarjan <- function(u) {
    #print(u)
    dft[u] <<- low[u] <<- count
    count <<- count + 1
    visited[u] <<- 1
    Stack <<- append(Stack, u)
    for (v in which(L[u, ] != 0)) {
      if (! visited[v]) {
        tarjan(v)
        low[u] <<- min(low[u], low[v])
      } else if (v %in% Stack) {
        low[u] <<- min(low[u], dft[v])
      }
    }
    if (dft[u] == low[u]) {
      while(length(Stack) > 0) {
        if (Stack[[length(Stack)]] != u) {
          #print(Stack[[length(Stack)]])
          Stack[[length(Stack)]] <<- NULL
          visited[u] <<- visited[u] + 1
        } else {
          Stack[[length(Stack)]] <<- NULL
          #print(u)
          #print(" ")
          break
        }
      }
    }
  }
  
  
  iteronce <- rep(0, p)
  if (iterreduction && p>1){
    # determine the nodes that only need 1 iteration
    # has no siblings and not in a cycle
    visited = rep(0, p)
    dft = rep(0, p)
    low = rep(0, p)
    count = 1
    Stack = list()
    for (u in 1:p) {
      if (! visited[u]) {
        tarjan(u)
      }
    }
    for (i in 1:p){
      if (sum(O[i, ]) == 1 && dft[i] == low[i] && visited[i] == 1) iteronce[i] = 1
    }
    
  }
  
  #recording the nodes that need iterations, for anderson acceleration
  iternodes <- which(iteronce == 0)
  
  accklist <- list()
  accklist$x = list()
  accklist$t = list()
  
  
  Lcur <- Linit; Ocur <- Oinit
  iter <- 1
  acckind <- 1
  repeat {
    
    ## record values of x^(k)
    accnodeind <- 1
    for (node in iternodes){
      pa <- which(L[, node] != 0)
      sp <- which(O[, node] != 0)
      sp <- sp[sp != node]
      # an extra space for w_ii,-i (RSS_avg)
      if (iter == 1){
        # every node in iternodes corresponds to 2 (len(pa)+len(sp))*K matrices
        # store the K vectors (x^(k), T(x^(k)=x^(k+1))) of length (lenL+lenO) by K columns
        accklist$x[[accnodeind]] <- matrix(0, length(pa)+length(sp)+1, K)
        accklist$t[[accnodeind]] <- matrix(0, length(pa)+length(sp)+1, K)
      }
      
      ## record values of x^(k)
      accklist$x[[accnodeind]][1:(length(pa)+length(sp)), acckind] <- c(Lcur[pa, node], Ocur[sp, node])
      if (iter == 1){
        accklist$x[[accnodeind]][length(pa)+length(sp)+1, acckind] <- log(Ocur[node, node] - Ocur[node, -node] %*% solve(Ocur[-node, -node]) %*% Ocur[-node, node])
      }
      else {
        accklist$x[[accnodeind]][length(pa)+length(sp)+1, acckind] <- accklist$t[[accnodeind]][length(pa)+length(sp)+1, ((acckind-2) %% K + 1)]
      }
      accnodeind <- accnodeind + 1
    }
    
    if(p>1){
      
      accnodeind = 1
      for (i in 1:p) {
        if (iter > 1) {
          if (iterreduction && iteronce[i]) next
        }
        pa <- which(L[, i] != 0)
        n.pa <- length(pa)
        sp <- which(O[, i] != 0)
        sp <- sp[sp != i]
        n.sp <- length(sp)
        len <- n.pa + n.sp
        IB <- diag(p) - Lcur
        Elessi <- IB[, -i]
        
        if (kappa(Ocur[-i, -i]) > maxkap) {
          stop(paste("The condition number of Ocur[-i, -i] is too large for node", i))
        }
        Zlessi <- Elessi %*% solve(Ocur[-i, -i])
        Zhelp <- S %*% Zlessi
        # The following line gets the indices of Zlessi corresponding to spouses
        zsp <- c(sp[sp < i], sp[sp > i] - 1)
        
        if (len > 0) {
          if (n.pa > 0) {
            ## DETERMINE a AND a0
            ## a: coefficient of B[i,pa(i)]
            ## a0: const independent to B[i,pa(i)]
            a <- rep(0, n.pa)
            for (k in 1:n.pa) {
              temp <- Lcur
              temp[pa[k], i] <- 2
              det2 <- Det(temp)
              temp[pa[k], i] <- 1
              det1 <- Det(temp)
              a[k] <- det2 - det1
            }
            temp[pa, i] <- 0
            a0 <- Det(temp)
            
            
            ind.pos <- which(a != 0)
            pa.pos <- pa[ind.pos]
            ## general case
            if (length(ind.pos) > 0) {
              ## can be simplied if no siblings
              if (n.sp == 0){
                M <- S[pa, pa]
                m <- S[pa, i]
                alpha <- solve(S[pa, pa], S[pa, i])
                y0 <- S[i, i] - S[i, pa] %*% alpha
                coef <- alpha + as.numeric(y0 / (a0 + t(a) %*% alpha)) * solve(S[pa, pa], a)
                Lcur[pa, i] <- coef
              }
              ## n.sp > 0 case
              else{
                ZtZ <- (t(Zlessi) %*% Zhelp)[zsp, zsp, drop = F]
                ZtY <- t(Zhelp)[zsp, pa, drop = F]
                YtY <- S[pa, pa]
                #if ( ncol(ZtY)==0 ){
                #  stop(paste(ncol(ZtY), "not equal to" ,ncol(YtY)))
                #}
                M <- rbind(cbind(ZtZ, ZtY), cbind(t(ZtY), YtY))
                m <- rbind(t(Zhelp)[zsp, i, drop = F], S[pa, i, drop = F])
                alpha <- solve(M, m)
                y0 <- S[i, i] - t(alpha) %*% m
                coef <- alpha + as.numeric(y0 / (a0 + t(a) %*% alpha[(n.sp+1):len])) * solve(M, as.matrix(c(rep(0,n.sp), a)))
                Ocur[i, sp] <- Ocur[sp, i] <- coef[1:n.sp]
                Lcur[pa, i] <- if(n.pa > 0){coef[(n.sp + 1):(n.sp + n.pa)]} else{c()}
              }
            }
            
            ## a = 0, reduced to Linear Regression
            else {
              ## Perform linear regression for Betas and (possibly) Omegas
              if (n.sp > 0){
                ZtZ <- (t(Zlessi) %*% Zhelp)[zsp, zsp, drop = F]
                ZtY <- t(Zhelp)[zsp, pa, drop = F]
                YtY <- S[pa, pa]
                M <- rbind(cbind(ZtZ, ZtY), cbind(t(ZtY), YtY))
                m <- rbind(t(Zhelp)[zsp, i, drop = F], S[pa, i, drop = F])
                coef <- solve(M, m)
                y0 <- S[i, i] - t(coef) %*% m
                if (any(is.na(coef))) {
                  stop(paste("Collinearity observed for node", i))
                }
                Ocur[sp, i] <- Ocur[i, sp] <- coef[1:n.sp]
                Lcur[pa, i] <- coef[(n.sp + 1):len]
              }
              else{
                M <- S[pa, pa]
                m <- S[pa, i]
                coef <- solve(M, m)
                y0 <- S[i, i] - t(coef) %*% m
                Lcur[pa, i] <- coef[(n.sp + 1):len]
              }
              
            }
          }
          else {
            ## no parents, only siblings
            ## Solve a simplified version with no Betas (Only Omegas)
            M <- ZtZ <-  (t(Zlessi) %*% Zhelp)[zsp, zsp, drop = F]
            m <- t(Zhelp)[zsp, i, drop = F]
            coef <- solve(M, m)
            Ocur[sp, i] <- Ocur[i, sp] <- coef
            y0 <- S[i, i] - t(coef) %*% m
          }
          
          ## Find the variance omega_{ii}
          RSS_avg <- S[i, i] - 2 * t(as.matrix(coef)) %*% m + t(as.matrix(coef)) %*% M %*% as.matrix(coef)
          if (kappa(Ocur[-i, -i]) > maxkap) {
            stop(paste("The Condition number of Ocur[-i, -i] is too large for node", i))
          }
          
          if (acc == "Anderson" && iteronce[i] == 0){
            accklist$t[[accnodeind]][1:(length(pa)+length(sp)), acckind] <- c(Lcur[pa, i], Ocur[sp, i])
            accklist$t[[accnodeind]][length(pa)+length(sp)+1, acckind] <- RSS_avg 
            accnodeind <- accnodeind + 1
          }
          
          
          ### print ###
          
          Ocur[i,i] <- RSS_avg + Ocur[i, -i] %*% solve(Ocur[-i, -i]) %*% Ocur[-i, i]
          
        }
        else{
          Ocur[i, i] = as.matrix(S[i, i])
        }
        
        
      }
      
      ## compute Anderson acceleration when iter %% K == 0
      if (acc == "Anderson" && iter %% K == 0){
        templist <- lapply(1:length(iternodes), 
                           function(i) t(accklist$t[[i]]-accklist$x[[i]]) %*% 
                             (accklist$t[[i]]-accklist$x[[i]])
        )
        UtU <- purrr::reduce(templist, `+`)
        if (kappa(UtU) <= maxkap){
          e <- rep(1, K)
          c_K <- solve(UtU, e) / sum(solve(UtU))
          
          accnodeind <- 1
          for (node in iternodes){
            pa <- which(L[, node] != 0)
            sp <- which(O[, node] != 0)
            sp <- sp[sp != node]
            ## record values of T(x^(k))
            accklist$t[[accnodeind]][ ,acckind] <- accklist$x[[accnodeind]] %*% c_K
            Lcur[pa, node] <- accklist$t[[accnodeind]][1:length(pa), acckind]
            Ocur[sp, node] <- Ocur[node, sp] <- accklist$t[[accnodeind]][(length(pa)+1):(length(pa)+length(sp)), acckind]
            if (kappa(Ocur[-node, -node]) > maxkap) {
              stop(paste("The Condition number of Ocur[-node, -node] is too large for node", node))
            }
            Ocur[node, node] <- exp(accklist$t[[accnodeind]][length(pa)+length(sp)+1, acckind]) + 
              Ocur[node, -node] %*% solve(Ocur[-node, -node]) %*% Ocur[-node, node]
            accnodeind <- accnodeind + 1
          }
        }
      }
      
      
      sigcur <- t(solve(diag(p) - Lcur)) %*% Ocur %*% solve(diag(p) - Lcur)
      bhat <- diag(p)-t(Lcur)
      if (iter == 1) {
        if (out == "true" || out == "all"){cat(iter, "\n")}
        if (maxiter == 1) {
          if (out == "true" || out == "all" || out == "final") {
            cat("Sigmahat \n")
            print(sigcur)
            cat("\nBhat \n")
            print(bhat)
            cat("\nOmegahat \n")
            print(Ocur)
            cat("\nLambdahat \n")
            print(Lcur)
            cat("\niterations \n")
            print(iter)
          }
          break
        }
      }
      else if (iter > 1) {
        dsig <- mean(abs(sigcur - sigpast))
        dLO <- sum(abs(Lcur - Lpast) + abs(Ocur - Opast)) / (sum(L) + sum(O))
        if (out == "true" || out == "all") {
          dsig6 <- format(round(dsig, 6), nsmall = 6)
          dLO6 <- format(round(dLO, 6), nsmall = 6)
          cat(iter, " Avg Change in L & O: ", dLO6, "| Avg Change in Sigma: ", dsig6, "\n")
        }
        
        if ((sigconv && dsig < tol) || (!sigconv && dLO < tol) || iter >= maxiter) {
          if (out == "true" || out == "all" || out == "final") {
            cat("Sigmahat \n")
            print(sigcur)
            cat("\nBhat \n")
            print(bhat)
            cat("\nOmegahat \n")
            print(Ocur)
            cat("\nLambdahat \n")
            print(Lcur)
            cat("\niterations \n")
            print(iter)
          }
          break
        }
      }
      sigpast <- sigcur; Lpast <- Lcur; Opast <- Ocur
      acckind <- iter %% K + 1
      iter <- iter + 1
    }
    
    else {
      
      sigcur = as.matrix(S[p, p])
      bhat = as.matrix(0)
      Ocur = as.matrix(S[p, p])
      Lcur = as.matrix(0)
      iter = 1
      break
    }
    
  }
  return(list(Sigmahat = sigcur, Bhat = bhat, Omegahat = Ocur, Lambdahat = Lcur,
              iterations = iter, converged = (iter < maxiter), andersonlist = accklist, iternodes=iternodes))
}





# Optional method for calling ricf, using a mag object instead of L and O.
# Note: undirected edges not allowed here.
ricfmag <- function(mag, S, Linit = NULL, Oinit = NULL, sigconv=TRUE, tol=10^(-6), maxiter=5000, out="none", maxkap = 10^9)
{
  if (!is.matrix(mag) || !is.matrix(S))
    stop("mag and X must be matrices!")
  if (nrow(mag) != ncol(mag))
    stop("mag must be a square matrix!")
  L <- mag %% 100
  O <- mag %/% 100 + diag(nrow(mag))
  return(ricf(L,O,S,Linit,Oinit,sigconv,tol,maxiter,out,maxkap))
}

# example 2 p=5
L=matrix(c(0,0.5,0,0,0,
           0,0,0.3,0.5,0,
           0.2,0,0,0,0.4,
           0,0,0,0,0.4,
           0,0,0,0,0),5,5)
L=t(L)

O=diag(5)
O[2,5]=O[5,2]=0.2
S=t(solve(diag(5)-L)) %*% O %*% solve(diag(5)-L)
L=matrix(c(0,0,1,0,0,
           1,0,0,0,0,
           0,1,0,0,0,
           0,1,0,0,0,
           0,0,1,1,0),5,5)
O[2,5]=O[5,2]=1
ricf_R_acc(L=L, O=O, S=S, maxiter = 4, iterreduction = TRUE)
ricf_R_acc(L=L, O=O, S=S, iterreduction = FALSE)


# example 1, p=3
L=matrix(c(0,0,0.5,1,0,0,0,0.5,0),3,3)
L[3,1]=0
O=diag(3)
O[1,3]=O[3,1]=0.4
S=t(solve(diag(3)-L)) %*% O %*% solve(diag(3)-L)
L[2,3]=1
O[1,3]=O[3,1]=1
ricf_R(L=L, O, S)

# example 2 p=5
L=matrix(c(0,0.5,0,0,0,
           0,0,0.3,0.5,0,
           0.2,0,0,0,0.4,
           0,0,0,0,0.4,
           0,0,0,0,0),5,5)

O=diag(5)
O[2,5]=O[5,2]=0.2
S=t(solve(diag(5)-L)) %*% O %*% solve(diag(5)-L)
L=matrix(c(0,1,0,0,0,
           0,0,1,1,0,
           1,0,0,0,1,
           0,0,0,0,1,
           0,0,0,0,0),5,5)
O[2,5]=O[5,2]=1
ricf_R(L=L, O=O, S=S, maxiter = 3000, iterreduction = TRUE)
ricf_R(L=L, O=O, S=S, iterreduction = FALSE)
