
DiscML <- function (x, phy, CI = FALSE, model = "ER", const =NULL, reversible =FALSE,  kappa = 1,
                    p = TRUE , ip = 0.1, alpha=FALSE, ialpha= 0.5, gammak = 8, zerocorrection=FALSE,
                    rootprobability=FALSE, irootprobability="RANDOM", ivnum =2, delta = 1, automatrixadj= TRUE,
                    mu= FALSE, simplify = FALSE ) 
{
  plotalpha <- FALSE
  showprogress <- FALSE
  if(is.logical(irootprobability))
    stop("irootprobability cannot be set to a logical value")
  ismu <- mu  
  rateparameters <- p
  if(!identical(T, TRUE))
    stop("The object T is not set to a logical value of TRUE. Before using this function, please change the value of T to TRUE by typing T <- TRUE")
  if(!identical(F, FALSE))
    stop("The object F is not set to a logical value of FALSE. Before using this function, please change the value of F to FALSE by typing F <- FALSE")
  if(identical( model , "GTR"))
    reversible <- TRUE
  if(identical(reversible, TRUE))
    rootprobability<-TRUE  
  templevels = NULL
  type ="discrete"
  method = "ML"
  .getSEs <- function(out)
  {
    h <- out$hessian
    if (any(diag(h) == 0)) {
      warning("The likelihood gradient seems flat in at least one dimension (gradient null):\ncannot compute the standard-errors of the transition rates.\n")
      se <- rep(NaN, nrow(h))
    }
    else {
      se <- suppressWarnings(sqrt(diag(solve(h))))
    }
    se
  }
  
  if(is.matrix(model)|| toupper(model) =="GTR"|| toupper(model) == "ER" || toupper(model) == "BDBI"||toupper(model) =="ARD" || toupper(model) == "SYM"|| toupper(model) == "BDER" || toupper(model) == "BDSYM"||  toupper(model) =="BDARD"||  toupper(model) =="BDIER" || toupper(model) == "BDISYM"|| toupper(model) =="BDIARD")
  {}
  else{
    stop("You must input right 'model'. Please check your spelling")
  }  
  if(identical(alpha,FALSE) )
    talpha=1
  if(is.list(x))
  {
    
    vectorCheck <- logical()
    for(i in 1:length(x))
    {
      vectorCheck[i] <- is.vector(x[[i]])
    }    
    if(all(vectorCheck))
    {
      listToMatrix<- matrix(0,length(x),length(x[[1]]) )
      for(i in 1: length(x))
        listToMatrix[i,] <- x[[i]]
      x<- listToMatrix
    }
    else
      stop("Each element of a list shall be a vector")    
    TIPS=1:ncol(x[[1]])
    NROW=nrow(x[[1]])
    NCOL=ncol(x[[1]])    
  }  
  if(is.data.frame(x))
  {
    x <- suppressWarnings(matrix(as.matrix(x), nrow(x), ncol(x)))
    if(is.character(x))
    {   
      tempx <- x
      x<- suppressWarnings(matrix(as.integer(tempx), nrow(tempx), ncol(tempx)))
      x<- tempx[!is.na(rowSums(x)),] 
    }    
  }
  if(is.matrix(x))
  {
    if(is.character(x))
    {   
      tempx <- x
      x<- suppressWarnings(matrix(as.integer(tempx), nrow(tempx), ncol(tempx)))
      x<- tempx[!is.na(rowSums(x)),] 
    } 
    TIPS <-1:ncol(x)
    NROW <-nrow(x)
    NCOL <-ncol(x)
  }
  if(is.vector(x))
  {
    TIPS=1:length(x)
    NROW=1
    NCOL=length(x)
  }
  if(identical(simplify, TRUE))
  {
    if(is.matrix(x))
      x<- matrix( as.numeric(!(x==0)), nrow(x), ncol(x) )
    if(is.vector(x))
      x<-as.numeric(!(x==0))
  }
  phy2= phy
  phy2 <- reorder(phy, "cladewise")
  phy2$edge[order(phy2$edge[,2]),2][TIPS]=NA
  phy2$edge=phy2$edge[!is.na(phy2$edge[,2]),]
  likslist= list()
  zlikslist= list()
  if (!inherits(phy, "phylo")) 
    stop("object \"phy\" is not of class \"phylo\"")
  if (is.null(phy$edge.length)) 
    stop("tree has no branch lengths")
  type <- match.arg(type, c("continuous", "discrete"))
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  if (nb.node != nb.tip - 1) 
    stop("\"phy\" is not rooted AND fully dichotomous.")
  likelivector = numeric(0)
  if(is.list(x))
    listCount <- length(x)
  if(!is.list(x))
  {
    listCount <- 1
    xx <- x
  }  
  obj <- list()
  if (kappa != 1) 
    phy$edge.length <- phy$edge.length^kappa
  for( jj in 1:NROW )
  {
    if(is.matrix(xx))
    {
      x <-  xx[jj,]
    }
    if (length(x) != nb.tip) 
      stop("length of phenotypic and of phylogenetic data do not match.")
    if (!is.null(names(x))) {
      if (all(names(x) %in% phy$tip.label)) 
        x <- x[phy$tip.label]
      else
        warning("the names of 'x' and the tip labels of the tree do not match: the former were ignored in the analysis.")
    }  
    if(is.matrix(xx))
    {
      if (is.numeric(automatrixadj))
      {
        if(jj == 1)
        {
          automatrixadj <- sort( automatrixadj)
          check <- logical(length(xx))
        }        
        if(!all(x %in% automatrixadj)) 
          stop("auotomatixadj should contain all possible discrete integer characters.")
        if(jj ==1)
        {
          nl <- length(automatrixadj)
          if (!is.factor(automatrixadj)) 
            automatrixadj2<- factor(automatrixadj)
          lvls <- levels(automatrixadj2)
          tt<- 1:nl
          names(tt)<-lvls
        }
        if(0 %in% automatrixadj) 
          x<-tt[x]
        else
          x <- tt[x+1]
        names(x) <-NULL
      }
      if(identical(automatrixadj,TRUE))
      {
        if(jj==1)
        {
          if (!is.factor(xx)) 
            kk <- factor(xx)
          lvls <- levels(kk)
          nl <- length(lvls)
        }
        if (!is.factor(x)) 
          xtemp <- factor(c(x,xx))
        x <- as.integer(xtemp[1:length(x)])
      }
    }
    if(is.vector(xx))
    {
      if (is.numeric(automatrixadj))
      {
        automatrixadj = sort( automatrixadj)
        check = logical(length(xx))
        if(!all(x %in% automatrixadj)) 
          stop("automatixadj should contain all possible discrete integer characters.")
        nl <- length(automatrixadj)
        if (!is.factor(automatrixadj)) 
          automatrixadj2<- factor(automatrixadj)
        lvls <- levels(automatrixadj2)
        tt<- 1:nl
        names(tt)<-lvls
        if(0 %in% automatrixadj) 
          x<-tt[x]
        else
          x <- tt[x+1]
        names(x) <-NULL
      }
      if(identical(automatrixadj,TRUE))
      {
        if (!is.factor(x)) 
          x <- factor(x)
        nl <- nlevels(x)
        lvls <- levels(x)
        x <- as.integer(x)
      }
    }
    if(jj == 1)
    { 
      if( identical(zerocorrection,TRUE) && !any(as.integer(lvls)==0))
        stop("Remember: If zerrocorection is set to TRUE, unobservable character state should be denoted by '0'.
             Please make sure that character state 0 is included in all possible character states.")
      if (method != "ML") 
        stop("only ML estimation is possible for discrete characters.")
      if (any(phy$edge.length <= 0)) 
        stop("some branches have length zero or negative")
      
      if (!is.matrix(model)) {
        rate <- matrix(NA, nl, nl)
        switch(model, ER = np <- rate[] <- 1, ARD = {
          np <- nl * (nl - 1)
          rate[col(rate) != row(rate)] <- 1:np
        }, SYM = {
          np <- nl * (nl - 1)/2
          sel <- col(rate) < row(rate)
          rate[sel] <- 1:np
          rate <- t(rate)
          rate[sel] <- 1:np
        }, BDER = {
          np <- 1
          rate[row(rate)!=col(rate)] = np+1
          selg=row(rate)==(col(rate)-1)
          sell=row(rate)==(col(rate)+1)
          rate[selg]=1
          rate[sell]=1
          rate[is.na(rate)]=2
        }, BDARD = {
          
          np <- 2*(nl-1)
          rate[row(rate)!=col(rate)] = np+1
          
          for(i in 1:(nl-1))
          {
            rate[i, i+1]=i
            rate[i+1, i]= nl-1+i
          }
          rate[is.na(rate)]=np+1
          
        }, BDSYM = {
          np <- (nl-1)
          rate[row(rate)!=col(rate)] = np+1
          for(i in 1:(nl-1))
          {
            rate[i, i+1]=i
            rate[i+1, i]=i
          }
          rate[is.na(rate)]=np+1
          
        }, BDBI = {
          np <- 2
          rate[row(rate)!=col(rate)] = np+1
          for(i in 1:(nl-1))
          {
            rate[i, i+1]=1
            rate[i+1, i]=2
          }
          
          rate[is.na(rate)]=np+1
        }, GTR = {
          np <- nl * (nl - 1)/2
          sel <- col(rate) < row(rate)
          rate[sel] <- 1:np
          rate <- t(rate)
          rate[sel] <- 1:np        
        },  BDISYM = {
          np <- 2
          rate[row(rate)!=col(rate)] = np+1
          rate[1, 2] =1
          rate[2, 1] =1
          for(i in 2:(nl-1))
          {
            rate[i, i+1]=2
            rate[i+1, i]=2
          }
          rate[is.na(rate)]=np+1
          
        }, BDIER = {
          np <- 1
          rate[row(rate)!=col(rate)] = np+1
          selg=row(rate)==(col(rate)-1)
          sell=row(rate)==(col(rate)+1)
          rate[selg]=1
          rate[sell]=1
          rate[is.na(rate)]=2
          
        }, BDIARD = {
          np <- 4
          rate = matrix(np+1,nl,nl)
          rate[1, 2] =1
          rate[2, 1] =2
          for(i in 2:(nl-1))
          {
            rate[i, i+1]=3
            rate[i+1, i]=4
          }
          rate[is.na(rate)]=np+1    
        })
      }
      else {
        rate = matrix(NA,nl,nl)
        if (ncol(model) != nrow(model)) 
          stop("the matrix given as 'model' is not square")
        if (ncol(model) != nl) 
          stop("the matrix 'model' must have as many rows as the number of categories in 'x'")
        rate[] <- suppressWarnings(as.integer(model))
        rate[is.na(rate)] <- 0
        
        np <- sum(suppressWarnings(!is.na( as.integer( model[col(model)!= row(model)] ))))
        diag(rate) <- np+1
      }
      if( identical(reversible,TRUE) || model =="GTR")
      {
        tmodel <- suppressWarnings(as.integer(model))
        ttmodel <- (tmodel == t(tmodel))
        test <- any(!ttmodel)
        if(is.matrix(model))
        {
          if(any(is.na(test)))
            stop("If 'reversible' is set to be TRUE, you must use symmetric metrix as a model. It is because when 'reversible' is set to be TRUE,
                 DiscML converts Symmetric matrix into reversible matrix. (Please see description documents).")
          if(!any(test))
            stop("If 'reversible' is set to be TRUE, you must use symmetric metrix as a model. It is because when 'reversible' is set to be TRUE,
                 DiscML converts Symmetric matrix into reversible matrix. (Please see description documents).")
        }
        else if ( model == "ARD"|| model =="BDARD" || model =="BDIARD" )
          stop("If 'reversible' is set to be TRUE, you must use symmetric metrix as a model. It is because when 'reversible' is set to be TRUE,
               DiscML converts Symmetric matrix into reversible matrix. (Please see description documents).")
        }
      index.matrix <- rate
      tmp <- cbind(1:nl, 1:nl)
      index.matrix[tmp] <- NA
      rate[tmp] <-  np + 1
      TIPS <- 1:nb.tip
      } 
    liks <- matrix(0, nb.tip + nb.node, nl)
    liks[cbind(TIPS, x)] <- 1
    likslist[[jj]] <- liks
  }
  phy <- reorder(phy, "postorder")
  Q <- matrix(0, nl, nl)
  e1 <- phy$edge[, 1]
  e2 <- phy$edge[, 2]
  EL <- phy$edge.length
  zliks <- matrix(0, NROW, nl)
  if(identical(rateparameters,FALSE) )
    stop("rateparameters cannot be set to FALSE. Hosever,they can be set to TRUE's or numeric values")
  if(is.numeric(rootprobability) )
  {
    if(length(rootprobability)!= nl)
    {
      
      stop("The number of elements of the root probability you have imputed does not match the number of possible character states, which is ", nl)
      
    }
  }
  
  if(is.numeric(irootprobability)&&identical(rootprobability,TRUE))
  {
    if(length(irootprobability)!=nl)
      stop("The number of elements in 'irootprobability' does not match the number of possible character states, whiche is ", nl, "\n")
  }
  dev <- function(r, output.liks = FALSE, output.zliks= FALSE, output.div=FALSE, turnRootProbabilityToFalse = FALSE, turnAlphaToFalse=FALSE, manualExceptAlpha=FALSE, output.Q=FALSE) {
    if(turnAlphaToFalse)
      alpha <-FALSE
    if(turnRootProbabilityToFalse)
      rootprobability <- FALSE
    if (any(is.nan(r)) || any(is.infinite(r))) 
      return(1e+50)
    if(!identical(manualExceptAlpha,FALSE))
    {
      length <- 0
      if(identical(rateparameters,TRUE))
      {
        p <- manualExceptAlpha[(length+1):(length+np)]
        length <- length + np
      }
      if(!identical(rateparameters,TRUE))
      {
        p <- rateparameters
        length <- length + np
      }
      if(identical(rootprobability, TRUE)) 
      {
        rprob <- manualExceptAlpha[(length+1):(length+nl)]
        if(sum(rprob)<=0)
          return(1e50)
        else
          rprob <- rprob/sum(rprob)
        length <- length + nl
      }
      else if(identical(rootprobability, FALSE))
      {
        rprob <- rep(1/nl,nl)
        length <- length + nl
      }
      else
      {
        rprob <- rootprobability/sum(rootprobability)
        nlength <- length + nl
        
      }        
      talpha <- r
      if(is.numeric(alpha))
        talpha <- alpha
      length <- length + 1
    }
    else
    {
      length <- 0
      if(identical(rateparameters,TRUE))
      {
        p <- r[(length+1):(length+np)]
        length <- length + np
      }
      if(!identical(rateparameters,TRUE))
      {
        p <- rateparameters
        length <- length + np
      }
      
      if(identical(rootprobability, TRUE))
      {
        rprob <- r[(length+1):(length+nl)]
        if(sum(rprob)<=0)
          return(1e50)
        else
          rprob <- rprob/sum(rprob)
        length <- length + nl
      } 
      else if(identical(rootprobability, FALSE))
      {
        rprob <- rep(1/nl,nl)
        length <- length + nl
      }
      else
      {
        rprob <- rootprobability/sum(rootprobability)
        length <- length + nl
      }      
      if(!identical(alpha,FALSE))
      {
        if(identical(alpha, TRUE))
          talpha <- r[length+1]
        else if(is.numeric(alpha))
          talpha <- alpha
        
        length <- length + 1
      }  
    }   
    if(!identical(rateparameters,TRUE)&&length(rateparameters)==np)
      p <- rateparameters
    if(!identical(rateparameters,TRUE)&&!(length(rateparameters)==np) )
      stop("The number of rate parameters does not fit into the rate matrix, which is," ,np)
    if(identical(rateparameters,FALSE))
      stop("rateparameters cannot be set to FALSE") 
    Q[] <- c(p, 0)[rate]  
    if (is.matrix(model))
    {     
      if(is.matrix(const))
        for(i in 1:length(const[,1]))
          Q[const[i,1] == model] <- as.character(const[i,2])      
      if(is.list(const))
        for(i in 1:length(const))
          Q[ names(const)[i] ==model] <- as.character(const[[i]])      
      Q <- matrix(suppressWarnings(as.numeric(Q)), nl,nl)
    }     
    diag(Q) <- rep(0, nl)   
    if(  reversible || model== "GTR" )
      Q <- t(t(Q)*rprob)
    diag(Q) <- -rowSums(Q)
    if(any(is.na(rprob)))
      return(1e50)
    tQ <- Q    
    div <- sum(diag(Q)*rprob)    
    if(div==0)
      return(1e50)    
    if(!identical(alpha,FALSE))
    {
      Q <- Q/(-div)
      #gamma rate 
      cat=gammak
      rateht<-matrix(NA, 1,cat)
      for(i in 1:cat)
      {
        num1=(i-0.5)/cat
        rateht[i]=qgamma(num1, shape=talpha, rate=talpha)
      } 
      rateht=rateht/sum(rateht)*cat      
      if(any(is.na(rateht)))
        return(1e+50)      
      #gamma rate end
    }   
    decompo <- eigen(Q)
    lambda <- decompo$values
    GAMMA <- decompo$vectors
    tempInvGAMMA <- try(solve(GAMMA))
    if( class(tempInvGAMMA)=="try-error"   )
      return(1e50)
    else
      invGAMMA<- tempInvGAMMA    
    for(jj in 1:NROW)
    {
      liks=likslist[[jj]]      
      if(identical(alpha,FALSE))
      {
        comp <- numeric(nb.tip + nb.node+1)
        for (i in seq(from = 1, by = 2, length.out = nb.node)) 
        {
          j <- i + 1L
          anc <- e1[i]
          des1 <- e2[i]
          des2 <- e2[j]
          v.l <- GAMMA %*% diag(exp(lambda * EL[i])) %*% 
            invGAMMA %*% liks[des1, ]
          v.r <- GAMMA %*% diag(exp(lambda * EL[j])) %*% 
            invGAMMA %*% liks[des2, ]
          v <- v.l * v.r
          comp[anc] <- sum(v)       
          liks[anc, ] <- v/comp[anc]          
        }
        likslist[[jj]] <- liks
        comp[nb.tip+nb.node +1] <- sum(liks[nb.tip+1, ]  * rprob)        
        if(!zerocorrection)
        {
          likeli =  sum(log(comp[-TIPS]))
        }       
        if(zerocorrection==TRUE)
        {
          zliks = cbind(matrix(1, nrow(liks),1),matrix(0, nrow(liks), ncol(liks)-1))
          for (i in seq(from = 1, by = 2, length.out = nb.node)) 
          {
            j <- i + 1L
            anc <- e1[i]
            des1 <- e2[i]
            des2 <- e2[j]            
            v.l <- GAMMA %*% diag(exp(lambda * EL[i])) %*% 
              invGAMMA %*% zliks[des1, ]
            v.r <- GAMMA %*% diag(exp(lambda * EL[j])) %*% 
              invGAMMA %*% zliks[des2, ]
            v <- v.l * v.r
            zliks[anc, ] <- v            
          }          
          zlikslist[[jj]] <- zliks
          zlikeli<- sum(zliks[nb.tip+1, ]*rprob)          
          if(zlikeli==1 || is.na(zlikeli))
          {
            return(1e+50)
          }
          likeli <- sum(log(comp[-TIPS]))-log(1-zlikeli)
        }
        likelivector[jj]= likeli
      }
      if(!identical(alpha,FALSE))
      {
        mu0 <- -div
        likeli <-0
        for(w in 1:cat)
        {
          mu <- mu0*rateht[w]         
          for (i in seq(from = 1, by = 2, length.out = nb.node))           {
            j <- i + 1L
            anc <- e1[i]
            des1 <- e2[i]
            des2 <- e2[j]
            v.l <- GAMMA %*% diag(exp(lambda *mu* EL[i])) %*% 
              invGAMMA %*% liks[des1, ]
            v.r <- GAMMA %*% diag(exp(lambda *mu* EL[j])) %*% 
              invGAMMA %*% liks[des2, ]
            v <- v.l * v.r
            liks[anc, ] <- v
          }          
          likeli <- likeli + sum(liks[nb.tip+1, ]*rprob)/gammak          
        }        
        likslist[[jj]] = liks        
        if(!zerocorrection)
          likeli <- log(likeli)        
        if(zerocorrection)
        {
          zlikeli<- 0
          zliks = cbind(matrix(1, nrow(liks),1),matrix(0, nrow(liks), ncol(liks)-1))          
          for(w in 1: cat)
          {  
            mu <- mu0*rateht[w]            
            for (i in seq(from = 1, by = 2, length.out = nb.node)) 
            {
              j <- i + 1L
              anc <- e1[i]
              des1 <- e2[i]
              des2 <- e2[j]
              
              v.l <- GAMMA %*% diag(exp(lambda *mu* EL[i])) %*% 
                invGAMMA %*% zliks[des1, ]
              v.r <- GAMMA %*% diag(exp(lambda *mu* EL[j])) %*% 
                invGAMMA %*% zliks[des2, ]
              v <- v.l * v.r
              zliks[anc, ] <- v              
            }            
            zlikeli<- zlikeli+ sum(zliks[nb.tip+1, ] * rprob)/gammak
          }
          zlikslist[[jj]] <- zliks          
          likeli <- ( log(likeli) -log(1-zlikeli)  )          
        }        
      }      
      likelivector[jj] <-  likeli      
    }       
    output = -2*sum(likelivector)    
    if(output.liks)
      return(likslist)
    if(output.zliks)
      return(zlikslist) 
    if(output.div)
      return(div)     
    if(output.Q)
      return(tQ)      
    if (is.na(output)||is.infinite(output)) 
      10e50
    else Re(output)
  }   
  length <- 0
  rateParametersIndicator <-0
  rootProbabilityIndicator <-0
  alphaIndicator <-0
  if(identical(rateparameters,TRUE))
  {
    rateParametersStart <- length+1
    rateParametersEnd <- length+np
    rateParametersIndicator <- 1
    length <- length + np
  }
  if(!identical(rateparameters,TRUE))
  {
    rateParametersStart <- length+1
    rateParametersEnd <- length+np
    rateParametersIndicator <- 0
    length <- length + np
  }
  if(identical(rootprobability, TRUE))
  {
    rootProbabilityStart <- length+1
    rootProbabilityEnd <- length+nl
    rootProbabilityIndicator <- 1
    length <- length + nl
  }
  if(!identical(rootprobability, TRUE))
  {
    rootProbabilityStart <- length+1
    rootProbabilityEnd <- length+nl
    rootProbabilityIndicator <-0
    length <- length + nl
  }
  if(identical(irootprobability , FALSE) || identical(irootprobability, TRUE))
    irprob <- rep(1/nl, nl)
  else if(is.numeric(irootprobability))
    irprob <- irootprobability
  else if(identical(toupper(irootprobability), "RANDOM"))
    irprob <- rep(1/nl, nl)
  if(identical(alpha,FALSE))
  {
    out<-list()
    if(!identical(rootprobability, TRUE))
    {
      initialPoint <- c(rep(ip,np))
      out <- nlminb(initialPoint, function(r)dev(r,turnAlphaToFalse=TRUE), lower= rep(0,np), upper = rep(1e50,np) )
    }
    if(identical(rootprobability, TRUE))
    {
      if(identical(irootprobability, TRUE) )
      {
        if(delta<nl)
        {
          multichoose <-  function(n, k) 
          {      
            if(k < 0 || n <= 0)
              return(NULL);      
            if(k == 0) {
              return(list(numeric(n)))
            }      
            if(n == 1) {
              return(list(k));
            }      
            arr <- multichoose(n - 1, k);
            arr2 <- multichoose(n, k - 1);      
            out <- list();      
            for(i in 1:(length(arr)) ) {
              arr[[i]]<- c(0, arr[[i]])
              out <-  c(out, arr[i]);
            }
            for(i in 1:(length(arr2))  ) {
              arr2[[i]][1] <- arr2[[i]][1]+1
              out <- c(out,( arr2[i]))
            }
            return(out);
          }
          minimumScore<-1e50
          minimumScoreRprob<-numeric(nl)
          minimumScoreRateparameters <-numeric()
          multichooseTemp<- multichoose(nl, delta)
          initialRootprobabilityPoint<- rep(1/nl,nl)
          initialPoint<-c(rep(ip,np), initialRootprobabilityPoint )
          temp<- nlminb(initialPoint, function(r)dev(r,turnAlphaToFalse=TRUE), lower= rep(0,nl+np), upper = rep(1e50,nl+np)  )
          minimumScore <- temp$objective
          minimumScoreRprob <-temp$par[rootProbabilityStart:rootProbabilityEnd]/sum(temp$par[rootProbabilityStart:rootProbabilityEnd] )
          minimumScoreRateParameters<- temp$par[1:np]   
          out$par <- c(minimumScoreRateParameters, minimumScoreRprob)
          out$objective <- minimumScore
          if(showprogress)
          {
            show("So far, the highest log likelhiood achieved:")
            show((-1/2)*minimumScore)
            cat("[1] Prior root probability estimated when its initial point is set to","(", initialRootprobabilityPoint,"):", "\n" )
            if(identical(rootprobability,TRUE))
              show(minimumScoreRprob )
            else if (identical(rootprobability,FALSE))
              show(rep(1/nl,nl))
            else
              show(rootprobability/sum(rootprobability))                       
            cat("The above data will be updated if DiscML has found a prior root probability that would result in
                higher log likelihood. Click here to see if there is any update.\n") 
          }
          }
        
        for(i in 1:length(multichooseTemp))
        {       
          initialRootprobabilityPoint <- multichooseTemp[[i]]/delta
          initialPoint <- c(rep(ip,np), initialRootprobabilityPoint)
          temp<-nlminb(initialPoint, function(r)dev(r, turnAlphaToFalse=TRUE), lower= rep(0,nl+np), upper = rep(1e50,nl+np) )
          if(temp$objective < minimumScore)
          {
            minimumScore <- temp$objective
            minimumScoreRprob <-temp$par[rootProbabilityStart:rootProbabilityEnd]/sum(temp$par[rootProbabilityStart:rootProbabilityEnd] )
            minimumScoreRateParameters<- temp$par[1:np]
            out$par <- c(minimumScoreRateParameters, minimumScoreRprob)
            out$objective <- minimumScore    
            if(showprogress)
            {
              show("#####################PROGRESS####################")
              show("So far, the highest log likelhiood achieved:")
              show((-1/2)*minimumScore)
              cat("[1] Prior root probability estimated when its initial point is set to","(", initialRootprobabilityPoint,"):", "\n" )
              
              if(identical(rootprobability,TRUE))
                show(minimumScoreRprob )
              else if (identical(rootprobability,FALSE))
                show(rep(1/nl,nl))
              else
                show(rootprobability/sum(rootprobability))    
              cat("The above data will be updated if DiscML has found a prior root probability that would result in
                  higher log likelihood. Click here to see if there is any update.\n")   
            }
            }
          }
        }
      if(!identical(irootprobability, TRUE))
      {
        if(identical(irootprobability, FALSE) || is.numeric(irootprobability))
        {
          initialPoint<-c(rep(ip,np), irprob )
          temp<- nlminb(initialPoint, function(r)dev(r,turnAlphaToFalse=TRUE), lower= rep(0,nl+np), upper = rep(1e50,nl+np)  )
          minimumScore <- temp$objective
          minimumScoreRprob <-temp$par[rootProbabilityStart:rootProbabilityEnd]/sum(temp$par[rootProbabilityStart:rootProbabilityEnd] )
          minimumScoreRateParameters<- temp$par[1:np]       
          out$par <- c(minimumScoreRateParameters, minimumScoreRprob)
          out$objective <- minimumScore 
        }
        else if (identical(toupper(irootprobability), "RANDOM"))
        {
          minimumScore <- 1e50
          randomVectors <- matrix(0, ivnum, nl)
          for(i in 1:ivnum)
          { 
            
            randomVectors[i,] <-  sample(1:100,nl, replace=TRUE)
            while(sum(randomVectors[i,])==0)
              randomVectors[i,] <-  sample(1:100,nl, replace=TRUE)
            randomVectors[i,] <-  randomVectors[i,]/sum( randomVectors[i,])
            initialRootprobabilityPoint <- randomVectors[i,]
            initialPoint <- c(rep(ip,np), initialRootprobabilityPoint)
            temp<-nlminb(initialPoint, function(r)dev(r, turnAlphaToFalse=TRUE), lower= rep(0,nl+np), upper = rep(1e50,nl+np) )
            if(temp$objective < minimumScore)
            {
              minimumScore <- temp$objective
              minimumScoreRprob <-temp$par[rootProbabilityStart:rootProbabilityEnd]/sum(temp$par[rootProbabilityStart:rootProbabilityEnd] )
              minimumScoreRateParameters<- temp$par[1:np]
              out$par <- c(minimumScoreRateParameters, minimumScoreRprob)
              out$objective <- minimumScore    
              if(identical(showprogress,TRUE))
              {
                show("#####################PROGRESS####################")
                show("So far, the highest log likelhiood achieved:")
                show((-1/2)*minimumScore)
                cat("[1] Prior root probability estimated when its initial point is set to","(", initialRootprobabilityPoint,"):", "\n" )
                
                if(identical(rootprobability,TRUE))
                  show(minimumScoreRprob )
                else if (identical(rootprobability,FALSE))
                  show(rep(1/nl,nl))
                else
                  show(rootprobability/sum(rootprobability))    
                cat("The above data will be updated if DiscML has found a prior root probability that would result in
                    higher log likelihood. Click here to see if there is any update.\n")   
              }
              }
            }
          }
        else
          stop("Please make sure you typed correct option for irootprobabilty.")
        }
        }
  }
  if(!identical(alpha,FALSE))
  {
    gammaOut <- list()
    out <- list()
    alphaStart<-length+1
    if(!identical(rootprobability, TRUE))
    {  
      initialPoint<-c(rep(ip,np), irprob, ialpha) 
      preFinal<- nlminb(initialPoint, function(r)dev(r), lower= c(rep(0,np+nl), 0.0001), upper = c(rep(1e50,np+nl), 100000)  )
      minimumScore <- preFinal$objective
      minimumScoreRprob <-preFinal$par[rootProbabilityStart:rootProbabilityEnd]/sum(preFinal$par[rootProbabilityStart:rootProbabilityEnd] )
      if(identical(rootprobability, FALSE))
        minimumScoreRprob<-rep(1/nl,nl)
      minimumScoreRateParameters<- preFinal$par[1:np]   
      out$par[1:length] <- preFinal$par[1:length]
      gammaOut$par <- preFinal$par[length+1]
      gammaOut$objective <- preFinal$objective
    }
    if(identical(rootprobability, TRUE))
    {
      multichoose <-  function(n, k) 
      {      
        if(k < 0 || n <= 0)
          return(NULL);      
        if(k == 0) {
          return(list(numeric(n)))
        }      
        if(n == 1) {
          return(list(k));
        }      
        arr <- multichoose(n - 1, k);
        arr2 <- multichoose(n, k - 1);      
        out <- list();      
        for(i in 1:(length(arr)) ) {
          arr[[i]]<- c(0, arr[[i]])
          out <-  c(out, arr[i]);
        }
        for(i in 1:(length(arr2))  ) {
          arr2[[i]][1] <- arr2[[i]][1]+1
          out <- c(out,( arr2[i]))
        }
        return(out);
      }
      minimumScore<-1e50
      minimumScoreRprob<-numeric(nl)
      minimumScoreRateparameters <-numeric()
      multichooseTemp<- multichoose(nl, delta)
      if(identical(irootprobability, TRUE) )
      {  
        gammaOut <- list()
        out <- list()
        if(delta<nl)
        {
          initialRootprobabilityPoint<- rep(1/nl,nl)
          initialPoint<-c(rep(ip,np), initialRootprobabilityPoint , ialpha)
          temp<- nlminb(initialPoint, function(r)dev(r), lower= c(rep(0,nl+np), 0.0001), upper = c(rep(1e50,nl+np), 100000)  )
          minimumScore <- temp$objective
          minimumScoreRprob <-temp$par[rootProbabilityStart:rootProbabilityEnd]/sum(temp$par[rootProbabilityStart:rootProbabilityEnd] )
          minimumScoreRateParameters<- temp$par[1:np] 
          out$par[1:length] <- temp$par[1:length]
          gammaOut$par <- temp$par[length+1]
          gammaOut$objective <- temp$objective
          if(identical(showprogress, TRUE))
          {
            show("So far, the highest log likelhiood achieved:")
            show((-1/2)*minimumScore)
            cat("[1] Prior root probability estimated when its initial point is set to","(", initialRootprobabilityPoint,"):", "\n" )
            if(identical(rootprobability,TRUE))
              show(minimumScoreRprob )
            else if (identical(rootprobability,FALSE))
              show(rep(1/nl,nl))
            else
              show(rootprobability/sum(rootprobability))                
            cat("The above data will be updated if DiscML has found a prior root probability that would result in
                higher log likelihood. Click here to see if there is any update.\n")             
          }
          }
        for(i in 1:length(multichooseTemp))
        {       
          initialRootprobabilityPoint <- multichooseTemp[[i]]/delta
          initialPoint <- c(rep(ip,np), initialRootprobabilityPoint, ialpha)
          temp<-nlminb(initialPoint, function(r)dev(r), lower= c(rep(0,nl+np), 0.0001), upper = c(rep(1e50,nl+np), 100000 ))
          if(temp$objective < minimumScore)
          {
            minimumScore <- temp$objective
            minimumScoreRprob <-temp$par[rootProbabilityStart:rootProbabilityEnd]/sum(temp$par[rootProbabilityStart:rootProbabilityEnd] )
            minimumScoreRateParameters<- temp$par[1:np]
            out$par[1:length] <- temp$par[1:length]
            gammaOut$par <- temp$par[length+1]
            gammaOut$objective <- temp$objective
            if(identical(showprogress, TRUE))
            {
              show("#####################PROGRESS####################")
              show("So far, the highest log likelhiood achieved:")
              show((-1/2)*minimumScore)
              cat("[1] Prior root probability estimated when its initial point is set to","(", initialRootprobabilityPoint,"):", "\n" )
              
              if(identical(rootprobability,TRUE))
                show(minimumScoreRprob )
              else if (identical(rootprobability,FALSE))
                show(rep(1/nl,nl))
              else
                show(rootprobability/sum(rootprobability))
              cat("The above data will be updated if DiscML has found a prior root probability that would result in
                  higher log likelihood. Click here to see if there is any update.\n")      
            }
            }
          }
        }
      if(identical(irootprobability, FALSE)|| is.numeric(irootprobability))
      {
        initialPoint<-c(rep(ip,np), irprob, ialpha )
        temp<- nlminb(initialPoint, function(r)dev(r,turnAlphaToFalse=TRUE), lower= c(rep(0,nl+np+1), 0.0001 ), upper = rep(1e50,nl+np) , 100000)
        minimumScore <- temp$objective
        minimumScoreRprob <-temp$par[rootProbabilityStart:rootProbabilityEnd]/sum(temp$par[rootProbabilityStart:rootProbabilityEnd] )
        minimumScoreRateParameters<- temp$par[1:np]  
        
        gammaOut <- list()
        out <- list()
        out$par[1:length] <- temp$par[1:length]
        gammaOut$par <- temp$par[length+1]
        gammaOut$objective <- temp$objective
      }
      else if (identical(toupper(irootprobability), "RANDOM"))
      {
        randomVectors <- matrix(0, ivnum, nl)
        for(i in 1:ivnum)
        {       
          randomVectors[i,] <-  sample(1:100,nl, replace=TRUE)
          while(sum(randomVectors[i,])==0)
            randomVectors[i,] <-  sample(1:100,nl, replace=TRUE)
          randomVectors[i,] <-  randomVectors[i,]/sum( randomVectors[i,])
          initialRootprobabilityPoint <- randomVectors[i,]
          initialPoint <- c(rep(ip,np), initialRootprobabilityPoint,ialpha)
          temp<-nlminb(initialPoint, function(r)dev(r), lower= c(rep(0,nl+np), 0.0001), upper = c(rep(1e50,nl+np), 100000 ))
          if(temp$objective < minimumScore)
          {
            minimumScore <- temp$objective
            minimumScoreRprob <-temp$par[rootProbabilityStart:rootProbabilityEnd]/sum(temp$par[rootProbabilityStart:rootProbabilityEnd] )
            minimumScoreRateParameters<- temp$par[1:np]
            out$par[1:length] <- temp$par[1:length]
            gammaOut$par <- temp$par[length+1]
            gammaOut$objective <- temp$objective
            if(identical(showprogress, TRUE))
            {
              show("#####################PROGRESS####################")
              show("So far, the highest log likelhiood achieved:")
              show((-1/2)*minimumScore)
              cat("[1] Prior root probability estimated when its initial point is set to","(", initialRootprobabilityPoint,"):", "\n" )
              
              if(identical(rootprobability,TRUE))
                show(minimumScoreRprob )
              else if (identical(rootprobability,FALSE))
                show(rep(1/nl,nl))
              else
                show(rootprobability/sum(rootprobability))
              cat("The above data will be updated if DiscML has found a prior root probability that would result in
                  higher log likelihood. Click here to see if there is any update.\n")      
            }
            }
          }
        }
      else
        stop("Please make sure you typed correct option for irootprobabilty.")
      }
      }
  if(!identical(alpha,FALSE))
    totalOutObjective <- gammaOut$objective
  else
    totalOutObjective <- out$objective
  if(is.numeric(rootprobability))
  {
    rprob<- rootprobability/sum(rootprobability)
    out$par[rootProbabilityStart:rootProbabilityEnd]<-rprob
  }
  if(identical(rootprobability, TRUE))
  {
    rprob <- out$par[rootProbabilityStart:rootProbabilityEnd]
    rprob <- rprob/sum(rprob)
  }
  else if(identical(rootprobability, FALSE))
  {
    rprob <- rep(1/nl, nl)
    out$par[rootProbabilityStart:rootProbabilityEnd]<-rprob
  }
  else
  {
    rprob <- rootprobability
    rprob <- rprob/sum(rprob)
    out$par[rootProbabilityStart:rootProbabilityEnd]<-rprob
  }
  if(identical( rateparameters,TRUE))
  {
    gamma<-0
    if(!identical(alpha,FALSE))
      gamma<-1
    if(identical(alpha,FALSE))
      gammaOut<-NULL
    totalOut<- c(out$par,  rep(gammaOut$par, gamma))
    lmu <- -dev(totalOut , output.div=TRUE)
  }
  else
  {
    gamma<-0
    if(!identical(alpha,FALSE))
      gamma<-1
    if(identical(alpha,FALSE))
      gammaOut<-NULL
    out$par[1:np]<-minimumScoreRateParameters<- rateparameters
    totalOut<- c(rateparameters, rprob, rep(gammaOut$par, gamma))
    lmu <- -dev(totalOut , output.div=TRUE)
  }
  obj$loglik <- -totalOutObjective/2
  obj$rates <-out$par[1:np]
  if(reversible && !is.matrix(model) )
  {
    if( toupper(substring(model, 1,2)) =="BD")
    {
      QQ <- dev(totalOut, output.Q =TRUE)
      np2 <- 2*(nl-1)
      rate2 <- rate
      rate2[row(rate2)!=col(rate2)] <- np2+1
      for(i in 1:(nl-1))
      {
        rate2[i, i+1]=i
        rate2[i+1, i]= nl-1+i
      }
      rate2[is.na(rate2)] <- np2+1
      diag(rate2) <- rep(np2+1, nl)
      index.matrix <- rate2
      obj$rates <- c( QQ[ col(QQ) == row(QQ)+1 ] , QQ[ col(QQ) == row(QQ)-1 ]    )
      rootProbabilityStart<-np2+1
      rootProbabilityEnd <- np2 +nl
    }
    else
    {
      QQ <- dev(totalOut, output.Q =TRUE)
      np2 <- nl * (nl - 1)
      par = numeric()
      par= QQ[row(QQ) != col(QQ)]
      rate2= matrix(NA,nl,nl)
      rate2[col(rate2) != row(rate2)] <- 1:np2
      index.matrix <- rate2
      obj$rates <- as.numeric(QQ[row(QQ)!=col(QQ)])
      rootProbabilityStart<-np2+1
      rootProbabilityEnd <- np2 +nl
    }
  }
  if(  is.matrix(model) )
  {
    QQ <- dev(totalOut, output.Q =TRUE)
    np2 <- nl * (nl - 1)
    par = numeric()
    par= QQ[row(QQ) != col(QQ)]
    rate2= matrix(NA,nl,nl)
    rate2[col(rate2) != row(rate2)] <- 1:np2
    index.matrix <- rate2
    obj$rates <- as.numeric(QQ[row(QQ)!=col(QQ)])
    rootProbabilityStart <- np2+1
    rootProbabilityEnd <- np2+nl
  }
  if(identical(mu, TRUE))
    obj$rates <-obj$rates/lmu
  else
    obj$rates <-obj$rates
  oldwarn <- options("warn")
  options(warn = -1)
  if(identical(rootprobability, TRUE))
  {
    out.nlm <- try(nlm(function(r) dev(r, turnAlphaToFalse=TRUE), p = out$par, 
                       iterlim = 1, stepmax = 0, hessian = TRUE))
    options(oldwarn)
    obj$se <- if (class(out.nlm) == "try-error") {
      warning("model fit suspicious: gradients apparently non-finite")
      rep(NaN, length)
    }
    else if( class(try(.getSEs(out.nlm))) =="try-error")
      rep(NaN, length)
    else
      suppressWarnings(.getSEs(out.nlm))
  }
  else
  {
    out.nlm <- try(nlm(function(r) dev(r, turnAlphaToFalse=TRUE, turnRootProbabilityToFalse = TRUE), p = out$par[1:np], 
                       iterlim = 1, stepmax = 0, hessian = TRUE))
    options(oldwarn)
    obj$se <- if (class(out.nlm) == "try-error") {
      warning("model fit suspicious: gradients apparently non-finite")
      c(rep(NaN, np))
    }
    else if( class(try(.getSEs(out.nlm))) =="try-error")
      c(rep(NaN, np))
    else
      .getSEs(out.nlm)
    
    obj$se <- c(obj$se, rep(0,nl))
  }
  if(!identical(alpha,FALSE)){
    gammaOut.nlm <- try(nlm(function(r) dev(r, manualExceptAlpha = out$par), p = gammaOut$par, 
                            iterlim = 1, stepmax = 0, hessian = TRUE))
    options(oldwarn)
    gammaObjSe <- if (class(gammaOut.nlm) == "try-error") {
      warning("model fit suspicious: gradients apparently non-finite")
      NaN
    }
    else if( try(suppressWarnings(.getSEs(gammaOut.nlm))) =="try-error")
      NaN
    else
      suppressWarnings(.getSEs(gammaOut.nlm))
    obj$se <- c(obj$se , gammaObjSe)
    
  }
  if(reversible)
  {
    index.matrix[ is.na(index.matrix)] <- np2+1
    index.matrix[index.matrix == np2+1] <- rep(NA, length(index.matrix[index.matrix == np2+1]))
  }
  else
  {
    index.matrix[ is.na(index.matrix)] <- np+1
    index.matrix[index.matrix == np+1] <- rep(NA, length(index.matrix[index.matrix == np+1]))
    if(is.matrix(model))
    {
      index.matrix[ is.na(index.matrix)] <- np+1
      index.matrix[row(index.matrix) == col(index.matrix)] <- rep(NA, nl)
    }
  }
  obj$index.matrix <- index.matrix
  tempQQ <- matrix(NA, nl ,nl)
  rate3<- rate
  rate3[rate3==0] <- rep(np+1 , length(rate3[rate3==0]) )
  tempQQ[] <-  c(obj$se[1:np], 0)[rate3]
  if( !(toupper(substring(model, 1,2)) == "BD") && (if(!is.matrix(model)){ toupper(model) == "ARD" }  || reversible ) )
    obj$se <- c(tempQQ[col(tempQQ) != row(tempQQ)], obj$se[(np+1):length(obj$se)])
  if( toupper(substring(model, 1, 2) ) =="BD" && reversible || toupper(model) =="BDARD")
  {
    obj$se <- c(tempQQ[col(tempQQ) ==1 + row(tempQQ)], tempQQ[col(tempQQ) == -1 + row(tempQQ)]   , obj$se[(np+1):length(obj$se)])
    rootProbabilityStart <- length(  c(tempQQ[col(tempQQ) ==1 + row(tempQQ)], tempQQ[col(tempQQ) == -1 + row(tempQQ)] ) )+1
    rootProbabilityEnd <- rootProbabilityStart +nl-1
  }
  if (CI)
  {
    obj$tobj = matrix(0,NROW, nl)
    obj$lik.anc <- dev(totalOut,  output.liks=TRUE)
    for(i in 1:NROW)
    {
      colnames(obj$lik.anc[[i]]) <- lvls
      obj$tobj[i,] <- obj$lik.anc[[i]][nb.tip+1,]
      obj$lik.anc[[i]]=obj$lik.anc[[i]]
      names(obj$lik.anc)[i] = paste("Data ", i, " :")
    }
    cc=numeric()
    colnames(obj$tobj) = lvls
    for( i in 1:NROW)
      cc[i]=paste("Data", i, ": ") 
    rownames(obj$tobj) <- cc
    obj$tobj=obj$tobj/rowSums(obj$tobj)
  }
  if(identical(rootprobability, TRUE))
  {
    srprob <- obj$se[rootProbabilityStart:rootProbabilityEnd]
  }
  else if(identical(rootprobability, FALSE))
  {
    srprob <- rep(0, nl)
  }
  else if(is.numeric(rootprobability))
  {
    srprob <- rep(0,nl)
  }
  obj$call <- match.call()
  class(obj) <- "DiscML"
  prpdata<- data.frame(characters = as.integer(lvls), estimates= rprob, "std-err" =srprob)
  obj$prp = prpdata
  
  
  if(!identical(alpha,FALSE))
  {
    if(identical(alpha,TRUE))
    {
      lalpha <- gammaOut$par
      salpha <- obj$se[alphaStart]
    }
    else if(identical(alpha,FALSE))
    {
      lalpha <- 1
      salpha <- 0
    }
    else
    {
      lalpha <- alpha
      salpha <- 0
    }
  }
  if(plotalpha && identical(alpha,TRUE))
  {
    if(lalpha > 10000)
      xx= seq(0.001, 10000, by = 10)
    if(lalpha <= 10000)
      xx = seq(0.001, lalpha*1.4, by= lalpha/100)
    newf<- function(xx)
    {
      tt= numeric()
      for(i in 1:length(xx))
      {
        tt[i] = -1/2*dev( c(totalOut[-length(totalOut)] ,xx[i])   )
      }
      return(tt)
    }
    plot(xx, newf(xx), xlab= "alpha", ylab= "log likelihood, given gamma(shape = alpha,scale=alpha)", type="l") 
  }
  obj$lvls <- lvls
  obj$alpha <- alpha
  if(!identical(alpha, FALSE))
  {
    obj$lalpha <- lalpha
    obj$salpha <- salpha
  }
  obj$mu <- lmu
  obj$ismu <- ismu
  obj$prp <- prpdata
  return(obj)
  }
