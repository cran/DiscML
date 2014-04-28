print.DiscML <- function (x, digits = 4, ...) 
{
  cat("\n  DiscML (Maximum Likelihood Method for Discrete Character States)\n\n")
  cat("Call: ")
  print(x$call)
  cat("\n")
  if (!is.null(x$loglik)) 
    cat("    Log-likelihood:", x$loglik, "\n\n")
  if (!is.null(x$resloglik)) 
    cat("    Residual log-likelihood:", x$resloglik, "\n\n")
  ratemat <- x$index.matrix
  if (is.null(ratemat)) {
    class(x) <- NULL
    x$resloglik <- x$loglik <- x$call <- NULL
    print(x)
  }
  else {
    dimnames(ratemat)[1:2] <- dimnames(x$lik.anc)[2]
    cat("Rate index matrix:\n")
    colnames(ratemat)<-x$lvls
    rownames(ratemat)<-x$lvls
    print(ratemat, na.print = ".", quote=FALSE)
    cat("\n")
    npar <- length(x$rates)
    estim <- data.frame(1:npar, round(x$rates, digits), round(x$se[1:npar],   digits))
    cat("Parameter estimates:(type '...$rates' to get this table): \n")
    names(estim) <- c("rate index", "estimate", "std-err")
    print(estim, row.names = FALSE)
    if(identical(x$ismu, TRUE))
      cat("\nThe value of mu estimated (type '...$mu' to get this.): ", x$mu, "\n\n")
    if (!is.null(x$lik.anc)) {
      writeLines("\nScaled likelihoods at the root for each site (type '...$lik.anc' to get them for unscaled version of all nodes):\n")
      print(round(Re(x$tobj), 7) )
      x$tobj=NULL
    }
    cat("\nPrior root probability estimated/used (type '...$prp' to get this table.): \n")
    show(x$prp)
    cat("\n\n")
  }
  if(!identical(x$alpha,FALSE)){
    if(identical(x$alpha,TRUE)) 
    {
      if(x$lalpha < 9999)
        cat("The value of alpha parameter of gamma used (or estimated) (type '...$alpha' to get this.):\n",
            x$lalpha, "with the standard error of", x$salpha,"\n\n")
      else
        cat("The value of alpha parameter of gamma estimated is over 10000,
          so the optimization of the alpha paramter may not be suitable \n\n")   
    }
    else
      cat("The value of alpha parameter of gamma used (or estimated) (type '...$alpha' to get this.):\n",
          x$lalpha, "with the standard error of", x$salpha,"\n\n")
  }
}