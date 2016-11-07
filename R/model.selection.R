### H function
calculate.H.function= function(chain, result, lower.upper=lower.upper, verbose = FALSE)
{
	if(class(result) == "epi") 
	{
		n = nrow(result$data)
	}else{
		cat("epi.object should belong to the epi class. \n ")
	}
	
	if(verbose)
	{
   		m <- paste(c("Model Selection using H.function for chain ", chain, " ... in progress:", floor(100 * chain/n), "%"), collapse="")
   		cat(m, "\r")
       	flush.console()
    }
	
	pr	<- 0
	for(i in 1:n)
	{
		Sigma	= result$ES[[chain]]
		Z		= result$Z[[chain]]
		pr <- mean(dmvnorm(Z, sigma= Sigma, log=TRUE)) - log(ptmvnorm(lower=lower.upper$lower[i, ], upper= lower.upper$upper[i, ], sigma= Sigma))
		pr <- pr + pr
	}
	return(pr)
}

model.selection = function( result, lower.upper, ebic.gamma = 0.5,  ncores = 1, EM.select = FALSE, verbose = TRUE)
{

	if(class(result) == "epi") 
	{
		p	 = ncol(result$path[[1]])
		nrho = length(result$rho)
		n	 = nrow(result$data)

	}else{
		cat("epi.object should belong to the epi class. \n ")
	} 

	if(EM.select) 
	{
		H = c()
		if(ncores == 1)
		{
			all.H <- lapply(1:nrho, function(i){calculate.H.function(chain=i, result=result, lower.upper=lower.upper, verbose = FALSE)})
		}else{
			cl <- makeCluster(ncores)
			all.H <- parLapply(cl, 1:nrho, get("calculate.H.function"), result=result, lower.upper=lower.upper, verbose = FALSE) 
			stopCluster(cl)
		}
		H <- unlist(all.H)
		loglik =  result$loglik - H
		rm(all.H, H)

	}else{

		loglik = result$loglik 
	}

	result$df <- sapply(1:nrho, function(x) sum(as.matrix(result$Theta[[x]])[upper.tri(result$Theta[[x]])] != 0 ))

	if(verbose)
	{
		cat("Calculating extended Bayesian information criterion (ebic) selection at the EM convergence ...")
		flush.console()
	}
	 
	result$ebic.scores	= - 2 * loglik + ( log(n) * result$df ) + ( 4 * ebic.gamma * log(p) * result$df )
	result$opt.index	= which.min(result$ebic.scores)
	result$opt.Theta	= result$Theta[[result$opt.index]]
	result$opt.path		= abs(sign(result$opt.Theta)) - diag(rep(1,p))
	result$opt.Sigma	= result$Sigma[[result$opt.index]]
	result$opt.rho		= result$rho[result$opt.index]

	if(verbose)
	{
		cat("done.\n")
		flush.console()
	}
	return(result)
}

