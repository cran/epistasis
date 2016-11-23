calculate.loglik_Y = function(chain, result, lower.upper=lower.upper, verbose = TRUE)
{
	if(class(result) == "epi") 
	{
		n = nrow(result$data)
	}else{
		cat("epi.object should belong to the epi class. \n ")
	}
	
	if(verbose)
	{
		m <- paste(c("Calculating observed log-likelihood in progress ... ", floor(100 * chain/n), "%"), collapse="")
   		cat(m, "\r")
       	flush.console()
    }
	
	loglik_Y <- 0
	for(i in 1:n)
	{
		Sigma	 <- result$ES[[chain]]
		loglik <- log(ptmvnorm(lowerx=lower.upper$lower[i, ], upperx= lower.upper$upper[i, ], sigma= Sigma))
		loglik_Y   <- loglik_Y + loglik
	}
	
	return(loglik_Y)		
}

model.selection = function( result, criterion, lower.upper, ebic.gamma = 0.5,  ncores = 1, loglik_Y=TRUE, verbose = TRUE)
{

	if(class(result) == "epi") 
	{
		p	 = ncol(result$path[[1]])
		nrho = length(result$rho)
		n	 = nrow(result$data)

	}else{
		cat("epi.object should belong to the epi class. \n ")
	} 

	if(loglik_Y)
	{
		loglik = c()
		if(ncores == 1)
		{
			all.loglik <- lapply(1:nrho, function(i){calculate.loglik_Y(chain=i, result=result, lower.upper=lower.upper, verbose = TRUE)})
		}else{
			cl <- makeCluster(ncores)
			all.loglik <- parLapply(cl, 1:nrho, get("calculate.loglik_Y"), result=result, lower.upper=lower.upper, verbose = TRUE) 
			stopCluster(cl)
		}
		loglik <- unlist(all.loglik)
		rm(all.loglik)
	}else{
		loglik = result$loglik 
	}

	result$df <- sapply(1:nrho, function(x) sum(as.matrix(result$Theta[[x]])[upper.tri(result$Theta[[x]])] != 0 ))

	if (criterion == "ebic")
	{
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
	}
	
	if (criterion == "aic")
	{
		if(verbose)
		{
			cat("Calculating Akaike information criterion (AIC) selection at the EM convergence ...")
			flush.console()
		}
		 
		result$aic.scores	= ( - 2 * loglik ) + ( 2 * result$df )
		result$opt.index	= which.min(result$aic.scores)
		result$opt.Theta	= result$Theta[[result$opt.index]]
		result$opt.path		= abs(sign(result$opt.Theta)) - diag(rep(1,p))
		result$opt.Sigma	= result$Sigma[[result$opt.index]]
		result$opt.rho		= result$rho[result$opt.index]
	}
	
	if( sum(result$opt.path ) < p - 2 )
	{
		ind <- unlist(lapply(1:length(result$rho), function(i) sum(result$path[[i]])))
		#result$opt.index	= which(ind >  ncol(result$data))[1] #map
		result$opt.index	= which(ind >  ncol(result$data))[1] + 1 # epistasis 
		if( (is.na(result$opt.index)) || (result$opt.index > length(result$rho))) result$opt.index = length(result$rho)
		result$opt.Theta	= result$Theta[[result$opt.index]]
		result$opt.path		= abs(sign(result$opt.Theta)) - diag(rep(1,p))
		result$opt.Sigma	= result$Sigma[[result$opt.index]]
		result$opt.rho		= result$rho[result$opt.index]
	}
	
	if(verbose)
	{
		cat("done.\n")
		flush.console()
	}
	return(result)
}

