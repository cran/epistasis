epistasis = function(y, method = "approx", rho = NULL, n.rho = NULL, rho.ratio = NULL, ncores = NULL, em.iter = 10, em.tol=.001, verbose = TRUE) 
{
	if(!is.matrix(y)) y <- as.matrix(y)
	if(is.null(ncores)) ncores <- detectCores() - 1
	if(is.null(em.tol)) em.tol = 0.001
	if(is.null(em.iter)) em.iter = 10
	n = nrow(y)
	p = ncol(y)
	result = list()
		
	if( method == "gibbs")
	{
		if((is.null(rho)) && (is.null(n.rho)) ) n.rho  = 10
		if(! is.null(rho)) n.rho  = length(rho)
			
		est = vector("list", n.rho)
		for(chain in 1 : n.rho ) 
		{
		if(verbose)
		{
			m <- paste(c("Graph estimation using Gibbs sampling within EMCopula ... in progress:", floor(100 * chain/n.rho), "%"), collapse="")
			cat(m, "\r")
			flush.console()
		}
			if( chain == 1)
			{
				est[[chain]] = vector("list", n.rho)
				Theta = sparseMatrix(i = 1:ncol(y), j = 1:ncol(y), x = 1)
				est[[chain]] = Gibbs_method(y, rho=rho, n_rho=n.rho, rho_ratio=rho.ratio, Theta = Theta, ncores = ncores, chain = chain, max.elongation = em.iter, em.tol=em.tol)
			}else{
				est[[chain]] = vector("list", n.rho)
				Theta = est[[(chain - 1)]]$Theta
				Theta = as(Theta, "dgTMatrix") 
				Theta = as(Theta, "sparseMatrix")
				est[[chain]] = Gibbs_method(y, rho=rho, n_rho=n.rho, rho_ratio=rho.ratio, Theta= Theta, ncores = ncores, chain = chain, max.elongation = em.iter, em.tol=em.tol)
			}
		}
		rm(Theta)
		gc()
	}
	
	if(method == "approx")
	{

		if( !is.null(rho) ) n.rho = length(rho) 
		if( is.null(n.rho) ) n.rho = 10
		
		ini = initialize(y, rho = rho, n_rho = n.rho, rho_ratio = rho.ratio, ncores=ncores )
		rho = ini$rho
		Z	= ini$Z
		ES	= ini$ES
		lower_upper = ini$lower_upper
		
		rm(ini)
		gc()
		
		est <- vector("list", n.rho)
		for(chain in 1 : n.rho) 
		{
		if(verbose)
			{
				m <- paste(c("Graph estimation using approximation within EMCopula ... in progress:", floor(100 * chain/n.rho), "%"), collapse="")
				cat(m, "\r")
				flush.console()
			}
			est[[chain]] <- vector("list", n.rho)
			est[[(chain)]] <- approx_method(y, Z, ES=ES, rho=rho, lower_upper=lower_upper, chain = chain, ncores = ncores, em_tol=em.tol, em_iter=em.iter)
		}

		rm(lower_upper)
		gc()
	}		
	result$Theta  = vector("list", n.rho)
	result$path   = vector("list", n.rho)
	result$Sigma  = vector("list", n.rho)
	result$ES	  = vector("list", n.rho)
	result$Z	  = vector("list", n.rho)
	result$rho	  = vector()
	result$loglik = vector()
	result$data	  = y
	rm(y)
	
	for(chain in 1:n.rho)
	{
		result$Theta[[chain]]	= est[[chain]]$Theta
		if(!is.null(colnames(result$data))) colnames(result$Theta[[chain]]) = colnames(result$data)
		result$path[[chain]]	= abs(sign(result$Theta[[chain]])) - Diagonal(p)
		result$Sigma[[chain]]	= est[[chain]]$Sigma
		result$ES[[chain]]		= est[[chain]]$ES
		result$Z[[chain]]		= est[[chain]]$Z
		result$rho[chain]		= est[[chain]]$rho
		result$loglik[chain]	= est[[chain]]$loglik
	}
	rm(est)
	
	if(method == "approx")
	{
		if(verbose)
		{
			cat("Graph estimation using approximation within EMCopula .... done.             \r\n")
			flush.console()
		}
	}
	
		if(method == "gibbs")
	{
		if(verbose)
		{
			cat("Graph estimation using Gibbs sampling within EMCopula .... done.             \r\n")
			flush.console()
		}
	}
	
	class(result) = "epi"
	return(result)
}

#-----------------------------------------------------#
#   		Plot for class "epi"      	              #
#-----------------------------------------------------#

plot.epi = function( x, n.markers=NULL , ...)
{
	if(length(x$rho) == 1 ) par(mfrow = c(1, 2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
	if(length(x$rho) == 2 ) par(mfrow = c(2, 2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
	if(length(x$rho) >= 3 ) par(mfrow = c(2, ceiling(length(x$rho)/2)), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
	
	for(chain in 1:length(x$rho))
	{
			if(length(x$rho) == 1 ) image(as.matrix(x$path[[chain]]), col = gray.colors(256), xlab= "", ylab="" ,main=paste( "rho ", x$rho[chain],  sep=""))
			if(length(x$rho) == 2 ) image(as.matrix(x$path[[chain]]), col = gray.colors(256), xlab= "", ylab="" ,main=paste( "rho ", x$rho[chain],  sep=""))
			
			adj = graph.adjacency(as.matrix(x$path[[chain]]), mode="undirected", diag=FALSE)
			if(is.null(n.markers)) 
			{
					memberships = 1
					vertex.color = "red" 
			}else{
				LG = length(n.markers)
				memberships = NULL
				i = 1
				while( i <= LG){
					grp <- rep(i, n.markers[i])
					memberships = c(memberships, grp)
					i = i + 1
				}
				if(chain == 1){
					color = sample(terrain.colors(max(memberships) + 10), max(memberships))
					cl = color[memberships]
				}
				vertex.color = cl
			}
			adj$layout	= layout.fruchterman.reingold
			plot(adj, vertex.color = vertex.color, edge.color='gray40', vertex.size = 7, vertex.label = NA , vertex.label.dist = NULL)
	}
	if(length(memberships) > 1) legend("bottomright", paste("group", 1:length(n.markers)), cex=0.7, col= color, pch=rep(20,10))		
}


#-----------------------------------------------------#
#   		Summary for class "epi"        		      #
#-----------------------------------------------------#

print.epi = function(x, ...){
	cat("Estimated a graph path for", length(x$rho), "penalty term(s)" , "\n")
	cat("Number of variables: p =", ncol(x$data), "\n")
	cat("Number of sample size: n =", nrow(x$data), "\n")
	cat("Number of levels in the data: k =", length(unique(sort(as.matrix(x$data)))), "\n")
	sparsLevel <- sapply(1:length(x$rho), function(i) sum(x$Theta[[i]])/ncol(x$data)/(ncol(x$data)-1))
    cat("Sparsity level:", sparsLevel,"\n")
	cat("To plot the graph path consider plot() function \n")
	cat("To select an optimal graph consider episelect() function \n")
}