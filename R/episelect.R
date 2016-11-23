episelect = function(epi.object, criteria = NULL, ebic.gamma = 0.5, loglik_Y = FALSE, ncores = NULL)
{
	if(is.null(ncores)) ncores <- detectCores() - 1
	if(is.null(criteria)) criteria <- "ebic"
	if(!loglik_Y) ncores = 1
	if(loglik_Y)
	{
		lower.upper = lower.upper(epi.object$data)
	}else{
		lower.upper = NULL
	}
	
	sel	= model.selection( epi.object, criterion = criteria, lower.upper=lower.upper, ebic.gamma=ebic.gamma, ncores = ncores, loglik_Y = loglik_Y, verbose = TRUE)
	
	class(sel) = "episelect"
	return(sel)
}

#-----------------------------------------------------#
#   		Plot for class "episelect"      		  #
#-----------------------------------------------------#
readkey <- function()
{
    cat ("Press [enter] to continue")
    line <- readline()
}

plot.episelect = function(x, n.markers = NULL, vertex.size = NULL, vertex.label = FALSE, ...)
{	
		if(! vertex.label) {
			vertex.label = NA
		}else{
			if(!is.null(colnames(x$opt.path) )) 
			{
				vertex.label = colnames(x$opt.path)
			}else{
				vertex.label= NA
			}
		}
		if(is.null(vertex.size)) vertex.size = 7
		
    	adj = graph.adjacency(as.matrix(x$opt.path), mode="undirected", diag=FALSE)
		if(is.null(n.markers)) 
		{
			memberships = 1
			vertex.color = "red"
		}else{
			LG = length(n.markers)
			memberships = NULL
			i = 1
			while( i <= LG)
			{
				grp <- rep(i, n.markers[i])
				memberships = c(memberships, grp)
				i = i + 1
			}
			color <- sample(terrain.colors(max(memberships)+10), max(memberships))
			vertex.color = color[memberships]
		}
		adj$layout	= layout.fruchterman.reingold 
	
		plot(adj, vertex.color= vertex.color , edge.color='gray40', vertex.size = vertex.size, vertex.label = vertex.label, vertex.label.dist = 0, main= "Selected graph")	 	  
		if(length(memberships) > 1) legend("bottomright", paste("group", 1:length(n.markers)), cex=0.7, col= color, pch=rep(20,10))
		readkey()
		image(as.matrix(x$opt.path), xlab= "", ylab="", col = gray.colors(256) ,main="Selected graph" )
}


