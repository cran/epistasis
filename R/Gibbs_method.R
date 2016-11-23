
Gibbs_method = function(y, rho = NULL, n_rho = NULL, rho_ratio = NULL, Theta=NULL, ncores = 4, chain = 1, max.elongation = 10, em.tol=0.001) 
{
  p <- ncol(y)
  n <- nrow(y)
  lower.upper = lower.upper(y)

  if(is.null(rho))
  {
	if(is.null(n_rho)) n_rho = 10
	if(is.null(rho_ratio)) rho_ratio = 0.2
	cr = cor(y, method="spearman")
	cr[is.na(cr)] <- 0
	S  = cr - diag(p)
	rho_max = max(max(S),-min(S))
	if(rho_max == 1) rho_max = 0.7
	rho_min = rho_ratio * rho_max
	rho = exp(seq(log(rho_max), log(rho_min), length = n_rho))
	rm(cr, S, rho_max, rho_min, rho_ratio)
  }

  #Main code
  Gibbs.method <- calculate_EM_Gibbs(chain, y, rho = rho[chain], Theta=Theta, lower.upper = lower.upper, em.tol = em.tol, em.iter = max.elongation, gibbs.iter = 500, mc.iter = 1000, ncores = ncores)
  invisible(return(Gibbs.method))
}  