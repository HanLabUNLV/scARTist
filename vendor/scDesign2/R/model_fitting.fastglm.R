#' Function which prints a message using shell echo; useful for printing messages from inside mclapply when running in Rstudio
message_parallel <- function(...){
  system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
}

# a wrapper around mclapply. For each iteration it uses withCallingHandlers to catch and print warnings and errors; note that invokeRestart("muffleWarning") is required in order to continue FUN exection and return the result. Printing is done via message_parallel function which uses shell echo to print messages to R console
safe_mclapply <- function(X, FUN, mc.cores, stop.on.error=T, ...){
  fun <- function(x){
    res_inner <- tryCatch({
      withCallingHandlers(
        expr = {
          FUN(x, ...)
        }, 
        warning = function(e) {
          message_parallel(trimws(paste0("WARNING [element ", x,"]: ", e)))
          # this line is required to continue FUN execution after the warning
          invokeRestart("muffleWarning")
        },
        error = function(e) {
          message_parallel(trimws(paste0("ERROR [element ", x,"]: ", e)))
        }
      )},
      error = function(e){
        # error is returned gracefully; other results of this core won't be affected
        return(e)
      }
    )
    return(res_inner)
  }
  res <- mclapply(X, fun, mc.cores=mc.cores)
  failed <- sapply(res, inherits, what = "error")
  if (any(failed == T)){
    error_indices <- paste0(which(failed == T), collapse=", ")
    error_traces <- paste0(lapply(res[which(failed == T)], function(x) x$message), collapse="\n\n")
    error_message <- sprintf("Elements with following indices failed with an error: %s. Error messages: \n\n%s", 
                             error_indices,
                             error_traces)
    if (stop.on.error)
      stop(error_message)
    else
      warning(error_message, "\n\n### Errors will be ignored ###")
  }
  return(res[!failed])
}

  

#' Fit the marginal distributions for each row of a count matrix
#'
#' @keywords internal
#'
#' @param x            A matrix of shape p by n that contains count values.
#' @param marginal     Specification of the types of marginal distribution.
#'                     Default value is 'auto_choose' which chooses between ZINB, NB, ZIP
#'                     and Poisson by a likelihood ratio test (lrt) and whether there is
#'                     underdispersion.
#'                     'zinb' will fit the ZINB model. If there is underdispersion, it
#'                     will choose between ZIP and Poisson by a lrt. Otherwise, it will try to
#'                     fit the ZINB model. If in this case, there is no zero at all or an error
#'                     occurs, it will fit an NB model instead.
#'                     'nb' fits the NB model that chooses between NB and Poisson depending
#'                     on whether there is underdispersion.
#'                     'poisson' simply fits the Poisson model.
#' @param pval_cutoff  Cutoff of p-value of the lrt that determines whether
#'                     there is zero inflation.
#' @param epsilon      Threshold value for preventing the transformed quantile
#'                     to collapse to 0 or 1.
#' @param jitter       Logical, whether a random projection should be performed in the
#'                     distributional transform.
#' @param DT           Logical, whether distributional transformed should be performed.
#'                     If set to FALSE, the returned object u will be NULL.
#'
#' @return             a list with the following components:
#'\describe{
#'  \item{params}{a matrix of shape p by 3. The values of each column are: the ZI proportion,
#'  the dispersion parameter (for Poisson, it's Inf), and the mean parameter.}
#'  \item{u}{NULL or a matrix of the same shape as x, which records the transformed quantiles,
#'  by DT.}
#'}
fit_marginals_orig <- function(x, marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                          pval_cutoff = 0.05, epsilon = 1e-5,
                          jitter = TRUE, DT = TRUE){
  p <- nrow(x)
  n <- ncol(x)

  marginal <- match.arg(marginal)
  if(marginal == 'auto_choose'){
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v){
        mle_Poisson <- glm(gene ~ 1, family = poisson, method="fastglm")
        tryCatch({
          mle_ZIP <- zeroinfl(gene ~ 1|1, dist = 'poisson')
          chisq_val <- 2 * (logLik(mle_ZIP) - logLik(mle_Poisson))
          pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
          if(pvalue < pval_cutoff)
            c(plogis(mle_ZIP$coefficients$zero), Inf, exp(mle_ZIP$coefficients$count))
          else
            c(0.0, Inf, m)
        },
        error = function(cond){
          c(0.0, Inf, m)})
      }else{
        mle_NB <- glm.nb(gene ~ 1, method="fastglm")
        if(min(gene) > 0)
          c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
        else
          tryCatch({
            mle_ZINB <- zeroinfl(gene ~ 1|1, dist = 'negbin')
            chisq_val <- 2 * (logLik(mle_ZINB) - logLik(mle_NB))
            pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
            if(pvalue < pval_cutoff)
              c(plogis(mle_ZINB$coefficients$zero), mle_ZINB$theta, exp(mle_ZINB$coefficients$count))
            else
              c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          },
          error = function(cond){
            c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          })
      }
    }))
  }else if(marginal == 'zinb'){
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v)
      {
        mle_Poisson <- glm(gene ~ 1, family = poisson, method="fastglm")
        tryCatch({
          mle_ZIP <- zeroinfl(gene ~ 1|1, dist = 'poisson')
          chisq_val <- 2 * (logLik(mle_ZIP) - logLik(mle_Poisson))
          pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
          if(pvalue < pval_cutoff)
            c(plogis(mle_ZIP$coefficients$zero), Inf, exp(mle_ZIP$coefficients$count))
          else
            c(0.0, Inf, m)
        },
        error = function(cond){
          c(0.0, Inf, m)})
      }
      else
      {
        if(min(gene) > 0)
        {
          mle_NB <- glm.nb(gene ~ 1, method="fastglm")
          c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
        }
        else
          tryCatch({
            mle_ZINB <- zeroinfl(gene ~ 1|1, dist = 'negbin')
            c(plogis(mle_ZINB$coefficients$zero), mle_ZINB$theta, exp(mle_ZINB$coefficients$count))
          },
          error = function(cond){
            mle_NB <- glm.nb(gene ~ 1, method="fastglm")
            c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          })
      }
    }))
  }else if(marginal == 'nb'){
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v){
        c(0.0, Inf, m)
      }else{
        mle_NB <- glm.nb(gene ~ 1, method="fastglm")
        c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
      }
    }))
  }else if(marginal == 'poisson'){
    params <- t(apply(x, 1, function(gene){
      c(0.0, Inf, mean(gene))
    }))
  }

  if(DT){
    u <- t(sapply(1:p, function(iter){
      param <- params[iter, ]
      gene <- unlist(x[iter,])
      prob0 <- param[1]
      u1 <- prob0 + (1 - prob0) * pnbinom(gene, size = param[2], mu = param[3])
      u2 <- (prob0 + (1 - prob0) * pnbinom(gene - 1, size = param[2], mu = param[3])) *
        as.integer(gene > 0)
      if(jitter)
        v <- runif(n)
      else
        v <- rep(0.5, n)
      r <- u1 * v + u2 * (1 - v)
      idx_adjust <- which(1-r < epsilon)
      r[idx_adjust] <- r[idx_adjust] - epsilon
      idx_adjust <- which(r < epsilon)
      r[idx_adjust] <- r[idx_adjust] + epsilon

      r
    }))
  }else{
    u <- NULL
  }

  return(list(params = params, u = u))
}

#' Fit the marginal distributions for each row of a count matrix
#'
#' @keywords internal
#'
#' @param x            A matrix of shape p by n that contains count values.
#'
#' @return             a list with the following components:
#'\describe{
#'  \item{params}{a matrix of shape p by 4. The values of each column are: the ZI proportion,
#'  the dispersion parameter (for Poisson, it's Inf), and the mean parameter, and the logLik.}
#'}
fit_marginals_nb_zinb <- function(x, marginals, ncores = 1) { 
  p <- nrow(x)
  n <- ncol(x)

  if (marginals == 'nb') {
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v){
        c(0.0, Inf, m, -Inf)
      }else{
        mle_NB <- glm.nb(gene ~ 1, method="fastglm")
        c(0.0, mle_NB$theta, exp(mle_NB$coefficients), logLik(mle_NB))
      }
    }))
  } else if (marginals == 'zinb') {
    params <- safe_mclapply(1:nrow(x), function(iter){
      gene = x[iter,]
      m <- mean(gene)
      v <- var(gene)
      
      tryCatch({
        mle_ZINB <- zeroinfl(gene ~ 1|1, dist = 'negbin')
        c(plogis(mle_ZINB$coefficients$zero), mle_ZINB$theta, exp(mle_ZINB$coefficients$count), logLik(mle_ZINB))
      },
      error = function(cond){
        c(0.0, Inf, m, -Inf)
      })
    }, mc.cores = ncores)
    params = matrix(unlist(params), ncol = 4, byrow = TRUE)
    rownames(params) = rownames(x) 
  }
  return(list(params = params))
}



#' Fit the marginal distributions for each row of a count matrix
#'
#' @keywords internal
#'
#' @param x            A matrix of shape p by n that contains count values.
#' @param marginal     Specification of the types of marginal distribution.
#'                     Default value is 'auto_choose' which chooses between ZINB, NB, ZIP
#'                     and Poisson by a likelihood ratio test (lrt) and whether there is
#'                     underdispersion.
#'                     'zinb' will fit the ZINB model. If there is underdispersion, it
#'                     will choose between ZIP and Poisson by a lrt. Otherwise, it will try to
#'                     fit the ZINB model. If in this case, there is no zero at all or an error
#'                     occurs, it will fit an NB model instead.
#'                     'nb' fits the NB model that chooses between NB and Poisson depending
#'                     on whether there is underdispersion.
#'                     'poisson' simply fits the Poisson model.
#' @param pval_cutoff  Cutoff of p-value of the lrt that determines whether
#'                     there is zero inflation.
#' @param epsilon      Threshold value for preventing the transformed quantile
#'                     to collapse to 0 or 1.
#' @param jitter       Logical, whether a random projection should be performed in the
#'                     distributional transform.
#' @param DT           Logical, whether distributional transformed should be performed.
#'                     If set to FALSE, the returned object u will be NULL.
#'
#' @return             a list with the following components:
#'\describe{
#'  \item{params}{a matrix of shape p by 3. The values of each column are: the ZI proportion,
#'  the dispersion parameter (for Poisson, it's Inf), and the mean parameter.}
#'  \item{u}{NULL or a matrix of the same shape as x, which records the transformed quantiles,
#'  by DT.}
#'}
fit_marginals <- function(x, marginal = 'nb',
                          pval_cutoff = 0.05, epsilon = 1e-5,
                          jitter = TRUE, DT = TRUE, ncores = 1){
  p <- nrow(x)
  n <- ncol(x)
  dist = c('nb', 'zinb')
  parallel_param <- safe_mclapply(1:length(dist), function(iter){
    fit_marginals_nb_zinb(x, marginals = dist[iter], ncores=ncores)
  }, mc.cores = 4)
  params_nb <- parallel_param[[1]]$params
  params_zinb <- parallel_param[[2]]$params
  nb_logLik = params_nb[,4]
  zinb_logLik = params_zinb[,4]

  chisq_val <- 2 * (zinb_logLik - nb_logLik)
  pvalue <- as.numeric(1 - pchisq(chisq_val, 1)) 
  params = matrix(,nrow=p, ncol=3)
  params[pvalue < pval_cutoff,] = params_zinb[pvalue < pval_cutoff,1:3]
  params[pvalue >= pval_cutoff,] = params_nb[pvalue >= pval_cutoff,1:3]

  if(DT){
    u <- t(sapply(1:p, function(iter){
      param <- params[iter, ]
      gene <- unlist(x[iter,])
      prob0 <- param[1]
      u1 <- prob0 + (1 - prob0) * pnbinom(gene, size = param[2], mu = param[3])
      u2 <- (prob0 + (1 - prob0) * pnbinom(gene - 1, size = param[2], mu = param[3])) *
        as.integer(gene > 0)
      if(jitter)
        v <- runif(n)
      else
        v <- rep(0.5, n)
      r <- u1 * v + u2 * (1 - v)
      idx_adjust <- which(1-r < epsilon)
      r[idx_adjust] <- r[idx_adjust] - epsilon
      idx_adjust <- which(r < epsilon)
      r[idx_adjust] <- r[idx_adjust] + epsilon

      r
    }))
  }else{
    u <- NULL
  }

  return(list(params = params, u = u))
}


#' Fit a Gaussian copula model for a count matrix of a single cell type
#'
#' @inheritParams fit_marginals
#' @param zp_cutoff            The maximum propotion of zero allowed for a gene to be included
#'                             in the joint copula model.
#' @param min_non_zero_num     The minimum number of non-zero values required for a gene to be
#'                             fitted a marginal model.
#'
#' @return The genes of \code{x} will be partitioned into three groups. The first group contains
#' genes whose zero proportion is less than \code{zp_cutoff}. The second group contains genes
#' whose zero proportion is greater than \code{zp_cutoff} but still contains at least
#' \code{min_non_zero_num} non-zero values. The third and last group contains the rest of the
#' genes. For the first group, a joint Gaussian copula model will be fitted. For the second group,
#' only the marginal distribution of each gene will be fitted. For the last group, no model will
#' be fitted and only the index of these genes will be recorded. A list that contains the above
#' fitted model will be returned that contains the following components.
#' \describe{
#' \item{cov_mat}{The fitted covariance (or equivalently in this case, correlation) matrix of the
#' Gaussin copula model.}
#' \item{marginal_param1}{A matrix of the parameters for the marginal distributions of genes in
#' group one.}
#' \item{marginal_param2}{A matrix of the parameters for the marginal distributions of genes in
#' group two.}
#' \item{gene_sel1}{A numeric vector of the row indices of the genes in group one.}
#' \item{gene_sel2}{A numeric vector of the row indices of the genes in group two.}
#' \item{gene_sel3}{A numeric vector of the row indices of the genes in group three.}
#' \item{zp_cutoff}{Same as the input.}
#' \item{min_non_zero_num}{Same as the input.}
#' \item{sim_method}{A character string that says 'copula'. To be distinguished with the
#' (w/o copula) model.}
#' }
#'
# @export
fit_Gaussian_copula <- function(x, marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                                jitter = TRUE, zp_cutoff = 0.8,
                                min_nonzero_num = 2, ncores = 1){
  marginal <- match.arg(marginal)
  n <- ncol(x)
  p <- nrow(x)

  gene_zero_prop <- apply(x, 1, function(y){
    sum(y < 1e-5) / n
  })

  gene_sel1 <- which(gene_zero_prop < zp_cutoff)
  gene_sel2 <- which(gene_zero_prop < 1.0 - min_nonzero_num/n &
                       gene_zero_prop >= zp_cutoff)
  gene_sel3 <- (1:p)[-c(gene_sel1, gene_sel2)]

  if(length(gene_sel1) > 0){
    marginal_result1 <- fit_marginals(x[gene_sel1, , drop = FALSE], marginal, jitter = jitter, DT = TRUE, ncores=ncores)
    quantile_normal <- qnorm(marginal_result1$u)
    cov_mat <- cor(t(quantile_normal))
  }else{
    cov_mat = NULL
    marginal_result1 = NULL
  }

  if(length(gene_sel2) > 0){
    marginal_result2 <- fit_marginals(x[gene_sel2, , drop = FALSE], marginal, DT = FALSE, ncores=ncores)
  }else{
    marginal_result2 = NULL
  }
  return(list(cov_mat = cov_mat, marginal_param1 = marginal_result1$params,
              marginal_param2 = marginal_result2$params,
              gene_sel1 = gene_sel1, gene_sel2 = gene_sel2, gene_sel3 = gene_sel3,
              zp_cutoff = zp_cutoff, min_nonzero_num = min_nonzero_num,
              sim_method = 'copula', n_cell = n, n_read = sum(x)))
}



#' Fit a (w/o copula) model for a count matrix of a single cell type
#'
#' This function only fits the marginal distribution for each gene.
#'
#' @inheritParams fit_marginals
#' @inheritParams fit_Gaussian_copula
#' @return The genes of \code{x} will be partitioned into two groups. The first group contains
#' genes with at least \code{min_non_zero_num} non-zero values. The second group contains the
#' other genes. For the first group, the marginal distribution of each gene will be fitted.
#' For the last group, no model will be fitted and only the index of these genes will be
#' recorded. A list that contains the above fitted model will be returned that contains the
#' following components.
#' \describe{
#' \item{marginal_param1}{A matrix of the parameters for the marginal distributions of genes in
#' group one.}
#' \item{gene_sel1}{A numeric vector of the row indices of the genes in group one.}
#' \item{gene_sel2}{A numeric vector of the row indices of the genes in group two.}
#' \item{min_non_zero_num}{Same as the input.}
#' \item{sim_method}{A character string that says 'ind', short for 'independent'. To be
#' distinguished with the copula model.}
#' }
#'
# @export
fit_wo_copula <- function(x, marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                          jitter = TRUE, min_nonzero_num = 2){
  marginal <- match.arg(marginal)
  n <- ncol(x)
  p <- nrow(x)

  gene_zero_prop <- apply(x, 1, function(y){
    sum(y < 1e-5) / n
  })

  gene_sel1 <- which(gene_zero_prop < 1.0 - min_nonzero_num/n)
  gene_sel2 <- (1:p)[-gene_sel1]

  if(length(gene_sel1) > 0){
    marginal_result1 <- fit_marginals(x[gene_sel1, ], marginal, jitter = jitter, DT = FALSE)
  }else{
    marginal_result1 = NULL
  }

  return(list(marginal_param1 = marginal_result1$params,
              gene_sel1 = gene_sel1, gene_sel2 = gene_sel2,
              min_nonzero_num = min_nonzero_num, sim_method = 'ind',
              n_cell = n, n_read = sum(x)))
}



#' Fit models for a count matrix
#'
#' @param data_mat      A matrix of shape p by n that contains count values. Each of its
#'                      column names should be the cell type names of that cell. Its column
#'                      names should also match \code{cell_type_sel}.
#' @param cell_type_sel A character vector that contains the selected cell types for which a
#'                      model will be fitted.
#' @param sim_method    Specification of the type of model.
#'                      Default value is 'copula', which selects the copula model.
#'                      'ind' will select the (w/o copula) model.
#' @inheritParams fit_Gaussian_copula
#' @param ncores        A numeric value that indicates the number of parallel cores for model
#'                      fitting. One core for each cell type.
#' @return A list with the same length as \code{cell_type_sel} that contains the fitted model
#' as each of its element.
#'
#' @export
fit_model_scDesign2 <- function(data_mat, cell_type_sel, sim_method = c('copula', 'ind'),
                                marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                                jitter = TRUE, zp_cutoff = 0.8,
                                min_nonzero_num = 2, ncores = 1){
  sim_method <- match.arg(sim_method)
  marginal <- match.arg(marginal)

  if(sum(abs(data_mat - round(data_mat))) > 1e-5){
    warning('The entries in the input matrix are not integers. Rounding is performed.')
    data_mat <- round(data_mat)
  }

  if(sim_method == 'copula'){
    param <- safe_mclapply(1:length(cell_type_sel), function(iter){
      X <- data_mat[, colnames(data_mat) == cell_type_sel[iter]]
      fit_Gaussian_copula(X, marginal,
                          jitter = jitter, zp_cutoff = zp_cutoff,
                          min_nonzero_num = min_nonzero_num, ncores=ncores)
    }, mc.cores = ncores)
  }else if(sim_method == 'ind'){
    param <- safe_mclapply(1:length(cell_type_sel), function(iter){
      fit_wo_copula(data_mat[, colnames(data_mat) == cell_type_sel[iter]], marginal,
                    jitter = jitter,
                    min_nonzero_num = min_nonzero_num)
    }, mc.cores = ncores)
  }

  names(param) <- cell_type_sel
  param
}


