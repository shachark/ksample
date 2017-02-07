# R wrapper for the DynSlicing implementation by Chao Ye and Bo Jiang.

# In the C functions:
# R_y is assumed to be sorted in ascending order and
# contain quantiles of F0 at the raw observations.
# R_x is assumed to hold sample identifiers in [0, ..., K-1] ordered according
# to ascending pooled-sample observation order.

# TODO: GOF tests interface
#SEXP R_C_dynslicing_1eqp(SEXP R_y, SEXP R_lambda)
#SEXP R_C_dynslicing_one(SEXP R_y, SEXP R_lambda, SEXP R_alpha)

# TODO simplified (sped-up) K-sample test
#SEXP R_C_dynslicing_keqp(SEXP R_x, SEXP R_dim, SEXP R_lambda, SEXP R_slices_wanted)

dynSlicing.ks.test = function(obs, g, lambda = NULL, alpha = 0.05, nr.perms = 0, perm.stats.wanted = F, slices.wanted = F) {
  if (!is.numeric(obs) && !is.ordered(obs)) {
    stop('obs is expected to be a numeric or ordered vector')
  }
  if (!is.numeric(g) && !is.integer(g) && !is.factor(g) && !is.logical(g)) {
    stop('g is expected to be a numeric, integer, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  if (is.factor(g)) {
    g = as.numeric(levels(g))[g]
  } else if (is.logical(g)) {
    g = as.numeric(g)
  }
  if (any(g != round(g))) {
    stop('g is expected to be a numeric, integer, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  
  K = as.integer(max(g) + 1)
  if (K <= 1) {
    stop('g is expected to be a numeric, integer, logical, or factor vector with values in {0, 1, ... K-1}, with K >= 2')
  }
  
  if (is.null(lambda)) {
    if (K == 2 && (sum(g == 0) == sum(g == 1))) {
      # This may be a crude approximation, should load it from their table:
      # http://www.people.fas.harvard.edu/~junliu/DS/lambda-table.html
      lambda = log(1/alpha) / log(length(g)) + 0.45
      #cat('Using lambda =', lambda, '\n')
    } else {
      stop('Did not implement yet a general automatic derivation of lambda')
    }
  } else {
    lambda = as.double(lambda)
  }
  
  slices_wanted = as.integer(slices.wanted)
  
  ret = list()

  g = as.integer(g)
  x = g[order(obs)]
  res = .Call('R_C_dynslicing_k', x, K, lambda, slices_wanted)
  ret$obs.stat = res[1]

  if (perm.stats.wanted) {
    ret$perm.stats = rep(NA, nr.perms)
  }
  
  # NOTE that in the author's wrapper implementation they are a little more liberal with p-values
  # and do not count '>=', rather only '>'.
  
  if (nr.perms > 0) {
    cnt = 1
    
    for (permi in 1:nr.perms) {
      x = sample(g)
      res = .Call('R_C_dynslicing_k', x, K, lambda, 0L)
      cnt = cnt + (res[1] >= ret$obs.stat)
      if (perm.stats.wanted) {
        ret$perm.stats[permi] = res[1]
      }
    }

    ret$p.value = cnt / (nr.perms + 1)
  }

  return (ret)
}
