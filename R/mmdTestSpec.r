# This is an R port for the MMD spectral variant code written by Arthur Gretton (2008)

rbf.dot = function(patterns1, patterns2, deg) {
  # Radial basis function inner product
  
  # Pattern input format : [pattern1 ; pattern2 ; ...]
  # Output : [p11*p21 p11*p22 ... ; p12*p21 ... ]
  # deg is kernel size

  size1 = dim(patterns1)
  size2 = dim(patterns2)

  # new vectorised version

  G = rowSums((patterns1 * patterns1))
  H = rowSums((patterns2 * patterns2))
    
  Q = repmat(G, 1, size2[1])
  R = repmat(H, size1[1], 1)
    
  H = Q + R - 2 * patterns1 %*% t(patterns2)

  H = exp(-H/2/deg ^ 2)
  
  return (H)
}
           
mmdTestSpec = function(X, Y, alpha, params) {
  # This function implements the MMD two-sample test using the Gram
  # matrix spectrum to estimate the coefficients of the infinite sum
  # of chi^2 (which constitutes the null distibution)
  
  # Inputs: 
  #   X contains dx columns, m rows. Each row is an i.i.d sample
  #   Y contains dy columns, m rows. Each row is an i.i.d sample
  #   alpha is the level of the test
  #   params is a list with the following elements
  #     sig is kernel size. If -1, use median-distance heuristic.
  #     numEigs is number of eigenvalues used to compute the null distribution. If it is -1, then use 2*m - 2
  #     numNullSamp is number of samples from null spectral MMD distribution used to estimate CDF
  
  # Outputs: 
  #   thresh: test threshold for level alpha test
  #   testStat: test statistic: m * MMD_b (biased)
  
  # Set kernel size to median distance between points, if no kernel specified
  
  m = nrow(X)
  
  if (params$numEigs == -1) {
    params$numEigs = 2 * m - 2
  }
  
  if (params$sig == -1) {
    # Set kernel size to median distance between points in aggregate sample
    Z = rbind(X, Y) # pool all samples together

    size1 = nrow(Z)
    if (size1 > 100) {
      Zmed = as.matrix(Z[1:100, ])
      size1 = 100
    } else {
      Zmed = Z
    }
    
    G = rowSums(Zmed * Zmed)
    Q = repmat(G, 1, size1)
    R = repmat(G, size1, 1)
    
    dists = Q + R - 2 * Zmed %*% t(Zmed)
    dists = dists - tril(dists)
    dists = reshape(as.matrix(dists), size1 ^ 2, 1)
    params$sig = sqrt(0.5 * median(dists[dists > 0])) #rbf_dot has factor two in kernel
  }
  
  K  = rbf.dot(X, X, params$sig)
  L  = rbf.dot(Y, Y, params$sig)
  KL = rbf.dot(X, Y, params$sig)
  
  # MMD statistic. Here we use biased v-statistic. 
  # NOTE: this is m * MMD_b
  testStat = 1/m * sum(sum(K + L - KL - t(KL)))
  
  # Draw samples from null distribution
  Kz = rbind(cbind(K, KL), cbind(t(KL), L))
  H = diag(2 * m) - 1 / (2 * m) * ones(2 * m, 2 * m)
  Kz = H %*% Kz %*% H
  
  opts = list(retvec = F, 
      ncv = min(nrow(Kz), max(2 * params$numEigs + 1, 20)) # (no idea why the documentation says this is the default but it isn't)
    )
  
  kEigs = eigs_sym(Kz, params$numEigs, opts = opts) # NOTE: this retains only largest magnitude eigenvalues
  
  # Empirical eigenvalues scaled by 1/2/m: see p. 2 Shawe-Taylor et al. (2005)
  kEigs = 1/2/m * abs(kEigs$values) 
  numEigs = length(kEigs)  
  
  nullSampMMD = zeros(1, params$numNullSamp)
  
  for (whichSamp in 1:params$numNullSamp) {
    nullSampMMD[whichSamp] = 2 * sum(kEigs * rnorm(params$numEigs) ^ 2)
  }
  
  nullSampMMD = sort(nullSampMMD)
  thresh = nullSampMMD[round((1 - alpha) * params$numNullSamp)]
  
  return (list(testStat = testStat, thresh = thresh, params = params))
}
