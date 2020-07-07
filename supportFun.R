if (!("fdapace" %in% rownames(installed.packages()))) {
  install.packages("fdapace")
}

#' type = 
#' 100: \int f(t)dt
#' 201: \int f(s,t)dt
#' 211: \int f(s,t)g(t)dt
#' 222: \int f(s,w)g(w,t)dw

integral = function(f, g = NULL, domain, type){
  if (type == 100){
    f = as.vector(f)
    len.f = length(f)
    result = sum((f[-1] + f[-len.f]) * diff(domain))/2
    return(result)
  }
  
  if (type == 201){
    nrow.f = dim(f)[1]
    ncol.f = dim(f)[2]
    gap.mat = matrix(diff(domain), nrow = nrow.f, ncol = length(domain)-1, byrow = T)
    result = matrix((rowSums(f[, -1] * gap.mat) + rowSums(f[, -ncol.f] * gap.mat))/2, nrow = nrow.f, ncol = 1)
    return(result)
  }
  
  if (type == 211){
    nrow.f = dim(f)[1]
    ncol.f = dim(f)[2]
    gap.mat = matrix(diff(domain), nrow = nrow.f, ncol = length(domain)-1, byrow = T)
    result = ((f[, -1] * gap.mat) %*% as.matrix(g[-1]) + (f[, -ncol.f] * gap.mat) %*% as.matrix(g[-ncol.f]))/2
    return(as.matrix(result))
  }
  
  if (type == 222){
    nrow.f = dim(f)[1]
    ncol.f = dim(f)[2]
    gap.mat = matrix(diff(domain), nrow = nrow.f, ncol = length(domain)-1, byrow = T)
    result = ((f[, -1] * gap.mat) %*% g[-1, ] + (f[, -ncol.f] * gap.mat) %*% g[-ncol.f, ])/2
    return(result)
  }
}

# Gram-Schmidt orthonormalization w.r.t. ker
GSortho = function(basis.origi, ker, domain.x){
  
  if (is.vector(basis.origi)){
    integ1 = integral(ker, basis.origi, domain.x, type = 211)
    integ2 = integral(integ1 * basis.origi, domain = domain.x, type = 100)
    if (integ2 > .Machine$double.eps){
      basis.ortho = basis.origi / integ2^.5
    }else{
      basis.ortho = 0
    }
    return(as.matrix(basis.ortho))
  }else{
    basis.ortho = array(NA, dim = dim(basis.origi))
  }
  
  p = dim(basis.origi)[2]
  for (i in 1:p){
    psi.tmp = as.matrix(basis.origi[, i])
    if (i > 1){
      for (j in 1:(i-1)){
        integ1 = integral(ker, psi.tmp, domain.x, type = 211)
        integ2 = integral(integ1 * as.matrix(basis.ortho[, j]), domain = domain.x, type = 100)
        psi.tmp = psi.tmp - as.matrix(basis.ortho[, j]) * integ2
      }
    }
    integ1 = integral(ker, psi.tmp, domain.x, type = 211)
    integ2 = integral(integ1 * psi.tmp, domain = domain.x, type = 100)
    if (integ2 > .Machine$double.eps){
      basis.ortho[, i] = psi.tmp / integ2^.5
    }else{
      basis.ortho[, i] = 0
    }
  }
  return(basis.ortho)
}

GCVLwls1D = function (yy, tt, win = NULL, kernel, npoly, nder, dataType) {
  if (is.vector(yy) && is.vector(tt) && !is.list(yy) && !is.list(tt)) {
    if (length(tt) < 21) 
      stop("You are trying to use a local linear weight smoother in a vector with less than 21 values.\n")
    myPartition = c(1:10, sample(10, length(tt) - 10, replace = TRUE))
    yy = split(yy, myPartition)
    tt = split(tt, myPartition)
    dataType = "Sparse"
  }
  t = unlist(tt)
  y = unlist(yy)[order(t)]
  t = sort(t)
  N = length(t)
  r = t[N] - t[1]
  
  if (dataType == "Sparse") {
    dstar = fdapace:::Minb(t, npoly + 2)
    if (dstar > r * 0.25) {
      dstar = dstar * 0.75
      warning(c("The min bandwidth choice is too big, reduce to ", dstar, "!\n"))
    }
    h0 = 2.5 * dstar
  } else if (dataType == "DenseWithMV")
    h0 = 2 * fdapace:::Minb(t, npoly + 1)
  else
    h0 = 1.5 * fdapace:::Minb(t, npoly + 1)
  if (is.nan(h0)) {
    if (kernel == "gauss")
      h0 = 0.2 * r
    else
      stop("The data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n")
  }
  h0 = min(h0, r)
  q = (r/(4 * h0))^(1/9)
  bwCandidates = sort(q^(0:9) * h0)
  
  idx = pracma::uniq(t)$n
  k0_candidates = list(quar = 0.9375, epan = 0.75, rect = 0.5, 
    gausvar = 0.498677, gausvar1 = 0.598413, gausvar2 = 0.298415, other = 0.398942)
  if (any(names(k0_candidates) == kernel))
    k0 = as.numeric(k0_candidates[kernel])
  else
    k0 = as.numeric(k0_candidates$other)
  gcvScores = c()
  for (i in 1:length(bwCandidates)) {
    newmu = fdapace:::Lwls1D(bwCandidates[i], kernel_type = kernel, 
      npoly = npoly, nder = nder, xin = t, yin = y, win = win, xout = sort(unique(t)))[idx]
    cvsum = sum((newmu - y)^2)
    gcvScores[i] = cvsum/(1 - (r * k0)/(N * bwCandidates[i]))^2
  }
  if (all((is.infinite(gcvScores)))) {
    bwCandidates = seq(max(bwCandidates), r, length.out = length(bwCandidates))
    for (i in 1:length(bwCandidates)) {
      newmu = fdapace:::Lwls1D(bwCandidates[i], kernel_type = kernel, 
        npoly = npoly, nder = nder, xin = t, yin = y, win = win, xout = sort(unique(t)))[idx]
      cvsum = sum((newmu - y)^2)
      gcvScores[i] = cvsum/(1 - (r * k0)/(N * bwCandidates[i]))^2
    }
  }
  if (all((is.infinite(gcvScores)))) 
    stop("The data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n")
  bInd = which(gcvScores == min(gcvScores))
  bScr = gcvScores[bInd][1]
  bOpt = max(bwCandidates[bInd])
  if (bOpt == r) 
    warning("data is too sparse, optimal bandwidth includes all the data!You may want to change to Gaussian kernel!\n")
  bOptList = list(bOpt = bOpt, bScore = bScr)
  
  return(bOptList)
}

getGCVscores = function (bw, kern, xin, yin, win = NULL, regGrid, verbose = FALSE) {
  if (is.null(win)) 
    win = rep(1, length(yin))
  fit = tryCatch(fdapace:::Lwls2D(bw, kern, xin = xin, yin = yin, win = win, xout1 = regGrid, xout2 = regGrid), 
    error = function(err) {
      if (verbose)
          message("Invalid bandwidth. Try enlarging the window size.\n")
      return(Inf)
    }
  )
  if (is.infinite(fit[1])) 
    return(Inf)
  if (any(is.nan(fit))) 
    return(Inf)
  obsFit = fdapace:::interp2lin(regGrid, regGrid, fit, xin[, 1], xin[, 2])
  res = sum((yin - obsFit)^2 * win)
  k0 = fdapace:::KernelAt0(kern)
  N = sum(win)
  r = diff(range(xin[, 1]))
  bottom = max(1 - 3 * (1/N) * (r * k0/bw)^2, 0)
  GCV = res/bottom^2
  return(GCV)
}

GCVLwls2D = function (t, obsGrid, regGrid, kern = 'epan', rcov, verbose = FALSE, win = NULL, 
                      dataType = 'Sparse', error = TRUE, h0 = NULL, ngrid = NULL) {
  r = diff(range(obsGrid)) * sqrt(2)
  
  minBW = fdapace:::GetMinb(t, dataType = 'Sparse', obsGrid = obsGrid)
  if (missing(h0)) 
    h0 = minBW
  if (is.null(h0)) 
    stop("The data is too sparse and no suitable bandwidth can be found! Try Gaussian Kernel instead!\n")
  h0 = min(h0, r/4)
  if (h0 < r/4)
    q = (r/(4 * h0))^(1/9)
  else if (h0 < r/2)
    q = (r/(2 * h0))^(1/9)
  else if (h0 < r) 
    q = (r/h0)^(1/9)
  else 
    stop("Data is too sparse. The minimal bandwidth is the range of data.")
  bw = (q^seq(0, 9, length.out = 10)) * h0
  opth = h0
  
  leave = FALSE
  iter = 0
  maxIter = 1
  minBWInvalid = FALSE
  while (!leave && iter < maxIter) {
    if (minBWInvalid) 
      minBW = bw[1]
    Scores = matrix(Inf, nrow = length(bw), ncol = 2)
    colnames(Scores) = c("SUM", "SD")
    for (i in rev(seq_along(bw))) {
      Scores[i, "SUM"] = getGCVscores(bw[i], kern = optns$kernel, xin = rcov$tPairs, yin = rcov$cxxn, win, regGrid, verbose = optns$verbose)
      if (is.infinite(Scores[i, "SUM"])) {
        minBWInvalid = TRUE
        if (i < length(bw)) {
          if (minBWInvalid) {
            minBW = bw[i + 1]
            minBWInvalid = FALSE
          }
        }
        break
      }
    }
    if (is.infinite(min(Scores[, "SUM"]))) {
      opth = max(bw)
      optgcv = Inf
    } else {
      ind = which.min(Scores[, "SUM"])
      opth = bw[ind]
      optgcv = Scores[ind, "SUM"]
    }
    if (opth >= r - .Machine$double.eps) {
      minBW = r
      leave = TRUE
      stop("Data is too sparse. The optimal bandwidth equals to the range of input time points. Try Gaussian kernel.")
    }
    if ((abs(opth - max(bw)) > .Machine$double.eps) && !is.infinite(optgcv)) 
      leave = TRUE
    else if (is.infinite(optgcv)) {
      if (verbose) 
        warning("Data is too sparse. Retry with larger bandwidths!")
      h0 = max(bw) * 1.01
    } else if ((abs(opth - max(bw)) > .Machine$double.eps)) {
      warning("The optimal bandwidth is not found in the candidate bandwidths. Retry with larger bandwidths")
      h0 = max(bw)
    }
    if (!leave) {
      newr = seq(0.5, 1, by = 0.05) * r
      ind = which(newr > h0)[1]
      q = (newr[ind]/h0)^(1/9)
      bw = q^(0:9) * h0
      if (verbose) {
        message("New bandwidth candidates:\n")
        print(bw)
      }
      iter = iter + 1
    }
  }
  
  ret = list(h = opth, gcv = optgcv, minBW = minBW)
  return(ret)
}

GetMean = function (x, t, outGrid, userBwMu = NULL, method = 'PACE') { 
  if (method == 'PACE') {
	win = rep(1, times = sum(lengths(t)))
  } 
  if (method == 'LH') {
	win = rep(1/lengths(t), times = lengths(t))
  }
  tin = unlist(t)
  xin = unlist(x)[order(tin)]
  win = win[order(tin)]
  tin = sort(tin)
  
  if (is.numeric(userBwMu)){
    if (userBwMu > 0) {
      bwMu = userBwMu
    }
  } else {
    bwMu = as.numeric(unlist(GCVLwls1D(yy = xin, tt = tin, win = win, kernel = 'epan', npoly = 1, nder = 0, dataType = 'Sparse'))[1])
    if (length(bwMu) == 0) {
      stop("The data is too sparse to estimate a mean function. Get more data!\n")
    }
  }

  mu = fdapace:::Lwls1D(bwMu, kernel_type = 'epan', npoly = 1, nder = 0, xin = tin, yin = xin, xout = outGrid, win = win)
  result = list(mu = mu, bwMu = bwMu, outGrid = outGrid)
  return(result)
}

GetPreAutoCov = function (x, t, dataType, error) {
  weight = rep(1/(lengths(t) * (lengths(t)-1)), times = lengths(t) * (lengths(t)-1))
  Ys = lapply(X = x, FUN = pracma::meshgrid)
  Xs = lapply(X = t, FUN = pracma::meshgrid)
  xx1 = unlist(do.call(rbind, lapply(Xs, "[", "X")))
  xx2 = unlist(do.call(rbind, lapply(Xs, "[", "Y")))
  yy1 = unlist(do.call(rbind, lapply(Ys, "[", "X")))
  yy2 = unlist(do.call(rbind, lapply(Ys, "[", "Y")))
  cyy = yy1 * yy2
  
  indx = rep(1:length(x), times = unlist(lapply(x, length))^2)
  tPairs = matrix(c(xx1, xx2), nrow = length(xx1), ncol = 2)
  if (error) {
    tneq = which(xx1 != xx2)
    teq = which(xx1 == xx2)
    indx = indx[tneq]
    diag = matrix(c(tPairs[teq, 1], cyy[teq]), ncol = 2)
    tPairs = tPairs[tneq, ]
    cxxn = cyy[tneq]
  } else {
    cxxn = cyy
  }
  result = list(tPairs = tPairs, cxxn = cxxn, weight = weight, indx = indx, cyy = cyy, 
                diag = diag, error = error, dataType = dataType)
  return(result)
}

GetAutoCov = function (x, t, obsGrid, regGrid, mu, userBwAutoCov = NULL) { 
	preAutoCov = GetPreAutoCov(x = x, t = t, dataType = optns$dataType, error = optns$error)

	if (is.numeric(userBwAutoCov)) {
  	if (userBwAutoCov > 0)
  		bwAutoCov = userBwAutoCov
	} else {
	  gcvObj = GCVLwls2D(t = t, obsGrid = obsGrid, regGrid = regGrid, kern = optns$kernel, 
		    rcov = preAutoCov, verbose = optns$verbose, win = preAutoCov$weight)
	  bwAutoCov = gcvObj$h
	}
	nonCentAutoCov = fdapace:::Lwls2D(bwAutoCov, kern = 'epan', xin = preAutoCov$tPairs, yin = preAutoCov$cxxn, 
						  win = preAutoCov$weight, xout1 = regGrid, xout2 = regGrid)
	autoCov = (nonCentAutoCov + t(nonCentAutoCov))/2 - tcrossprod(as.matrix(mu))

	preSigma2 = GetMean(x = lapply(x, function(x) x^2), t, outGrid = regGrid, userBwMu = NULL)
	bwSigma2 = preSigma2$bwMu
	sigma2 = integral(f = preSigma2$mu - mu^2 - diag(autoCov), g = NULL, domain = regGrid, type = 100) / diff(range(regGrid))
	if (sigma2 < 0) {
	  sigma2 = .Machine$double.eps
	}
    
	result = list(autoCov = autoCov, bwAutoCov = bwAutoCov, sigma2 = sigma2, bwSigma2 = bwSigma2, outGrid = regGrid)
	return(result)
}

GetCrossCov = function(y, x, t, outGrid, mu, userBwCrossCov = NULL, method = 'PACE') {
  xy = list()
  for (i in 1:length(y))
    xy[[i]] = x[[i]] * y[i]
  if (!is.null(userBwCrossCov)) {
    bwCrossCov = userBwCrossCov
    nonCentCrossCovObj = GetMean(xy, t, outGrid, userBwMu = bwCrossCov, method = method)
  } else {
    nonCentCrossCovObj = GetMean(xy, t, outGrid, userBwMu = NULL, method = method)
    bwCrossCov = nonCentCrossCovObj$bwMu
  }
  crossCov = nonCentCrossCovObj$mu - mu * mean(y)
  
  result = list(crossCov = crossCov, bwCrossCov = bwCrossCov, outGrid = outGrid)
  return(result)
}

PLEASS = function (yOld, xOld, tOld, optns, pMax = NULL, FVEthreshold, approxObj, yNew = NULL, xNew = NULL, tNew = NULL, yTrueNew = NULL, betaTrue = NULL, tuning = NULL) {
  if (!("fdapace" %in% rownames(installed.packages()))) 
    install.packages("fdapace")
  
  if (is.null(pMax)){
    if (!is.null(FVEthreshold))
      pMax = which.max(approxObj$cumFVE >= FVEthreshold) 
    else
      pMax = length(approxObj$eigenval)
  }
  
  ### local linear smoother
  regGrid = approxObj$regGrid
  mu = approxObj$mu
  autoCov = approxObj$autoCov
  sigma2 = approxObj$sigma2
  crossCov = approxObj$crossCov
  
  ### estimating beta
  # basis functions of Krylov subspace  
  KSBasis = matrix(NA, nrow = length(regGrid), ncol = pMax+1)
  KSBasis[, 1] = crossCov
  if (ncol(KSBasis) >= 2) {
    for (j in 2:ncol(KSBasis)) {
      KSBasis[, j] = integral(f = autoCov, g = KSBasis[, j-1], domain = regGrid, type = 211)
    }    
  }
  # orthonormalize Krylov basis and obtain APLS basis
  APLSBasis = GSortho(basis.origi = KSBasis, ker = autoCov, domain.x = regGrid)
  # beta.hat
  betaCoefVec = integral(f = t(APLSBasis), g = crossCov, domain = regGrid, type = 211)
  betaHat = matrix(NA, nrow = length(regGrid), ncol = pMax)
  betaHat[, 1] = APLSBasis[, 1] * betaCoefVec[1]
  if (pMax >= 2) {
    for (p in 2:pMax) {
      betaHat[, p] = APLSBasis[, p] * betaCoefVec[p] + betaHat[, p-1]
    }
  }
  
  ### data cleaning: exclude time points of regGrid
  regGrid.low = range(regGrid)[1]
  regGrid.upp = range(regGrid)[2]
  excludeOld = NULL
  for (i in 1:length(tOld)) {
    locTmp = (tOld[[i]] >= regGrid.low & tOld[[i]] <= regGrid.upp)
    tOld[[i]] = tOld[[i]][locTmp]
    xOld[[i]] = xOld[[i]][locTmp]
    if (sum(locTmp) < 2)
      excludeOld = c(excludeOld, i)
  }
  if (!is.null(excludeOld)){
    xOld = xOld[-excludeOld]
    tOld = tOld[-excludeOld]
    yOld = yOld[-excludeOld]
  }
  if ((!is.null(tNew)) & (!is.null(xNew))){
    excludeNew = NULL
    for (i in 1:length(tNew)) {
      locTmp = (tNew[[i]] >= regGrid.low & tNew[[i]] <= regGrid.upp)
      tNew[[i]] = tNew[[i]][locTmp]
      xNew[[i]] = xNew[[i]][locTmp]
      if (sum(locTmp) < 2)
        excludeNew = c(excludeNew, i)
    }
    if (!is.null(excludeNew)){
      xNew = xNew[-excludeNew]
      tNew = tNew[-excludeNew]
      if (!is.null(yNew)){
        yNew = yNew[-excludeNew]
      }
    }
  }
  xAll = c(xOld, xNew)
  tAll = c(tOld, tNew)
  
  ### prediction
  # cov between X and Xi's
  covXandXi = matrix(NA, nrow = length(regGrid), ncol = pMax)
  for (j in 1:pMax) {
    covXandXi[, j] = integral(f = autoCov, g = APLSBasis[, j], domain = regGrid, type = 211)
  }
  # interpolation
  obsGrid = sort(unique(unlist(tAll)))
  muObs = fdapace::ConvertSupport(fromGrid = regGrid, toGrid = obsGrid, mu = mu)
  autoCovObs = fdapace::ConvertSupport(fromGrid = regGrid, toGrid = obsGrid, Cov = autoCov)
  covXandXiObs = fdapace::ConvertSupport(fromGrid = regGrid, toGrid = obsGrid, phi = covXandXi)
  # eta.hat
  etaHatAll = matrix(NA, nrow = length(tAll), ncol = pMax)
  sdEtaHatAll = matrix(NA, nrow = length(tAll), ncol = pMax)
  for (i in 1:nrow(etaHatAll)){
    locCurt = obsGrid %in% tAll[[i]]
    xCurt = xAll[[i]]
    muCurt = muObs[locCurt]
    sigmaHatCurtInv = solve(sigma2 * diag(length(xCurt)) + autoCovObs[locCurt, locCurt])
    for (p in 1:pMax) {
      xiVec =  t(covXandXiObs[locCurt, 1:p]) %*% sigmaHatCurtInv %*% (xCurt - muCurt)
      etaHatAll[i, p] = mean(yOld) + crossprod(betaCoefVec[1:p], xiVec)
      sdEtaHatAll[i, p] = abs(
        t(as.matrix(betaCoefVec[1:p])) %*%
        (diag(p) - t(covXandXiObs[locCurt, 1:p]) %*% sigmaHatCurtInv %*% covXandXiObs[locCurt, 1:p]) %*% 
        as.matrix(betaCoefVec[1:p])
      )^.5
    }
  }
  
  ### Error measures
  isee = NULL; reisee = NULL
  mspe = NULL; remspe = NULL
  coverPr = NULL
  coverMissL = NULL
  coverMissR = NULL
  if (!is.null(betaTrue)){
    isee = apply(betaHat, 2, function(x) integral(f=(x-betaTrue)^2, g = NULL, domain=regGrid, type=100))
    reisee = isee/integral(f=betaTrue^2, g = NULL, domain=regGrid, type=100)
  }
  if (!is.null(yNew)){
    etaHatNew = as.matrix(etaHatAll[-(1:length(yOld)), ])
    mspe = apply(etaHatNew, 2, function(x) mean((x - yNew)^2))
    remspe = apply(etaHatNew, 2, function(x) mean((x - yNew)^2/mean((mean(yOld) - yNew)^2)))
    if (!is.null(yTrueNew)){
      sdEtaHatNew = as.matrix(sdEtaHatAll[-(1:length(yOld)), ])
      coverPr = colMeans(as.matrix(
        apply(etaHatNew-2*sdEtaHatNew, 2, function(x) yTrueNew >= x) *
          apply(etaHatNew+2*sdEtaHatNew, 2, function(x) yTrueNew <= x)
      ))
      
      coverMissL= colMeans(as.matrix(
        apply(etaHatNew-2*sdEtaHatNew, 2, function(x) yTrueNew < x)
      ))
      coverMissR= colMeans(as.matrix(
        apply(etaHatNew+2*sdEtaHatNew, 2, function(x) yTrueNew > x)
      ))
    }
  }
  
  if (is.null(tuning))
    tuning = F
  if (tuning){
    etaHatOld = as.matrix(etaHatAll[1:length(yOld), ])
    CVmat = array(NA, dim = dim(etaHatOld))
    for (i in 1:length(yOld)){
      CVmat[i, ] = (etaHatOld[i, ] - mean(yOld) + mean(yOld[-i]) - yOld[i])^2
    }
    pOpt = which.min(colMeans(CVmat))
    if (!is.null(betaTrue)){
      isee = isee[pOpt]
      reisee = reisee[pOpt]
    }
    if (!is.null(yNew)){
      mspe = mspe[pOpt]
      remspe = remspe[pOpt]
      if (!is.null(yTrueNew)){
        coverPr = coverPr[pOpt]
        coverMissL = coverMissL[pOpt]
        coverMissR = coverMissR[pOpt]
      }
    }
  }

  return(list(isee = isee, reisee = reisee, mspe = mspe, remspe = remspe, coverPr = coverPr, coverMissL = coverMissL, coverMissR = coverMissR))
}

PACE = function (yOld, xOld, tOld, optns, pMax = NULL, FVEthreshold, approxObj, yNew = NULL, xNew = NULL, tNew = NULL, yTrueNew = NULL, betaTrue = NULL, tuning = NULL) {
  if (!("fdapace" %in% rownames(installed.packages()))) 
    install.packages("fdapace")
  
  if (is.null(pMax)){
    if (!is.null(FVEthreshold))
      pMax = which.max(approxObj$cumFVE >= FVEthreshold) 
    else
      pMax = length(approxObj$eigenval)
  }
  
  ### local linear smoother
  regGrid = approxObj$regGrid
  mu = approxObj$mu
  autoCov = approxObj$autoCov
  sigma2 = approxObj$sigma2
  crossCov = approxObj$crossCov
  
  ### estimating beta
  # eigen basis functions
  if (pMax == 1){
    eigenBasis = as.matrix(approxObj$eigenfun[, 1])
  }else {
    eigenBasis = approxObj$eigenfun[, 1:pMax]
  }
  # beta.hat
  betaCoefVec = as.vector(integral(f = t(eigenBasis), g = crossCov, domain = regGrid, type = 211))/approxObj$eigenval[1:pMax]
  betaHat = matrix(NA, nrow = length(regGrid), ncol = pMax)
  betaHat[, 1] = eigenBasis[, 1] * betaCoefVec[1]
  if (pMax >= 2) {
    for (p in 2:pMax) {
      betaHat[, p] = eigenBasis[, p] * betaCoefVec[p] + betaHat[, p-1]
    }
  }

  ### data cleaning: exclude time points of regGrid
  regGrid.low = range(regGrid)[1]
  regGrid.upp = range(regGrid)[2]
  excludeOld = NULL
  for (i in 1:length(tOld)) {
    locTmp = (tOld[[i]] >= regGrid.low & tOld[[i]] <= regGrid.upp)
    tOld[[i]] = tOld[[i]][locTmp]
    xOld[[i]] = xOld[[i]][locTmp]
    if (sum(locTmp) < 2)
      excludeOld = c(excludeOld, i)
  }
  if (!is.null(excludeOld)){
    xOld = xOld[-excludeOld]
    tOld = tOld[-excludeOld]
    yOld = yOld[-excludeOld]
  }
  if ((!is.null(tNew)) & (!is.null(xNew))){
    excludeNew = NULL
    for (i in 1:length(tNew)) {
      locTmp = (tNew[[i]] >= regGrid.low & tNew[[i]] <= regGrid.upp)
      tNew[[i]] = tNew[[i]][locTmp]
      xNew[[i]] = xNew[[i]][locTmp]
      if (sum(locTmp) < 2)
        excludeNew = c(excludeNew, i)
    }
    if (!is.null(excludeNew)){
      xNew = xNew[-excludeNew]
      tNew = tNew[-excludeNew]
      if (!is.null(yNew)){
        yNew = yNew[-excludeNew]
      }
    }
  }
  xAll = c(xOld, xNew)
  tAll = c(tOld, tNew)
  
  ### prediction
  # interpolation
  obsGrid = sort(unique(unlist(tAll)))
  muObs = fdapace::ConvertSupport(fromGrid = regGrid, toGrid = obsGrid, mu = mu)
  autoCovObs = fdapace::ConvertSupport(fromGrid = regGrid, toGrid = obsGrid, Cov = autoCov)
  eigenBasisObs = fdapace::ConvertSupport(fromGrid = regGrid, toGrid = obsGrid, phi = eigenBasis)
  covXandXiObs = eigenBasisObs *
    matrix(approxObj$eigenval[1:pMax], nrow = length(obsGrid), ncol = pMax, byrow = T)
  # eta.hat
  etaHatAll = matrix(NA, nrow = length(tAll), ncol = pMax)
  sdEtaHatAll = matrix(NA, nrow = length(tAll), ncol = pMax)
  x.hat = list()
  for (i in 1:nrow(etaHatAll)){
    locCurt = obsGrid %in% tAll[[i]]
    xCurt = xAll[[i]]
    muCurt = muObs[locCurt]
    sigmaHatCurtInv = solve(sigma2 * diag(length(xCurt)) + autoCovObs[locCurt, locCurt])
    x.hat[[i]] = matrix(NA, nrow = length(locCurt), ncol = pMax)
    for (p in 1:pMax) {
      xiVec =  t(covXandXiObs[locCurt, 1:p]) %*% sigmaHatCurtInv %*% (xCurt - muCurt)
      etaHatAll[i, p] = mean(yOld) + crossprod(betaCoefVec[1:p], xiVec)
      if (p == 1){
        LambdaMat = as.matrix(approxObj$eigenval[1])
      }else{
        LambdaMat = diag(approxObj$eigenval[1:p])
      }
      sdEtaHatAll[i, p] = abs(
        t(as.matrix(betaCoefVec[1:p])) %*%
          (LambdaMat - t(covXandXiObs[locCurt, 1:p]) %*% sigmaHatCurtInv %*% covXandXiObs[locCurt, 1:p]) %*%
          as.matrix(betaCoefVec[1:p])
      )^.5
      
      x.hat[[i]][, p] = muCurt + eigenBasisObs[, 1:p] %*% xiVec
    }
  }
    
  ### Error measures
  isee = NULL; reisee = NULL
  mspe = NULL; remspe = NULL
  coverPr = NULL
  coverMissL = NULL
  coverMissR = NULL 
  if (!is.null(betaTrue)){
    betaTrue.norm.sq = integral(f = betaTrue^2, g = NULL, domain = regGrid, type = 100)
    isee = apply(betaHat, 2, function(x) integral(f=(x-betaTrue)^2, g = NULL, domain=regGrid, type=100))
    reisee = isee/betaTrue.norm.sq
  }

  if (!is.null(yNew)){
    etaHatNew = as.matrix(etaHatAll[-(1:length(yOld)), ])
    mspe = apply(etaHatNew, 2, function(x) mean((x - yNew)^2))
    remspe = apply(etaHatNew, 2, function(x) mean((x - yNew)^2/mean((mean(yOld) - yNew)^2)))
    if (!is.null(yTrueNew)){
      sdEtaHatNew = as.matrix(sdEtaHatAll[-(1:length(yOld)), ])
      coverPr = colMeans(as.matrix(
        apply(etaHatNew-2*sdEtaHatNew, 2, function(x) yTrueNew >= x) *
          apply(etaHatNew+2*sdEtaHatNew, 2, function(x) yTrueNew <= x)
      ))
      
      coverMissL= colMeans(as.matrix(
        apply(etaHatNew-2*sdEtaHatNew, 2, function(x) yTrueNew < x)
      ))
      coverMissR= colMeans(as.matrix(
        apply(etaHatNew+2*sdEtaHatNew, 2, function(x) yTrueNew > x)
      ))
      
    }
  }
  
  if (is.null(tuning))
    tuning = F
  if (tuning){
    AICmat = matrix(NA, nrow = length(xOld), ncol = pMax)
    for (i in 1:length(xOld)){
      for (p in 1:pMax){
        AICmat[i, p] = sum((xOld[[i]] - x.hat[[i]][, p])^2)/(2 * sigma2) +
          lengths(xOld)[i]/2 * (log(2*pi) + log(sigma2))
      }
    }
    pOpt = which.min(colSums(AICmat) + (1:pMax))
    
    if (!is.null(betaTrue)){
      isee = isee[pOpt]
      reisee = reisee[pOpt]
    }
    if (!is.null(yNew)){
      mspe = mspe[pOpt]
      remspe = remspe[pOpt]
      if (!is.null(yTrueNew)){
        coverPr = coverPr[pOpt]
        coverMissL = coverMissL[pOpt]
        coverMissR = coverMissR[pOpt]
      }
    }
  }else{
    if (!is.null(betaTrue)){
      isee = isee[pMax]
      reisee = reisee[pMax]
    }
    if (!is.null(yNew)){
      mspe = mspe[pMax]
      remspe = remspe[pMax]
      if (!is.null(yTrueNew)){
        coverPr = coverPr[pMax]
        coverMissL = coverMissL[pMax]
        coverMissR = coverMissR[pMax]
      }
    }
  }
  
  return(list(isee = isee, reisee = reisee, mspe = mspe, remspe = remspe, coverPr = coverPr, coverMissL = coverMissL, coverMissR = coverMissR))
}

#' integrate error measures

creatErrMat = function(colNames, errLst){
  for (i in 1:length(errLst)){
    if (i == 1)
      nrowErrMat = nrow(errLst[[1]])
    if (i > 1){
      if (nrowErrMat > nrow(errLst[[i]]))
        nrowErrMat = nrow(errLst[[i]])
    }
  }
  for (i in 1:length(errLst)){
    if (i == 1)
      errMat = errLst[[1]][1:nrowErrMat, ]
    if (i > 1){
      errMat = cbind(
        errMat, errLst[[i]][1:nrowErrMat, ]
      )
    }
  }
  colnames(errMat) = colNames
  return(errMat)
}

#'  boxplot of error measures

boxplotErr = function(nrowplot, errMat, type){
  if (!("ggplot2" %in% rownames(installed.packages()))) 
    install.packages("ggplot2")
  if (!("reshape2" %in% rownames(installed.packages()))) 
    install.packages("reshape2")
  library(ggplot2)
  library(reshape2)
  
  melton = melt(
    data.frame(errMat[1:nrowplot,], 
               Replication = 1:nrowplot), 
    id.vars = "Replication")
  bplot = ggplot(melton, aes(x = variable, y = value, colour = variable)) + 
    geom_boxplot(outlier.shape = NA, width = .5) +
    coord_cartesian(ylim = c(0, 
                             ifelse(
                               type == "Coverage Percentage", 
                               1, 
                               ifelse(type == "ReMSPE" | type == "ReISEE", 2, 10)
                               )
                             )
                    ) +
    labs(x = '', y = type, title = '') + 
    theme_bw() +
    theme(legend.position = "none", 
          panel.border = element_blank(), 
          plot.title = element_text(hjust = 0.5),
          text = element_text(size=20),
          axis.text.x = element_text(angle=0, hjust=0.5)
    )
  plot(bplot)
  
  type = ifelse(type == "Coverage Percentage", "Cover", type)
  file = ifelse(simu,
                paste0(type, '_',
                       'Simu', simuCase, '_',
                       nSubject, 'nSubject_',
                       J, 'eigen_',
                       FVEthreshold * 1e4, 'FVE_',
                       RR, 'repeats_',
                       optns$nRegGrid, 'nRegGrid_',
                       SNRy, 'SNR_',
                       propTrain * 100, 'train.pdf'),
                paste0(type, '_',
                       'Real', realCase, '_', 
                       FVEthreshold * 1e4, 'FVE_',
                       RR, 'repeats_', 
                       optns$nRegGrid, 'nRegGrid_',
                       propTrain * 100, 'train.pdf')
  )
  ggsave(file = file, 
         width = 6,
         height = 8,
         units = 'in',
         dpi = 300,
         path = "figure")
}

FACE = function (y = NULL, x, t, optns) {
  if (!("face" %in% rownames(installed.packages()))) 
    install.packages("face")
  # reformulate data
  argvals = numeric(0)
  subj = numeric(0)
  x.val = numeric(0)
  xy.val = numeric(0)
  for (i in 1:length(x)){
    for (j in 1:lengths(x)[i]){
      argvals = c(argvals, t[[i]][j])
      subj = c(subj, i)
      x.val = c(x.val, x[[i]][j])
      xy.val = c(xy.val, x[[i]][j] * y[i])
    }
  }
  data.x = data.frame(argvals = argvals, subj = subj, y = x.val)
  data.xy = data.frame(argvals = argvals, subj = subj, y = xy.val)
  regGrid = seq(from = min(unlist(t)), to = max(unlist(t)), length.out = optns$nRegGrid)
  
  # FACE
  faceObj.x = face::face.sparse(data = data.x, argvals.new = regGrid, pve = optns$FVEthreshold)
  faceObj.xy = face::face.sparse(data = data.xy, argvals.new = regGrid, pve = optns$FVEthreshold)
  
  # make sure autoCov is non-neg definite
  autoCov = (faceObj.x$Chat.new + t(faceObj.x$Chat.new))/2
  eig <- eigen(autoCov)
  positiveInd <- eig[["values"]] >= 0
  d <- eig[["values"]][positiveInd]
  eigenV <- eig[["vectors"]][, positiveInd]
  if (optns$maxK < length(d)) {
    d <- d[1:optns$maxK]
    eigenV <- eigenV[, 1:optns$maxK]
  }
  phi <- apply(eigenV, 2, function(x) {
    x <- x/sqrt(integral(f = x^2, g = NULL, domain = regGrid, type = 100))
  })
  lambda <- mean(diff(regGrid)) * d
  cumFVE = cumsum(lambda)/sum(lambda) * 100
  autoCov <- phi %*% diag(x = lambda, nrow = length(lambda)) %*% t(phi)
  
  # other outputs
  mu = faceObj.x$mu.new
  sigma2 = mean(faceObj.x$var.error.new)
  eigenfun = phi[, 1:min(length(d), optns$maxK)]
  eigenval = lambda[1:min(length(d), optns$maxK)]
  if (is.null(y)) {
    crossCov = NULL
  } else {
    crossCov = faceObj.xy$mu.new - mu * mean(y)
  }
  
  return(list(
      obsGrid = sort(unlist(argvals)),
      regGrid = regGrid,
      mu = mu,
      autoCov = autoCov,
      crossCov = crossCov,
      sigma2 = sigma2,
      eigenfun = eigenfun,
      eigenval = eigenval,
      cumFVE = cumFVE
    )
  )
}

llsFit = function (y = NULL, x, t, optns, method = 'PACE') {
  if (!("fdapace" %in% rownames(installed.packages()))) install.packages("fdapace")
  
  if (method == 'PACE') {
    llsObjRaw = fdapace::FPCA(x, t, optns)
    obsGrid = llsObjRaw$obsGrid
    regGrid = llsObjRaw$workGrid
    mu = llsObjRaw$mu
    autoCov = llsObjRaw$fittedCov
    sigma2 = llsObjRaw$sigma2
    eigenfun = llsObjRaw$phi
    eigenval = llsObjRaw$lambda
    cumFVE = llsObjRaw$cumFVE
  }
  if (method == 'LH') {
    obsGrid = sort(unique(unlist(t)))
    regGrid = seq(from = min(obsGrid), to = max(obsGrid), length.out = optns$nRegGrid)
    mu = GetMean(x, t, regGrid, method = 'LH')$mu
    autoCovObj = GetAutoCov(x, t, obsGrid, regGrid, mu = mu)
    eigObj = fdapace:::GetEigenAnalysisResults(smoothCov = autoCovObj$autoCov, regGrid, optns, mu = mu)
    autoCov = eigObj$fittedCov
    sigma2 = autoCovObj$sigma2
    eigenfun = eigObj$phi
    eigenval = eigObj$lambda
    cumFVE = eigObj$cumFVE
  }
  
  if (is.null(y)) {
    crossCov = NULL
  } else {
    crossCov = GetCrossCov(y, x, t, regGrid, mu, method = method)$crossCov
  }
  
  return(list(obsGrid = obsGrid, regGrid = regGrid, 
    mu = mu, autoCov = autoCov, crossCov = crossCov, 
    sigma2 = sigma2, eigenfun = eigenfun, eigenval = eigenval, cumFVE = cumFVE))
}

qlFit = function(y, xList, xyList, tList, init.x, init.xy, optns, optnsQl){
  if (!("fdapace" %in% rownames(installed.packages()))) install.packages("fdapace")
  
  # initiate
  y = yOld
  xList = xOld
  xyList = xyOld
  tList = tOld
  init.x = llsObj.x
  init.xy = llsObj.xy
  
  obsGrid = init.x$obsGrid
  regGrid = init.x$regGrid
  Kfold = optnsQl$Kfold
  hLength = optnsQl$hLength
  h1All = runif(n = hLength, min = 0, max = 2)
  h2All = runif(n = hLength, min = 0, max = 2)
  h3All = runif(n = hLength, min = 0, max = 2)
  h4All = runif(n = hLength, min = 0, max = 2)
  del = optnsQl$del
  
  mu0.x = init.x$mu
  mu0.xy = init.xy$mu
  autoCov0.x = init.x$autoCov
  autoCov0.xy = init.xy$autoCov
  sigma20.x = init.x$sigma2
  sigma20.xy = init.xy$sigma2

  sigmaFun0 = diag(autoCov0.x)^.5
  tRhoFun0 = (regGrid[2] - regGrid[1]) * (0:(length(regGrid)-1))
  rhoMat0 = diag(1/sigmaFun0) %*% autoCov0.x %*% diag(1/sigmaFun0)
  rhoFun0 = c(1, numeric(length(regGrid)-1))
  for (t in 1:(length(tRhoFun0)-1)){
    if (t == length(tRhoFun0)-1)
      rhoFun0[t+1] = (rhoMat0[t+1, 1] + rhoMat0[1, t+1])/2
    else
      rhoFun0[t+1] = mean(diag(rhoMat0[-tail(1:length(tRhoFun0), t), -(1:t)]) + diag(rhoMat0[-(1:t), -tail(1:length(tRhoFun0), t)]))
  }

  # transform the datatype frome 'list' to 'matrix'
  xMat = matrix(NA, nrow = length(xList), ncol = max(lengths(xList)))
  xyMat = matrix(NA, nrow = length(xyList), ncol = max(lengths(xyList)))
  tMat = matrix(NA, nrow = length(tList), ncol = max(lengths(tList)))
  for (i in 1:nrow(xMat)){
    xMat[i, 1:lengths(xList)[i]] = xList[[i]]
    xyMat[i, 1:lengths(xyList)[i]] = xyList[[i]]
    tMat[i, 1:lengths(tList)[i]] = tList[[i]]
  }
  
  # Kfold CV
  indexNew = sample.int(nrow(xMat), size = nrow(xMat), replace = F)
  CVkji.h1 = array(0, dim = c(Kfold, hLength, length(xList)))
  CVkji.h4 = array(0, dim = c(Kfold, hLength, length(xList)))
  CVkji.h23 = array(0, dim = c(Kfold, hLength, length(xList)))
  
  for (k in 1:Kfold){
    if (k < Kfold)
      vali = ((k-1) * round(length(indexNew)/Kfold) + 1):(k * round(length(indexNew)/Kfold))
    else 
      vali = ((k-1) * round(length(indexNew)/Kfold) + 1):length(indexNew)
    for (j in 1:hLength){
      mu = est.utf(ut0 = mu0.x, h1 = h1All[j], cov0 = autoCov0.x, lam20 = sigma20.x, tt = tMat[indexNew[-vali], ], Y = xMat[indexNew[-vali], ], del = del, tts = regGrid, stcov = regGrid)
      
      crossCovRaw = est.utf(ut0 = mu0.xy, h1 = h4All[j], cov0 = autoCov0.xy, lam20 = sigma20.xy, tt = tMat[indexNew[-vali], ], Y = xyMat[indexNew[-vali], ], del = del, tts = regGrid, stcov = regGrid)
      
      muObs = fdapace::ConvertSupport(fromGrid = regGrid, toGrid = obsGrid, mu = mu)
      crossCovRawObs = fdapace::ConvertSupport(fromGrid = regGrid, toGrid = obsGrid, mu = crossCovRaw)
      
      for (i in 1:length(vali)) {
        CVkji.h1[k, j, i] = mean((xList[[indexNew[vali[i]]]] - muObs[obsGrid %in% tList[[indexNew[vali[i]]]]])^2)
        CVkji.h4[k, j, i] = mean((xyList[[indexNew[vali[i]]]] - crossCovRawObs[obsGrid %in% tList[[indexNew[vali[i]]]]])^2)
        cat(c(k, j, i), '\n')
      }
    }
  }
  
  h1Opt = h1All[which.min(colSums(apply(CVkji.h1, 1:2, mean)))]
  h4Opt = h4All[which.min(colSums(apply(CVkji.h4, 1:2, mean)))]
  mu = est.utf(ut0 = mu0.x, h1 = h1Opt, cov0 = autoCov0.x, lam20 = sigma20.x, tt = tMat, Y = xMat, del = del, tts = regGrid, stcov = regGrid)
  crossCov = - mean(y) * mu + est.utf(ut0 = mu0.xy, h1 = h4Opt, cov0 = autoCov0.xy, lam20 = sigma20.xy, tt = tMat, Y = xyMat, del = del, tts = regGrid, stcov = regGrid)
  
  for (k in 1:Kfold){
    if (k < Kfold)
      vali = ((k-1) * round(length(indexNew)/Kfold) + 1):(k * round(length(indexNew)/Kfold))
    else 
      vali = ((k-1) * round(length(indexNew)/Kfold) + 1):length(indexNew)
    for (j in 1:hLength){
      covObj = estva(ut0 = mu, h2 = h2All[j], h3 = h3All[j], SIGT0 = sigmaFun0, RHOT0 = rhoFun0, lam20 = sigma20.x, tt = tMat[indexNew[-vali], ], Y = xMat[indexNew[-vali], ], del = del, tgrids = regGrid, tgridrho = tRhoFun0, tts = regGrid)
      rhoFun = covObj$rho[, 1]
      sigmaFun = covObj$sigma[, 1]
      sigma2.x = covObj$lam
      
      sigmaFunObs = fdapace::ConvertSupport(fromGrid = regGrid, toGrid = obsGrid, mu = sigmaFun)
      
      for (i in 1:length(vali)) {
        xCurrent = xList[[indexNew[vali[i]]]]
        tCurrent = tList[[indexNew[vali[i]]]]
        muCurrent = fdapace::ConvertSupport(fromGrid = regGrid, toGrid = tCurrent, mu = mu)
        sigmaFunCurrent = fdapace::ConvertSupport(fromGrid = regGrid, toGrid = tCurrent, mu = sigmaFun)
        for (l in 1:length(tCurrent)){
          for (m in 1:length(tCurrent)){
            CVkji.h23[k, j, i] = CVkji.h23[k, j, i] + (
                (xCurrent[l] - muCurrent[l])(xCurrent[m] - muCurrent[m]) - 
                sigmaFunCurrent[l] * sigmaFunCurrent[m] * fdapace::ConvertSupport(fromGrid = regGrid, toGrid = c(abs(l-m)), mu = sigmaFun) -
                sigma2.x * ifelse(l == m, 1, 0)
              )^2
          }
        }
        CVkji.h23[k, j, i] = CVkji.h23[k, j, i]/length(tCurrent)/(length(tCurrent)-1)
      }
    }
  }
  
  h2Opt = h2All[which.min(colSums(apply(CVkji.h23, 1:2, mean)))]
  h3Opt = h3All[which.min(colSums(apply(CVkji.h23, 1:2, mean)))]
  covObj = estva(ut0 = mu, h2 = h2Opt, h3 = h3Opt, SIGT0 = sigmaFun0, RHOT0 = rhoFun0, lam20 = sigma20.x, tt = tMat, Y = xMat, del = del, tgrids = regGrid, tgridrho = tRhoFun0, tts = regGrid)
  rhoFun = covObj$rho[, 1]
  sigmaFun = covObj$sigma[, 1]
  sigma2.x = covObj$lam
  
  rhoMat = array(NA, dim = dim(rhoMat0))
  for (l in 1:nrow(rhoMat)) {
    for (m in 1:ncol(rhoMat)) {
      rhoMat[l, m] = rhoFun[abs(l-m) + 1]
    }
  }
  autoCov = diag(sigmaFun) %*% rhoMat %*% diag(sigmaFun)
  eigenObj = fdapace:::GetEigenAnalysisResults(smoothCov = autoCov, regGrid = regGrid, optns = optns, muWork = mu) 
  eigenfun = eigenObj$phi
  eigenval = eigenObj$lambda
  
  return(list(mu = mu, crossCov = crossCov, autoCov = autoCov, sigma2 = sigma2.x, regGrid = regGrid, obsGrid = obsGrid, eigenfun = eigenfun, eigenval = eigenval))
}

orthoBasis = function(order, denseGrid, type, normalized = T){
  if (!("orthopolynom" %in% rownames(installed.packages()))) 
    install.packages("orthopolynom")
  library(orthopolynom)
  
  if (type == 'shiftedLegendre') {
    poly.value = polynomial.values(slegendre.polynomials(n = max(order), normalized = normalized), x = denseGrid)
    res = matrix(NA, nrow = length(order), ncol = length(denseGrid))
    for (i in 1:length(order)){
      res[i, ] = poly.value[[order[i]+1]]
    }
  }
  return(res)
}

test.ortho = function(basis.ortho, ker, domain.x, domain.y = NULL){
  if (is.null(domain.y) & is.null(ker)){
    p = ncol(basis.ortho)
    H = array(NA, dim = c(p, p))
    for (i in 1:p){
      for (j in 1:p){
        H[i, j] = integral(basis.ortho[, i] * basis.ortho[, j], domain = domain.x, type = 100)
      }
    }
  }
  
  if (is.null(domain.y) & !is.null(ker)){
    p = ncol(basis.ortho)
    H = array(NA, dim = c(p, p))
    for (i in 1:p){
      for (j in 1:p){
        integ = integral(ker, basis.ortho[, j], domain.x, type = 211)
        H[i, j] = integral(basis.ortho[, i] * integ, domain = domain.x, type = 100)
      }
    }
  }
  
  if (!is.null(domain.y) & is.null(ker)){
    p = dim(basis.ortho)[3]
    H = array(NA, dim = c(p, p))
    for (i in 1:p){
      for (j in 1:p){
        integ = integral(basis.ortho[, , i] * basis.ortho[, , j], domain = domain.y, type = 201)
        H[i, j] = integral(integ, domain = domain.x, type = 100)
      }
    }
  }
  
  if (!is.null(domain.y) & !is.null(ker)){
    p = dim(basis.ortho)[3]
    H = array(NA, dim = c(p, p))
    for (i in 1:p){
      for (j in 1:p){
        integ1 = integral(ker, basis.ortho[, , j], domain.x, type = 222)
        integ2 = integral(basis.ortho[, , i] * integ1, domain = domain.y, type = 201)
        H[i, j] = integral(integ2, domain = domain.x, type = 100)
      }
    }
  }
  
  return(H)
}

