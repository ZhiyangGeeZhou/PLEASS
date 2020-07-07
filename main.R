#' Put main.R & supportFun.R into the identical folder.
#' In the same folder, create two subfolders named figure' & 'Rimage', respectively.
rm(list = ls())
if (!("rstudioapi" %in% rownames(installed.packages()))) 
  install.packages("rstudioapi")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('supportFun.R')

#### Global parameters
set.seed(1)
options(warn = -1, digits = 4)
RR = 200 # number of replication
propTrain = .8 # proportion of training set
simu = T # T for simulation, F for real data
tuning = T
pMax = NULL
FVEthreshold = .95
llsMethod = 'PACE' # 'PACE' for Yao, Muller & Wang (2005); 'LH' for Li & Hsing (2010)
quasilk = F # quasi-likelihood estimation of cov function following Zhou, Lin & Liang (2018)
## parameters for simulation
simuCase = 3
J = 9 # number of basis functions
SNRx = 100 # SNR for errorX in simulation
SNRy = SNRx # SNR for errorY in simulation
## parameters for real case
sudoAnalyze = F # F: using testing+training (instead of training only) data to fit cov structure
realCase = 1 # 1 for PBC; 2 for DTI

#### Simulated data

if (simu) {
  nSubject = 1000
  tRange = c(0, 1)
  denseGrid = seq(from = tRange[1], to = tRange[2], length.out = 101)
  Lambda = c(seq(from = 100, by = -10, length.out = ceiling(J/3)), 
             seq(from = 10, by = -1, length.out = ceiling(J/3)),
             seq(from = 1, by = -.1, length.out = J-2*ceiling(J/3)))

  # xMuFun = function(t) t + exp(-t^2)
  yMuMat = matrix(0, nrow = nSubject, ncol = 1)
  xMuMat = matrix(0, nrow = nSubject, ncol = length(denseGrid))
  eigenFun = orthoBasis(order = 1:J, denseGrid, type = 'shiftedLegendre', normalized = T) 
  
  betaTrueCoef = numeric(J)
  if (simuCase == 1){
    betaTrueCoef = c(rep(1, times = ceiling(J/3)), rep(0, times = J - ceiling(J/3)))
  }
  if (simuCase == 2){
    betaTrueCoef = c(rep(0, times = floor((J - ceiling(J/3))/2)), rep(1, times = ceiling(J/3)), rep(0, times = ceiling((J - ceiling(J/3))/2)))
  }
  if (simuCase == 3){
    betaTrueCoef = c(rep(0, times = J - ceiling(J/3)), rep(1, times = ceiling(J/3)))
  }

  betaTrue = t(betaTrueCoef) %*% eigenFun
  
  tAll = list()
  xAll = list()
  yAll = list()
  yTrueAll = list()
  xyAll = list()
  
  for (R in 1:RR) {
    eigenScore = matrix(NA, nrow = nSubject, ncol = J)
    for (j in 1:J){
      eigenScore[, j] = Lambda[j]^.5 * rnorm(nSubject)
    }
    xTrueTmp = eigenScore %*% eigenFun + xMuMat
    yTrueTmp = integral(f = xTrueTmp - xMuMat, g = betaTrue, domain = denseGrid, type = 211) + yMuMat
    errorX = (sum(Lambda)/SNRx)^.5 * matrix(rnorm(nSubject * length(denseGrid)), nrow = nSubject, ncol = length(denseGrid))
    errorY = (var(yTrueTmp)/SNRy)^.5 * rnorm(nSubject)
    xContami = xTrueTmp + errorX
    yTmp = yTrueTmp + errorY
    
    tTmp = list()
    xTmp = list()
    xyTmp = list()
    for (i in 1:nSubject) {
      indexTmp = sort(sample(1:length(denseGrid), size = sample(3:6, size = 1)))
      tTmp[[i]] = denseGrid[indexTmp]
      xTmp[[i]] = xContami[i, indexTmp]
      xyTmp[[i]] = yTmp[i] * xTmp[[i]]
    }
    
    tAll[[R]] = tTmp
    xAll[[R]] = xTmp
    yTrueAll[[R]] = yTrueTmp
    yAll[[R]] = yTmp
    xyAll[[R]] = xyTmp
  }
}

#### Real data

if (!simu) {
  
  if (realCase == 1) {
    # PBC data (pbcseq in R package survival)
    
    if (!("survival" %in% rownames(installed.packages()))) {
      install.packages("survival")
    }
    
    dayMax = 3000
    tmp = survival::pbcseq[
        survival::pbcseq$day <= dayMax 
        &
        !is.na(survival::pbcseq$alk.phos)
      , ]
    cor(tmp[, c('bili', 'albumin', 'alk.phos', 'ast', 'platelet', 'protime')])
    
    attach(tmp)
    
    y = numeric()
    x = list()
    xy = list()
    t = list()
    idLeft = sort(unique(id))
    for (i in idLeft) {
      y[i] = tail(ast[id == i], n=1)
      x[[i]] = alk.phos[id == i]
      xy[[i]] = y[i] * x[[i]]
      t[[i]] = day[id == i]/dayMax
    }
    idLeft = idLeft[lengths(t) >= 2]
    y = (y[lengths(t) >= 2])
    x = x[lengths(t) >= 2]
    xy = xy[lengths(t) >= 2]
    t = t[lengths(t) >= 2]
    
    fdapace::CreateDesignPlot(t, isColorPlot = F, noDiagonal = F, addLegend = F)
    detach(tmp)
  } 
 
  if (realCase == 2) {
    if (!("refund" %in% rownames(installed.packages()))) install.packages("refund")
    
    tmp = refund::DTI
    
    attach(tmp)
    
    pasatFilter = pasat[!is.na(pasat)]
    idFilter = ID[!is.na(pasat)]
    visitFilter = visit[!is.na(pasat)]
    ccaFilter = cca[!is.na(pasat), ]
    rcstFilter = rcst[!is.na(pasat), ]
    
    y = numeric()
    x = list()
    xy = list()
    t = list()
    
    ids = unique(idFilter)
    for (i in 1:length(ids)) {
      if (sum(idFilter == ids[i]) >= 2) {
        visitCurt = visitFilter[idFilter == ids[i]]
        (visitMax = max(visitCurt))
        (visitMin = min(visitCurt))
        
        y[i] = pasatFilter[idFilter == ids[i] & visitFilter == visitMax]
        x[[i]] = ccaFilter[idFilter == ids[i] & visitFilter == visitMax, ][!is.na(ccaFilter[idFilter == ids[i] & visitFilter == visitMax, ])]
        t[[i]] = seq(from = 0, to = 1, length.out = ncol(ccaFilter))[!is.na(ccaFilter[idFilter == ids[i] & visitFilter == visitMax, ])]
        xy[[i]] = y[i] * x[[i]]
      }
    }
    
    y = (y[lengths(t) >= 2])
    x = x[lengths(t) >= 2]
    xy = xy[lengths(t) >= 2]
    t = t[lengths(t) >= 2]
    
    detach(tmp)
  }

}

#######################################################################
#######################################################################
#######################################################################

#### settings for analysis
optns = list(
  FVEthreshold = FVEthreshold,
  userBwCov = 0,
  methodBwCov = "GCV",
  userBwMu = 0,
  methodBwMu = "GCV",
  dataType = "Sparse",
  error = TRUE, 
  fitEigenValues = FALSE,
  kernel = "epan",
  lean = FALSE,
  userCov = NULL,
  userMu = NULL,
  userSigma2 = NULL,
  verbose = FALSE,
  nRegGrid = 1000,
  maxK = 20
)

optnsQl = list(
  Kfold = 5,
  del = 1e-1,
  hLength = 20
)

# Dangerous! clear existing result
timeLlsfit = numeric(); timeQlfit = numeric(); timeFacefit = numeric()
timePaceLls = numeric(); iseePaceLls = reiseePaceLls = mspePaceLls = remspePaceLls = coverPaceLls = coverMissPaceLls =NULL
timePaceFace = numeric(); iseePaceFace = reiseePaceFace = mspePaceFace = remspePaceFace = coverPaceFace = coverMissPaceFace =NULL
timePleaLls = numeric(); iseePleaLls = reiseePleaLls = mspePleaLls = remspePleaLls = coverPleaLls = coverMissPleaLls = NULL
timePleaFace = numeric(); iseePleaFace = reiseePleaFace = mspePleaFace = remspePleaFace = coverPleaFace = coverMissPleaFace = NULL

# Check current progress
checkPaceLls = length(timePaceLls)
checkPaceFace = length(timePaceFace)
checkPleaLls = length(timePleaLls)
checkPleaFace = length(timePleaFace)

################ Replica ##################
if (simu) {
  mins = NULL
  maxs = NULL
  for (R in 1:RR) {
    mins[R] = min(unlist(tAll[[R]]))
    maxs[R] = max(unlist(tAll[[R]]))
  }
  comGrid = seq(from = max(mins), to = min(maxs), length.out = length(denseGrid))
  betaTrueComGrid = approx(x = denseGrid, y = betaTrue, xout = comGrid, method = 'linear', rule = 2)$y
}else if (sudoAnalyze) {
  # local linear smoother (LLS)
  ptm0 = proc.time()[3]
  llsObj = llsFit(y = y, x = x, t = t, optns, method = llsMethod)
  ptm1 = proc.time()[3]
  timeLlsfit = ptm1 - ptm0

  # quasi-likelihood fit following Zhou, Lin & Liang (2018)
  if (quasilk) {
    xyOld = xy[sampIdx]
    ptm0 = proc.time()[3]
    llsObj.xy = llsFit(y = NULL, x = xy, t = t, optns, method = llsMethod)
    qlObj = qlFit(y = y, xList = x, xyList = xy, tList = t, init.x = llsObj, init.xy = llsObj.xy, optns, optnsQl)
    ptm1 = proc.time()[3]
    timeQlfit = ptm1 - ptm0
  }
  
  # FACE Xiao et al. (2018)
  ptm0 = proc.time()[3]
  faceObj = FACE(y = y, x = x, t = t, optns)
  ptm1 = proc.time()[3]
  timeFacefit = ptm1 - ptm0
}

for (R in 1:RR) {
  if (simu) {
    sampIdx = sample(1:nSubject, round(nSubject * propTrain))
    yOld = yAll[[R]][sampIdx]
    xOld = xAll[[R]][sampIdx]
    xyOld = xyAll[[R]][sampIdx]
    tOld = tAll[[R]][sampIdx]
    yNew = yAll[[R]][-sampIdx]
    xNew = xAll[[R]][-sampIdx]
    xyNew = xyAll[[R]][-sampIdx]
    tNew = tAll[[R]][-sampIdx]
    yTrueNew = yTrueAll[[R]][-sampIdx]
    
    # local linear smoother (LLS)
    if (R > checkPaceLls | R > checkPleaLls) {
      ptm0 = proc.time()[3]
      llsObj = llsFit(y = yOld, x = xOld, t = tOld, optns, method = llsMethod)
      ptm1 = proc.time()[3]
      timeLlsfit[R] = ptm1 - ptm0
    }
    
    # quasi-likelihood fit following Zhou, Lin & Liang (2018)
    if (quasilk) {
      if (R > checkPaceLls | R > checkPleaLls) {
        ptm0 = proc.time()[3]
        llsObj.xy = llsFit(y = NULL, x = xyAll[[R]], t = tAll[[R]], optns)
        qlObj = qlFit(y = yOld, xList = xOld, xyList = xyOld, tList = tOld, init.x = llsObj, init.xy = llsObj.xy, optns, optnsQl)
        ptm1 = proc.time()[3]
        timeQlfit[R] = ptm1 - ptm0
      }
    }
    
    # FACE
    if (R > checkPaceFace | R > checkPleaFace) {
      ptm0 = proc.time()[3]
      faceObj = FACE(y = yOld, x = xOld, t = tOld, optns)
      ptm1 = proc.time()[3]
      timeFacefit[R] = ptm1 - ptm0
    }
    
    ## regression
    # PACE + LLS
    if (R > checkPaceLls){
      ptm0 = proc.time()[3]
      resPaceLls = PACE(yOld, xOld, tOld, optns, pMax = pMax, FVEthreshold = FVEthreshold, approxObj = llsObj, yNew, xNew, tNew, yTrueNew, 
                        betaTrue = approx(x = comGrid, y = betaTrueComGrid, xout = llsObj$regGrid, method = 'linear', rule = 2)$y, tuning = tuning
                        )
      ptm1 = proc.time()[3]
      timePaceLls[R] = ptm1 - ptm0
      
      iseePaceLls = rbind(iseePaceLls, resPaceLls$isee)
      reiseePaceLls = rbind(reiseePaceLls, resPaceLls$reisee)
      remspePaceLls = rbind(remspePaceLls, resPaceLls$remspe)
      coverPaceLls = rbind(coverPaceLls, resPaceLls$coverPr)
      coverMissPaceLls = rbind(coverMissPaceLls, c(resPaceLls$coverMissL, resPaceLls$coverMissR))
    }
    
    # PACE + FACE
    if (R > checkPaceFace){
      ptm0 = proc.time()[3]
      resPaceFace = PACE(yOld, xOld, tOld, optns, pMax = pMax, FVEthreshold = FVEthreshold, approxObj = faceObj, yNew, xNew, tNew, yTrueNew, 
                        betaTrue = approx(x = comGrid, y = betaTrueComGrid, xout = llsObj$regGrid, method = 'linear', rule = 2)$y, tuning = tuning
                        )
      ptm1 = proc.time()[3]
      timePaceFace[R] = ptm1 - ptm0
      
      iseePaceFace = rbind(iseePaceFace, resPaceFace$isee)
      reiseePaceFace = rbind(reiseePaceFace, resPaceFace$reisee)
      remspePaceFace = rbind(remspePaceFace, resPaceFace$remspe)
      coverPaceFace = rbind(coverPaceFace, resPaceFace$coverPr)
      coverMissPaceFace = rbind(coverMissPaceFace, c(resPaceFace$coverMissL, resPaceLls$coverMissR))
    }
    
    # PLEASS + LLS
    if (R > checkPleaLls){
      ptm0 = proc.time()[3]
      resPleaLls = PLEASS(yOld, xOld, tOld, optns, pMax = pMax, FVEthreshold = FVEthreshold, approxObj = llsObj, yNew, xNew, tNew, yTrueNew, 
                          betaTrue = approx(x = comGrid, y = betaTrueComGrid, xout = llsObj$regGrid, method = 'linear', rule = 2)$y, tuning = tuning
                          )
      ptm1 = proc.time()[3]
      timePleaLls[R] = ptm1 - ptm0
      
      iseePleaLls = rbind(iseePleaLls, resPleaLls$isee)
      reiseePleaLls = rbind(reiseePleaLls, resPleaLls$reisee)
      remspePleaLls = rbind(remspePleaLls, resPleaLls$remspe)
      coverPleaLls = rbind(coverPleaLls, resPleaLls$coverPr)
      coverMissPleaLls = rbind(coverMissPleaLls, c(resPleaLls$coverMissL, resPaceLls$coverMissR))
    }
    
    # PLEASS + FACE
    if (R > checkPleaFace){
      ptm0 = proc.time()[3]
      resPleaFace = PLEASS(yOld, xOld, tOld, optns, pMax = pMax, FVEthreshold = FVEthreshold, approxObj = faceObj, yNew, xNew, tNew, yTrueNew,
                           betaTrue = approx(x = comGrid, y = betaTrueComGrid, xout = faceObj$regGrid, method = 'linear', rule = 2)$y, tuning = tuning
                           )

      ptm1 = proc.time()[3]
      timePleaFace[R] = ptm1 - ptm0
      
      iseePleaFace = rbind(iseePleaFace, resPleaFace$isee)
      reiseePleaFace = rbind(reiseePleaFace, resPleaFace$reisee)
      remspePleaFace = rbind(remspePleaFace, resPleaFace$remspe)
      coverPleaFace = rbind(coverPleaFace, resPleaFace$coverPr)
      coverMissPleaFace = rbind(coverMissPleaFace, c(resPleaFace$coverMissL, resPaceLls$coverMissR))
    }
  }
  
  if (!simu) {
    sampIdx = sample(1:length(x), round(length(x) * propTrain))
    yOld = y[sampIdx]
    xOld = x[sampIdx]
    tOld = t[sampIdx]
    yNew = y[-sampIdx]
    xNew = x[-sampIdx]
    tNew = t[-sampIdx]
    
    if (!sudoAnalyze){
      # local linear smoother (LLS)
      if (R > checkPaceLls | R > checkPleaLls) {
        ptm0 = proc.time()[3]
        llsObj = llsFit(y = yOld, x = xOld, t = tOld, optns, method = llsMethod)
        ptm1 = proc.time()[3]
        timeLlsfit[R] = ptm1 - ptm0
      }
  
      # quasi-likelihood fit following Zhou, Lin & Liang (2018)
      if (quasilk) {
        if (R > checkPaceLls | R > checkPleaLls) {
          xyOld = xy[sampIdx]
          ptm0 = proc.time()[3]
          llsObj.xy = llsFit(y = NULL, x = xyOld, t = tOld, optns)
          qlObj = qlFit(y = yOld, xList = xOld, xyList = xyOld, tList = tOld, init.x = llsObj, init.xy = llsObj.xy, optns, optnsQl)
          ptm1 = proc.time()[3]
          timeQlfit[R] = ptm1 - ptm0
        }
      }
      
      # FACE
      if (R > checkPaceFace | R > checkPleaFace) {
        ptm0 = proc.time()[3]
        faceObj = FACE(y = yOld, x = xOld, t = tOld, optns)
        ptm1 = proc.time()[3]
        timeFacefit[R] = ptm1 - ptm0
      }
    }
    
    ## regression
    # PACE + LLS
    if (R > checkPaceLls){
      ptm0 = proc.time()[3]
      resPaceLls = PACE(yOld, xOld, tOld, optns, pMax = pMax, FVEthreshold = FVEthreshold, approxObj = llsObj, yNew, xNew, tNew, yTrueNew = NULL, betaTrue = NULL)
      ptm1 = proc.time()[3]
      timePaceLls[R] = ptm1 - ptm0
      
      mspePaceLls = rbind(mspePaceLls, resPaceLls$mspe)
      remspePaceLls = rbind(remspePaceLls, resPaceLls$remspe)
    }
    # PACE + FACE
    if (R > checkPaceLls){
      ptm0 = proc.time()[3]
      resPaceFace = PACE(yOld, xOld, tOld, optns, pMax = pMax, FVEthreshold = FVEthreshold, approxObj = faceObj, yNew, xNew, tNew, yTrueNew = NULL, betaTrue = NULL)
      ptm1 = proc.time()[3]
      timePaceFace[R] = ptm1 - ptm0
      
      mspePaceFace = rbind(mspePaceFace, resPaceFace$mspe)
      remspePaceFace = rbind(remspePaceFace, resPaceFace$remspe)
    }
    # PLEASS + LLS
    if (R > checkPleaLls){
      ptm0 = proc.time()[3]
      resPleaLls = PLEASS(yOld, xOld, tOld, optns, pMax = pMax, FVEthreshold = FVEthreshold, approxObj = llsObj, yNew, xNew, tNew, yTrueNew = NULL, betaTrue = NULL, tuning = tuning)
      ptm1 = proc.time()[3]
      timePleaLls[R] = ptm1 - ptm0
      
      mspePleaLls = rbind(mspePleaLls, resPleaLls$mspe)
      remspePleaLls = rbind(remspePleaLls, resPleaLls$remspe)
    }
    # PLEASS + FACE
    if (R > checkPleaFace){
      ptm0 = proc.time()[3]
      resPleaFace = PLEASS(yOld, xOld, tOld, optns, pMax = pMax, FVEthreshold = FVEthreshold, approxObj = faceObj, yNew, xNew, tNew, yTrueNew = NULL, betaTrue = NULL, tuning = tuning)
      ptm1 = proc.time()[3]
      timePleaFace[R] = ptm1 - ptm0
      
      mspePleaFace = rbind(mspePleaFace, resPleaFace$mspe)
      remspePleaFace = rbind(remspePleaFace, resPleaFace$remspe)
    }
    
  }
  
  file = ifelse(simu,
                paste0('Rimage/',
                      'Simu', simuCase, '_',
                      nSubject, 'nSubject_',
                      J, 'eigen_',
                      FVEthreshold, 'FVE_',
                      RR, 'repeats_',
                      optns$nRegGrid, 'nRegGrid_',
                      SNRy, 'SNR_',
                      propTrain * 100, 'train.RData'),
                 paste0('Rimage/',
                       'Real', realCase, '_', 
                       FVEthreshold * 1e4, 'FVE_',
                       RR, 'repeats_', 
                       optns$nRegGrid, 'nRegGrid_',
                       propTrain * 100, 'train.RData')
  )
  
  if (R %% 10 == 0){
    if (sudoAnalyze == F)
      remove(llsObj, llsObj.xy, qlObj, faceObj) # to save storage space
    save.image(file = file)
  }
  if (R %% 40 == 0){
    cat(R, '\n')
  }else cat(R)
}

################ Plot errors ##################
source("supportFun.R")
if (tuning){
  ncolPlot = 1
  colNames = c(
    "PLEASS.L",
    "PACE.L",
    "PLEASS.F",
    "PACE.F"
  )
}else{
  ncolPlot = min(10, maxK)
  colNames = c(
    apply(expand.grid("PLEASS.L", 1:ncolPlot), 1, paste0, collapse='.'),
    apply(expand.grid("PACE.L", 1:ncolPlot), 1, paste0, collapse='.'),
    apply(expand.grid("PLEASS.F", 1:ncolPlot), 1, paste0, collapse='.'),
    apply(expand.grid("PACE.F", 1:ncolPlot), 1, paste0, collapse='.')
  )
}
if (simu){
  reiseeMat = creatErrMat(colNames,
                          list(
                            as.matrix(reiseePleaLls[, 1:ncolPlot]), 
                            as.matrix(reiseePaceLls[, 1:ncolPlot]), 
                            as.matrix(reiseePleaFace[, 1:ncolPlot]), 
                            as.matrix(reiseePaceFace[, 1:ncolPlot])
                            )
                          )
  boxplotErr(nrow(reiseeMat), reiseeMat, type = "ReISEE")
  coverMat = creatErrMat(colNames,
                         list(
                           as.matrix(coverPleaLls[, 1:ncolPlot]),
                           as.matrix(coverPaceLls[, 1:ncolPlot]),
                           as.matrix(coverPleaFace[, 1:ncolPlot]),
                           as.matrix(coverPaceFace[, 1:ncolPlot])
                           )
                         )
  boxplotErr(nrow(coverMat), coverMat, type = "Coverage Percentage")
  
}else{
  remspeMat = creatErrMat(colNames,
                          list(
                            as.matrix(remspePleaLls[, 1:ncolPlot]), 
                            as.matrix(remspePaceLls[, 1:ncolPlot]),
                            as.matrix(remspePleaFace[, 1:ncolPlot]), 
                            as.matrix(remspePaceFace[, 1:ncolPlot])
                          )
  )
  boxplotErr(nrow(remspeMat), remspeMat, type = "ReMSPE")
}

