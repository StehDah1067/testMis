globalVariables(c("stage"))

xPrint <- function(x, dec = 1) sprintf(paste("%.", dec, "f", sep = ""), x)

CondRelSurv_NE_Pop <- function(rs_0, ivall = c(95, 103)) { # rs_0 = x
  rs = se = years = cond_rs = cond_se = NULL
  tb <- rs_0[rs > 0 & se > 0 & cond_rs > 0]
  if (tb[, min(years)] == 1) { # min im 1.Jahr
    if (tb[1:min(10, nrow(tb))][, min(cond_rs + cond_se)] > 99 &
        tb[, max(rs)] < 110) {
      return(FALSE)
    }
  } # surv > ref.pop
  if (tb[, min(cond_rs)] > ivall[1] &
      tb[, max(cond_rs)] < ivall[2]) {
    return(FALSE)
  }
  if (nrow(tb) > 4) { # mindestens 4 Stützstellen für Spline
    uu <- rsSpline(tb[, .(years, rs)])
    i1 <- floor(min_Rs(uu))
    m1 <- max_Rs(uu, 1)
    if (tb[m1, rs - se] < 103 & tb[i1, rs + se] > 95) {
      return(FALSE)
    }
  } #  plotSpline(tb)
  return(TRUE)
}
min_Rs <- function(uu) {
  if (nrow(uu$mins) > 0) {
    return(min(uu$mins$x[1]))
  }
  if (uu$d012$y[1] <= uu$d012$y[nrow(uu$d012)]) { # stetig steigend
    return(1)
  }
  return(max(uu$d012$x)) # letztes Jahr
}
max_Rs <- function(uu, i1) { # i1=1
  y = x = NULL
  if (nrow(uu$maxs) > 0) {
    ix <- which(uu$maxs$x > i1)
    m1 <- max(uu$maxs[ix, y])
    if (uu$d012[nrow(uu$d012), y] > m1) { # maximum am Ende
      return(uu$d012[nrow(uu$d012), floor(x)])
    }
    return(uu$d012[y == m1, floor(x)])
  }
  return(uu$d012[nrow(uu$d012), x])
}

test_rsIncrease <- function(tab, rs_0, opt) { # tab=dp$tab
  maxRs_gt <- function(tab, uu, vv) {
    s1 <- max_Rs(uu, 1)
    tab[s1, rs - se] > vv
  }
  rs = se = years = NULL
  tab <- data.table(tab[rs > 0 & se > 0, ])
  uu <- rsSpline(tab[, .(years, rs)])
  if (maxRs_gt(tab, uu, 110)) {
    return(T)
  }
  # print(showTab(tab))
  #  plotSpline(tab[,.(stage,years,rs)],60)
  i1 <- floor(min_Rs(uu))
  m1 <- max_Rs(uu, i1)
  #    print(paste(m1,i1))
  if (m1 <= i1) { # max vor min
    return(F)
  }
  x <- tab[i1:nrow(tab), ]
  if (nrow(x) < 4) { # kein spline möglich
    return(F)
  }
  lo <- x[1, rs + 1.00 * se]
  i_max <- m1 - i1 + 1
  if (CondRelSurv_NE_Pop(x[c(1:i_max)], ivall = c(95, 103))) {
    if (x[, .N] == 0) {
      return(FALSE)
    }
    hi <- x[i_max, rs]
    return(hi > lo) #
  } else {
    return(FALSE)
  }
}

rsSpline <- function(tb) { # tb = tb[,.(years,rs)]
  findExtrema <- function(d012) {
    i <- 1
    ix <- 0
    n1 <- nrow(d012) - 1
    while (i < n1) {
      if ((d012$y1[i] > 0 & d012$y1[i + 1] < 0) |
          (d012$y1[i] < 0 & d012$y1[i + 1] > 0)) {
        ix <- c(ix, i)
      }
      i <- i + 1
    }
    return(ix[-1])
  }
  y2 = NULL
  tb$rs <- as.numeric(tb$rs)
  sp1 <- stats::smooth.spline(tb)
  Sx <- seq(1, nrow(tb), 0.1)
  d0 <- data.frame(stats::predict(sp1, Sx))
  d1 <- data.frame(stats::predict(sp1, Sx, deriv = 1))
  d2 <- data.frame(stats::predict(sp1, Sx, deriv = 2))
  setnames(d1, "y", "y1")
  setnames(d2, "y", "y2")
  d12 <- merge(d1, d2, by = "x")
  d012 <- setDT(merge(d0, d12, by = "x"))
  extrema <- findExtrema(d012)
  mins <- d012[extrema][y2 > 0]
  # if(d012$y[1] < mins$y[1])
  #   mins = rbind(d012[1],mins) #min. am linken Rand
  maxs <- d012[extrema][y2 < 0]
  #  plotSpline(tb[,.(years,rs)])
  return(list(d012 = d012, sp1 = sp1, mins = mins, maxs = maxs)) # ,Extrema=extrema))
}
#
# ddI.logRegPredict <- function(alive, am, opt) { # 	am=m2
#   alive$prob <- predict(am, alive, type = "response")
#   alive$zoutl <- ifelse(alive$prob > opt$pCut, 1, 0)
#   return(alive)
# }

ddI.predict <- function(dp, z1, opt, sp, rs_0) {
  prob <- zoutl <- death <- NULL
  dead <- z1[death == 1]
  dead[, zoutl := 0]
  dead[, prob := 0]
  alive <- z1[death == 0]
  if (!is.null(dp$m1)) {
    alive$prob <- stats::predict.glm(dp$m1, newdata = alive, type = "response")
    alive$zoutl <- ifelse(alive$prob > opt$pCut, 1, 0)
  } else {
    alive[, zoutl := 0]
    alive[, prob := 0]
  }
  dp$data <- rbind(alive, dead)
  zz <- subset(dp$data, zoutl == 0)
  dp$tab <- myRelSurv(zz, opt, sp)
  return(dp)
}


scoreFuAgeSurvy <- function(dp, z1, rs_0, sp, opt) {
  myScore <- function(alive) {
    a1 <- alive[order(fu.age, decreasing = TRUE), .(id, fu.age, fu.time)] # [1:1000,]
    a2 <- alive[order(fu.time, decreasing = TRUE), .(id)] # [1:1000,]
    # bb = a1[a2,on='id',all=FALSE]
    bb <- merge(a1, a2, by = "id", all.x = FALSE, all.y = FALSE)
    bb <- data.table(merge(a1, a2, by = "id", all.x = FALSE, all.y = FALSE))
    bb$norm <- sqrt(bb$fu.age^2 + bb$fu.time^2)
    bb <- bb[order(norm), ]
    bb$nr <- c(1:nrow(bb))
    return(bb)
  }
  zoutl <- prob <- death <- fu.time <- fu.age <- NULL
  dp$data[, zoutl := 0]
  dp$data[, prob := NA]
  dp$tab <- rs_0
  alive <- data.table(subset(dp$data, death == 0))
  deaths <- data.table(subset(dp$data, death == 1))
  bb <- myScore(alive)
  q1 <- 0
  while (test_rsIncrease(dp$tab, rs_0, opt)) {
    ii <- which(bb$nr >= stats::quantile(bb$nr, opt$quant - q1))
    alive[id %in% bb$id[ii], ]$zoutl <- 1
    dp$data <- rbind(alive, deaths)
    z2 <- dp$data[zoutl == 0, ]
    x1 <- myRelSurv(z2, opt, sp)
    dp$tab <- x1
    q1 <- q1 + 0.001
    if(dp$scoring==0) dp$scoring <- 1
    if(PercentMissedDeath(dp) > opt$md_cutoff) {
      break
    }
  }
  dp$quant <- opt$quant - q1
  return(dp)
}

PercentMissedDeath <- function(dp) {
  zoutl = NULL
  death <- nrow(subset(dp$data, death == 1))
  missed_death <- nrow(subset(dp$data, zoutl == 1))
  if (death > 0) pp <- missed_death * 100 / (death + missed_death) else pp <- 0
  return(pp)
}

myRelSurv <- function(z1, opt, sp) {
#  z1[, vitstat := death + 1]
  z1[,':=' (vitstat = z1$death+1)]
  rs_0 <- periodR::period(data.frame(z1),
                          k = opt$survYears,
                          sp$males, sp$females,
                          perbeg=opt$perbeg,
                          opt$lastYear,
                          method = "edererII")
  data.table(
    stage = unique(z1$stage),
    years = c(1:length(rs_0$rel.surv)),
    rs = rs_0$rel.surv,
    se = rs_0$rel.surv.stderr,
    cond_rs = rs_0$cond.rel.surv,
    cond_se = rs_0$cond.rel.surv.stderr
  )
}

init_dp <- function(z1, deaths, opt) {
  dp <- list(stage = as.character(unique(z1$stage)))
  dp$CondRelSurv_ne_Pop <- 0
  dp$md_gt_cutoff <- 0
  dp$scoring <- 0
  dp$rs_0Increase <- 0
  dp$reg_sig <- -1
  dp$p_reg_md <- -1
  dp$quant <- opt$quant # 1-0.005
  dp$ageCut <- as.numeric(stats::quantile(deaths$fu.age, dp$quant))
  dp$survCut <- as.numeric(stats::quantile(deaths$fu.time, opt$quant))
  return(dp)
}


findCutPoints <- function(z1, dp, opt, sp, rs_0, quant) { # q1=z1;quant = opt$quant-q1
  death <- NULL
  brglm_NotConverge <- function(dead) {
    tt <- tryCatch(
      fit <- brglm::brglm(formula = zoutl ~ fu.age + fu.time, family = "binomial", data = dead),
      error = function(e) e,
      warning = function(w) w
    )
    u <- grep("Iteration limit reached", tt$message)
    if (length(u) > 0) {
      return(TRUE)
    }
    #    print(tt$message)
    return(FALSE)
  }
  setOutlier <- function(dp, dead) {
    q.age <- ifelse(dead$fu.age > dp$ageCut, 1, 0)
    q.surv <- ifelse(dead$fu.time > dp$survCut, 1, 0)
    zoutl <- ifelse(q.surv == 1 & q.age == 1, 1, 0)
    return(zoutl)
  }
  if (test_rsIncrease(tab = rs_0, rs_0, opt)) {
    dp$rs_0Increase <- 1
    dead <- subset(z1, death == 1)
    dead$zoutl <- setOutlier(dp, dead)
    q2 <- 0.0005
    while (brglm_NotConverge(dead)) { # reduziert das Quantil bis Glm konvergiert
      quant <- quant - q2
      dp$ageCut <- as.numeric(stats::quantile(dead$fu.age, quant))
      dp$survCut <- as.numeric(stats::quantile(dead$fu.time, quant))
      dead$zoutl <- setOutlier(dp, dead)
    }
    oldw <- getOption("warn")
    options(warn = -1)
    dp$m1 <- brglm::brglm(formula = zoutl ~ fu.age + fu.time, family = "binomial", data = dead)
    options(warn = oldw)
  }
  dp <- ddI.predict(dp, z1, opt, sp, rs_0)
  #  print(table(dp$data$zoutl))
  return(dp)
}

iterateMissedDeaths <- function(dp, rs_0, z1, deaths, sp, opt) {
  showCut <- function(dp) paste("age", dp$ageCut, ", surv", dp$survCut, "quant:", dp$quant)
  #  dp = ddi.init_dp(z1,deaths,opt)
  q1 <- 1 - dp$quant
  delta <- 0.005
  dp$CondRelSurv_ne_Pop <- 1
  dp <- findCutPoints(z1, dp, opt, sp, rs_0, opt$quant)
  while (test_rsIncrease(dp$tab, rs_0, opt)) {
    #     	print(opt$quant-q1)
    dp$ageCut <- as.numeric(stats::quantile(deaths$fu.age, opt$quant - q1))
    dp$survCut <- as.numeric(stats::quantile(deaths$fu.time, opt$quant - q1))
    dp <- findCutPoints(z1, dp, opt, sp, rs_0, opt$quant - q1)
    q1 <- q1 + delta
    dp <- ddI.predict(dp, z1, opt, sp, rs_0)
    #  showTab(dp$tab) #dp$tab[,.(stage,years,rs,cond_rs,cond_se)]
    #     plotSpline(dp$tab[,.(stage,years,rs)])
    if (PercentMissedDeath(dp) > opt$md_cutoff) {
      break
    }
  }
  dp$quant <- dp$quant - (q1 - delta)
  dp$reg_sig = regCoefSigni(dp$m1)
  #  showCut(dp)
  #  dp$data[death==1,.N,zoutl]
  return(dp)
}

regCoefSigni <- function(m1){#m1 =dp$m1
  if(is.null(m1))
    return(-1)
  uu <- coefficients(summary(m1))
  setDT(as.data.frame(uu))
  bb <- uu[,'Pr(>|z|)']
  if(length(bb) < 3)
    return(-1)
  signi <- 0
  if(bb[2]  <= .05 | bb[3] <= .05)
    signi <- 1
  return(signi)
}

tabValues = function(dp,z1,rs_0){
  dp$tab[, icd := paste(sort(unique(z1$icd)),collapse = ',')]
  dp$tab[, sex:= paste(sort(unique(z1$sex)),collapse = ',')]
  dp$tab[, agr := paste(sort(unique(z1$agr)),collapse = ',')]
  dp$tab[, rs_before := xPrint(rs_0$rs)]
  #  dp$tab[, rs_all := xPrint(rs_0$rs)]
  dp$tab[, rs_before.se := xPrint(rs_0$se)]
  dp$tab[, rs_after := xPrint(rs)]
  dp$tab[, rs_after.se := xPrint(se)]
  dp$tab <- dp$tab[, .(icd,sex,agr,stage, years, rs_before, rs_before.se, rs_after, rs_after.se)]
  return(dp)
}

#z1=dat[stage=='I'];
dedectImmortals <- function(z1, opt, sp){#z1 = dat
  rs <- se <- rs_all <- rs_all.se <- NULL
  years <-  NULL
#  print(z1[,.N,stage])
  if (nrow(z1) == 0) {
    return(paste("keine Daten in z1"))
  }
  cls = c('id', 'icd','stage','sex','diagage','agr','dm','dy','fm','fy','death','fu.age','fu.time')
  cls = names(z1)
  z1 <- z1[,cls,with = FALSE]
  percentDeaths <- function(z1) sum(z1$death == 1) * 100 / nrow(z1)
  z1$zoutl<- 0
  deaths <- z1[with(z1,death == 1)]
  dp <- init_dp(z1, deaths, opt)
  rs_0 <- myRelSurv(z1, opt, sp)
  if (CondRelSurv_NE_Pop(rs_0)) {
    dp <- iterateMissedDeaths(dp, rs_0, z1, deaths, sp, opt)
    dp <- ddI.predict(dp, z1, opt, sp, rs_0)
    if((percentDeaths(z1) < 10.0)  |
       (PercentMissedDeath(dp) >= opt$md_cutoff) |
       (dp$reg_sig < 1)
    ){
      if(PercentMissedDeath(dp) > opt$md_cutoff)
        dp$p_reg_md <- PercentMissedDeath(dp)
      dp <- scoreFuAgeSurvy(dp, z1, rs_0, sp, opt)
    }
  } else {
     dp <- ddI.predict(dp, z1, opt, sp, rs_0)
  }
  dp <- tabValues(dp,z1,rs_0)
  return(dp)
}

result <- function(dp,opt){#dp = xx$dp
  getAnalysisParam <- function(d0){#d0 = dp[[1]]
    data.table(
      agr = unique(d0$data$agr),
      sex = unique(d0$data$sex),
      stage = d0$stage,
      CondRelSurv_ne_Pop = d0$CondRelSurv_ne_Pop,
      reg_sig = d0$reg_sig,
      p_reg_md = d0$p_reg_md,
      scoring = d0$scoring,
      md_gt_cutoff = d0$md_gt_cutoff,
      rs_0Increase = d0$rs_0Increase,
      ageCut = round(d0$ageCut,1),
      survCut = round(d0$survCut,1),
      quant = d0$quant
    )
  }
  getRegCoef <- function(dp,dd){#ii=2
    if(!is.null(dp$m1)){
      mod = summary(dp$m1)
      if(regCoefSigni(dp$m1)) si = 1 else si = 0
      rc <- as.data.table(mod$coefficients)[,c(1,2,4)]
      rc[,varnames := row.names(mod$coefficients)]
      rc[,signi := si]
    } else{
      rc = data.table(Estimate=0, 'Std. Error'=0,
                 'Pr(>|z|)'=1)
      rc[,varnames := 'NULL']
      rc[,signi := -1]
    }
    rc[,stage := unique(dp$data$stage)]
    rc[,agr := unique(dd$agr)]
    rc[,sex := unique(dd$sex)]
    rc <- rc[,c('sex','agr','stage','varnames',
          'Estimate','Std. Error','Pr(>|z|)')]#,'signi')]
    return(rc)
  }
  dd <- dp$data;
  oo <- getAnalysisParam(dp)
  tab <- data.frame(dp$tab)
  regCoef <- getRegCoef(dp,dd)
  dd[,vitstat := NULL]
  setnames(dd,'zoutl','missedDeath')
  dd[,death := ifelse(missedDeath == 1,1,death)]
  opt$quant <-NULL;opt$pCut <- NULL
  Param <- list( opt = opt, runTime=oo)
  dd = data.frame(dd)
#  print(Param)
  return(list(dat=dd,survTab = tab,Param=Param,regCoef=regCoef))
}


exitTime <- function(t1,opt){#t1 = z0
  Ende  <- as.Date(paste(opt$lastYear,12,31,sep='-'))
  t1[,fd := as.Date(paste(fy,fm,'15',sep='-'))]
#  t1[,fd := ifelse(fd > as.Date(paste(opt$lastYear,12,15,sep='-')),Ende,fd)]
  t1[death == 0,fd := ifelse(fd >= as.Date(paste(opt$lastYear,12,15,sep='-')),Ende,fd)]
  t1[death == 1,fd := ifelse(fd > as.Date(paste(opt$lastYear,12,15,sep='-')),Ende,fd)]
  t1[,dd := as.Date(paste(dy,dm,'15',sep='-'))]
  t1$fd <- ifelse(t1$fd == t1$dd,t1$fd+15,t1$fd)
  t1$death <- ifelse(t1$fd >= Ende,0,t1$death)
  #t1$vitstat1 <- t1$death+1
  class(t1$fd) <- "Date"
  t1[,fu.time := as.numeric(fd - dd)/365.24]
  t1[,fu.age := diagage + fu.time]
  return(t1[])
}
showMessages = function(dat){
  if(nrow(dat) == 0){
    e <- simpleError("No patient data, empty data.frame!")
    tryCatch(stop(e)#,  finally = print("")
        )
  }

  if(nrow(dat) < 500)
    withCallingHandlers({ warning("Too small stratum size"); 1+2 },
        warning = function(w) {
          print(paste('Too small stratum size! N =',nrow(dat)))
        }
    )
#  return()
}


#' Classify missed deaths in cancer registry data
#'
#' Classify potential missed deaths by comparing
#' survival times of living and deceased persons.
#' @details
#' This function compares length of follow-up time
#' and age at end of follow-up between living and
#' deceased persons in the input dataset (\code{dat}).
#' To increase the validity of these comparisons,
#' data should be stratified according to
#' prognostic factors. Currently, prognostic factors
#' reflect items typically available in epidemiological
#' cancer registry data: Sex, Age at diagnosis,
#' Diagnosis and Stage. The input dataset should include
#' only one combination of the values of these variables,
#' e.g., Women with stage I melanoma of the skin between
#' 50 and 69 years old at diagnosis. However, to strengthen
#' the validity of the classification, we recommend a
#' stratum size of at least 500 persons.
#'
#' Because the algorithm implemented here requires reliable
#' information on follow-up time for deceased persons, cases
#' with unknown date of diagnosis (e.g., Death Certificate-Only
#' [DCO] cases) should be excluded. Cases lost to follow-up must
#' be censored at the last date on which vital status was known.
#'
#' The input dataset should include one entry per person.
#' Therefore, if the original data include multiple cases
#' per person, only the most recent case should be included.
#' If there are multiple, simultaneously diagnosed cases, the
#' most lethal diagnosis (determined by the average prognosis
#' in the respective population or from the literature) should
#' be included.
#'
#' Within \code{classify}, relative survival is estimated by the
#' period approach (Gondos A, Brenner H, Holleczek B (2009)).Per
#' default we chose the last 5 years before follow-up end
#' \code{fu_end} as the analysis period. The width of the
#' analysis period can be changed by the parameter
#' \code{periodWidth}, and follow-up end is always set by
#' \code{fu_end}.
#'
#'
#' Whenever 1-year conditional relative survival is between 95\%
#' and 103\% in all follow-up years, we assume that survial in the
#' patient group is similar to survival in the general population.
#' In this case the search for missed deaths is not carried out.
#'
#'
#' @param dat \code{data.frame} consisting of 11 variables:\cr
#'  \code{id} - (int) unique identifier for each data line,\cr
#'  \code{sex} -  (int) 1 = male, 2 = female, \cr
#'  \code{icd} - (chr) identifier for type of cancer (e.g.,
#'  according to ICD-10),\cr
#'  \code{stage} - (factor) cancer stage (e.g., "I", "II",
#'  "III", "IV", "Missing" for UICC stages),\cr
#'  \code{diagage} - (int) age at diagnosis in years,\cr
#'  \code{agr} - (factor) age group at diagnosis,\cr
#'  \code{dm} - (int) month of diagnosis, a number in \code{1:12},\cr
#'  \code{dy} - (int) calendar year of diagnosis (e.g., 2010),\cr
#'  \code{fm} - (int) month of follow-up end, a number in \code{1:12}, \cr
#'  \code{fy} - (int) calendar year of follow-up end (e.g., 2012),\cr
#'  \code{death} - (int) death at follow-up end (0 = alive, 1 = dead).\cr
#' To improve classification, it is recommended to use only one
#' combination of icd, sex, age group and stage.
#'
#' @param sp \code{list} of two \code{data.frame}s named
#' \code{females} resp. \code{males}. Each \code{data.frame}
#' contains age- and calendar year-specific survival
#' probabilities of women resp. men in the underlying population.
#' Each column should represent one calendar year and should
#' consist of a numeric vector of length 100 providing  1-year
#' survival probabilities from ages 0 to 99 years. Each column
#' should be named as the year it represents (e.g., \code{"2004"}).
#' If the data in \code{dat} include one or more years for which
#' there is no corresponding column, the column corresponding to
#' the closest year prior to the missing year(s) is used.
#'
#' @param fu_end Integer, indicating the last year of follow-up
#' (eg. 2018). The function assumes that follow-up is complete
#' through December 31st of this year. Patients with date of
#' deaths after this date are regarded as survivors.
#'
#' @param survYears Integer, indicating maximum follow-up time
#' in years (e.g., \code{survYears = 30}). Cases diagnosed in the
#' time span \cr\code{fu_end - survYears - periodWidth + 1} to
#' \code{fu_end} \cr are included in relative survival calculations.
#'
#' @param md_cutoff Numeric, (default \code{value md_cutoff =
#' 3.0}) indicating the maximum accepted percentage of
#' classified missed deaths. When the percentage of missed deaths
#' classified by logistic regression exceeds \code{md_cutoff},
#' the scoring algorithm is applied. Scoring is also stopped
#' when \code{md_cutoff} is exceeded.
#'
#' @param periodWidth Numeric, (default value \code{periodWidth = 5})
#'  indicating the width of the analysis period for relative
#'  survival estimation in years.\cr
#'  For example, if \code{fu_end = 2016} and \code{periodWidth = 2}
#'  the analysis period starts on 2015/01/01 and ends on 2016/12/31.
#'
#' @return \code{classify} returns a \code{list} containing
#'  four objects (\code{dat}, \code{survTab}, \code{Param} and
#'  \code{regCoef})\cr

#'
#' \itemize{
#'    \item{\code{dat}:\cr \code{data.frame} with the input data and four}
#' additional columns:
#'  \itemize{
#'  \item{\code{fu.age}: Age at follow-up end in years}
#'  \item{\code{fu.time}: Follow-up time in years}
#'  \item{\code{missedDeath}: 0: not a missed deaths;
#'  1: a missed deaths}
#'  \item{\code{prob}: propability of beeing a
#'     missed deaths estimated by the logistic regression.
#'     If the scoring algorithm is appield \code{prob} is
#'     set to \code{NA}}
#'  }
#'    \item{\code{survTab}:\cr \code{data.frame}
#'    contaning estimates of relative survival before
#'    and after excluding 'missed deaths'. \code{survTab}
#'    is used by \code{plotSurv()} for plotting survival rates.}
#'    \item{\code{Param}:\cr \code{list} with input and
#'    runtime parameters. \code{Param} is displayed by \code{classify.Param()}
#'    }
#'    \item{\code{regCoef}:\cr \code{data.frame} with coefficients
#'    of the logistic regression. If \code{Pr(>|z|) > 0.05} for
#'    both fu.age and fu.time, the logistic regression
#'    was regarded as not significant and the scoring algorithm
#'    was applied for classification of missed deaths.}
#'}
#'
#' @export
#' @import data.table
#' @importFrom Rdpack reprompt
#'
#'
#' @examples
#' data(C61)
#' data(sp)
#' xx = classify(C61,sp,fu_end=2018,survYears=30,md_cutoff=4.0,periodWidth=4)
#' classify.summary(xx)
#' classify.plotSurv(xx)
#' classify.parameter(xx)
#' classify.regCoef(xx)
#' head(xx$dat)
#'
#' @references
#'  \insertRef{Dahm2023}{misseddeaths}\cr\cr
#'  \insertRef{Gondos2009}{misseddeaths}
#'
#'@seealso
#'  \code{\link{classify.parameter}}
#'  \code{\link{classify.plotSurv}}
#'  \code{\link{classify.regCoef}}
#'  \code{\link{classify.summary}}
#'

classify <- function(dat,sp,fu_end,survYears,
                     md_cutoff=3.0,periodWidth=5) { # dat = z0;fu_end=2019;survYears=35;md_cutoff=3.0;periodWidth=5
#warning about minimum stratum size
  showMessages(dat)
  opt <- list(survYears=survYears,
              lastYear = fu_end,
              md_cutoff = md_cutoff,
              periodWidth = as.integer(periodWidth),
              perbeg = fu_end-periodWidth+1,
              quant = 0.999,
              pCut = 0.9) #cut point logistic regression
  dat <- setDT(dat)
  dat <- exitTime(dat,opt)
  dp <- dedectImmortals(dat, opt, sp)
  rr = result(dp,opt)
  print(classify.summary(rr))
  return(rr)
}
