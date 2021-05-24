################################################################################
####################SOURCE CODE#################################################
################################################################################

library(dplyr)
library(sars)
library(compute.es)
library(RColorBrewer) #colour pallete
col <- brewer.pal(4, "Set1")
library(ggplot2)
library(gridExtra)


################################################################
####ANCOVA AND EFFECT SIZES############################
##########################################################

ANC <- function(dat, A, S, G, plot = FALSE){
  
  #First get z and sig for each of the three
  modA <- lm(log10(s + 0.1) ~ log10(a), data = A)
  modS <- lm(log10(s + 0.1) ~ log10(a), data = S)
  modG <- lm(log10(s + 0.1) ~ log10(a), data = G) 
  
  #return model fits if plot = TRUE
  if (plot){
    ll <- list(modG, modS, modA)
    return(ll)
  }
  
  zpA <- summary(modA)$coefficients[2, c("Estimate", "Pr(>|t|)")]
  zpS <- summary(modS)$coefficients[2, c("Estimate", "Pr(>|t|)")]
  zpG <- summary(modG)$coefficients[2, c("Estimate", "Pr(>|t|)")] 
  
  res1 <- round(c(zpA, zpS, zpG), 2)
  
  #get the data into format for ANCOVA - we just want
  #the S and G data
  ds <- filter(dat, factor != "a")
  
  #Then do the ancova
  mod1 <- aov(log10(s + 0.1) ~ log10(a) * factor, data=ds)
  mod2 <- aov(log10(s + 0.1) ~ log10(a) + factor, data=ds)
  anc <- anova(mod1,mod2)
  
  res2 <- round(c(anc$F[2], anc$`Pr(>F)`[2]), 2)
  
  ##FR alternative R method
  mod11 <- lm(log10(s + 0.1) ~ log10(a) * factor, data=ds)
  ss <- anova(mod11)[,2]
  ss_sum <- sum(ss)
  R2_ind <- ss[1:3] / ss_sum #calculate R2 for each indiv term
  #check this matches with model R2
 if (!round(sum(R2_ind), 4) == round(summary(mod11)$r.squared, 4)){
   stop("R2 does not match")
 }
  #check extracting R2 for interaction term
  if (!rownames(anova(mod11))[3] == "log10(a):factor"){
    stop("anova rownames do not match")
  }
  R2_int <- R2_ind[3] #the R2 of just interaction term
  R_2 <- sqrt(R2_int) #the alternative R value for Cohen D
  
  #check extracting R2 for covariate term
  if (!rownames(anova(mod11))[1] == "log10(a)"){
    stop("anova rownames do not match_B")
  }
  R2_cov <- R2_ind[1] #the R2 of just covariate term
  R_3 <- sqrt(R2_cov) #the alternative R value for Cohen D
  
  #4th Method: the multiple correlation of the mod1
  R_4 <- sqrt(summary(mod11)$r.squared)
  
  ##Cohen's D and P-value
  #use same data as in lm (so log10 with + 0.1 for richness)
  R <- cor(log10(ds$a), log10(ds$s + 0.1), method = "pearson")
  #R <- cor(A$a, A$s + 0.1), method = "pearson")
  n <- nrow(ds) / 2
  
  #using R as simply correlation between a and s
  ES <- a.fes(f = anc$F[2], n.1 = n, n.2 = n, R = R, q = 1, verbose = FALSE)
  res3 <- c(ES$d, ES$pval.d, ES$var.d)
  
  #using R as sqrt(R2) of interaction term
  ES_2 <- a.fes(f = anc$F[2], n.1 = n, n.2 = n, R = R_2, q = 1, verbose = FALSE)
  res3_2 <- c(ES_2$d, ES_2$pval.d, ES_2$var.d)

  #Using R as multiple correlation from interaction model
  ES_4 <- a.fes(f = anc$F[2], n.1 = n, n.2 = n, R = R_4, q = 1, verbose = FALSE)
  res3_4 <- c(ES_4$d, ES_4$pval.d, ES_4$var.d)
  
  #Using t rather than F to generate the effect size
  #calculate t-manually and check it matches with t-value of interaction
  #term in mod11
  slope_g <- summary(modG)[["coefficients"]][,"Estimate"][2]
  slope_s <- summary(modS)[["coefficients"]][,"Estimate"][2]
  se_g <- summary(modG)[["coefficients"]][,"Std. Error"][2]
  se_s <- summary(modS)[["coefficients"]][,"Std. Error"][2]
  t_obs <- (slope_s - slope_g) / sqrt(se_g^2+se_s^2)
  t_value <- summary(mod11)[["coefficients"]][,"t value"][4]
  if (round(t_obs, 3) != round(t_value, 3)) stop("brown paper bag")
  #use tes rather than a.fes to calculated effect size
  ES_5 <- tes(t_value, n, n, verbose = FALSE)
  res3_5 <- c(ES_5$d, ES_5$pval.d, ES_5$var.d)
  
  res <- as.vector(c(res1, res2, res3, res3_2, res3_4, res3_5))
  names(res) <- c("zA", "pA", "zS", "pS", "zG", "pG", "F", "Pf", 
                  "D_1", "Pd_1", "varD_1", "D_2", "Pd_2", "varD_2",
                  "D_4", "Pd_4", "varD_4", "D_5", "Pd_5", "varD_5")
  
  return(res)
}

############################################################
##calculate pooled mean effect size estimate (d+)
###########################################################

#direct weights defined as the inverse of the variance of d 

dPlus <- function(d, var){
  #d is the individual effect sizes; var is the individual ES variances
  
  inv_var <- 1 / var
  
  #formula in Field (2000)
  dplus <- weighted.mean(d, inv_var)
  
  #weighted mean func gives same as this, which is the actual formula in Field (2000)
  if (dplus != (sum(d / var)) / (sum(1/var))) stop("yoshemi")

  #standard deviation of mean and Z formulas taken from Field (2000)
  sdM <- sqrt((sum(inv_var)) ^ - 1)
  
  Z = dplus / sdM
  #P-value formula used inside metagen R package
  P = 2 * pnorm(abs(Z), lower.tail = FALSE)
  
  #standard error of dplus (which is a weighted mean)
  #formula taken from Wikipedia
  #https://en.wikipedia.org/wiki/Weighted_arithmetic_mean:
  #'The standard error of the weighted mean (unit input variances), 
  #can be shown via uncertainty propagation to be:' 
  SE <- (sqrt(sum(inv_var)))^-1
  
  #confidence interval formula from inside metagen R package (standard CI formula)
  lower  <- dplus - qnorm(1 - 0.05 / 2) * SE
  upper  <- dplus + qnorm(1 - 0.05 / 2) * SE
  
  resPlus <- round(c("dplus" = dplus, "Z" = Z, "P" = P, "CIL" = lower, "CIU" = upper), 2)
  return(resPlus)
}

######################################################################
###function to create the Figure 1 plots - ISARs for A, S and G######
######################################################################

ANC_plot <- function(dat, A, S, G, ti = "test", yl = ""){
  #ti = main title, yl = y-axis name
  
  #get regression models using ANC function
  fits <- ANC(dat = dat, A = A, S = S, G = G, plot = TRUE)
  
  #range of y-values for ylim
  ran <- range(c(log10(G$s + 0.1), log10(S$s + 0.1), log10(A$s + 0.1)))
  
  #build blank plot and add points for A, S and G
  par(adj = 0.5)
  plot(log10(s) ~ log10(a), data = A, 
       type='n', xlab="Log area", ylab= yl,  cex.axis=2, cex.lab=2,
       ylim = ran)
  par(adj=0.05)#allows title to be on left margin
  title(main=ti, col.main="black",cex.main=2.05, line=-1.5)
  points(log10(G$a), log10(G$s + 0.1), pch=20, col="red")
  points(log10(S$a), log10(S$s + 0.1), pch=20, col="blue")
  points(log10(A$a), log10(A$s + 0.1), pch=4)  
  
  #add regression lines
  abline(fits[[1]], lty = 1, col = "red", lwd = 2)
  abline(fits[[2]], lty = 1, col = "blue", lwd = 2)
  abline(fits[[3]],lty = 1, lwd = 2)
}


#####################################################################
##get non-linear power statistics (z, p-value (of z), z CIs, c, R2)
####################################################################

non_lin <- function(dat){
  np <- sar_power(dat)
  cc <- np$par[1]
  z <- np$par[2]
  p <- as.vector(np$sigConf[2, "Pr(>|t|)"])
  r2 <- np$R2
  ci <- np$sigConf[2, c("2.5%", "97.5%")]
  #if nls confint has errored it returns nothing and the ncols of the sigconf
  #matrix is six, so if this happens, we cant get the CIs from nls method
  if (ncol(np$sigConf) > 6){
  #compare sars par estimates with nls estimates
  if (round(np$sigConf[2, "Estimate"], 2) !=
      round(np$sigConf[2, "nls.Est."], 2)) stop ("nls par est. doesn't match")
  ci2 <- np$sigConf[2, c("nls.2.5%", "nls.97.5%")]
  } else {
    ci2 <- c(NA, NA)
  }
  res <- round(c(z, "P" = p, ci, ci2, cc, "R2" = r2),2)
  return(res)
}

#############################################################################
####MMI analyses - model comparison and 1st derivative of mmi curve########
############################################################################

##uses two methods for calculating mmi curve 1st deriv:
#1) numerical approximation used in the paper
#2) caclulates 1st derivative of each constituent model curve using the exact
#derivative functions now in sars. Then weights these.

#curve_cross argument for if function is used with deriv_cross function. If TRUE,
#this checks how long the sequences of area values created is - if over 10,000 it makes 
#a sequence of length 10,000 instead, as if not deriv_cross takes too long to run
mmi_fun <- function(dat, GS = "partial", NT = "none", HT = "none", HC = NULL,
                    curve_cross = FALSE, verb = FALSE, display = FALSE){
  
 mod <- c("power", "loga", "negexpo", "monod", "ratio", "heleg", "weibull3")
  
  s <- sar_average(data = dat, obj = mod, normaTest = NT, homoTest = HT,
                   homoCor = HC, crit = "AICc", grid_start = GS, verb = verb,
                   display = display)
  
  ###get delta AICc values
  AIC_vals <- vapply(s$details$fits, function(x){
    x$AICc
  }, 
  FUN.VALUE = numeric(1))
  
  min_AIC <- min(AIC_vals)
  del_AIC <- AIC_vals - min_AIC
  
  #check names match (wont work if models removed from checks)
  identical(names(s$details$fits), mod)
  
  ##create sequence of area values for smooth curve
  yt <- seq(min(dat$a), floor(max(dat$a)), 0.25)
  if (curve_cross){
    if (length(yt) > 10000){
      rt <- floor(max(dat$a)) - min(dat$a)
      rtt <- rt / 10000
      yt <- seq(min(dat$a), floor(max(dat$a)), rtt)
    }
  }
  
  #get weights and model names from summary table
  wei_tab <- summary(s)$Model_table[,1:2]
  
  #match order of names in table to order in fits
  mm <- match(names(s$details$fits), wei_tab$Model)
  
  #get weights in correct order (i.e. order of the fits in s)
  wei_ord <- wei_tab$Weight[mm]
  
  ##iterate across model fits in s to get calculated values using fitted par values and new
  #sequence of area values
  #returns a matrix
  fitt <- vapply(s$details$fits, function(x){
    pPars = x$par
    x$model$mod.fun(yt, pPars)
  }, FUN.VALUE = numeric(length(yt)))

  ##multiply each model's calculated values with its weight
  fitt2 <- matrix(NA, nrow = nrow(fitt), ncol = ncol(fitt))
  for (i in 1:ncol(fitt2)){
    fitt2[,i] <- fitt[,i] * wei_ord[i]
  }
  
  #sum each row to get final mmi values
  nyt <- rowSums(fitt2)
  
  dXt <- rowMeans(embed(yt,2)) # centers the X values for plotting
  dYt <- diff(nyt)/diff(yt) # the derivative
  
  da_difft <- as.data.frame(cbind(dXt,dYt))
  
  ##Repeat using weighted 1st derivatives directly from sars package
  
  ##iterate across model fits in s to get 1st derivative values 
  #returns a matrix
  der <- vapply(s$details$fits, function(x){
    pPars = x$par
    x$model$d1.fun(yt, pPars)
  }, FUN.VALUE = numeric(length(yt)))
  
  #very occasionally, using this method can return a NaN for
  #a model. If this happens, just replace the model der values
  #with those from the approximation method
  if (any(apply(der, 2, anyNA))){
    cat("New Deriv method returned NA for a model\n")
    wNA <- which(apply(der, 2, anyNA))
    der[, wNA] <- fitt[, wNA]
    if (colnames(der)[wNA] != colnames(fitt)[wNA]) stop("abnormal")
  }
  
  ##multiply each model's deriv values with its weight
  fitt3 <- matrix(NA, nrow = nrow(der), ncol = ncol(der))
  for (i in 1:ncol(fitt3)){
    fitt3[,i] <- der[,i] * wei_ord[i]
  }
  
  #get sum of weighted deriv values
  nyt2 <- rowSums(fitt3)
  
  da_difft2 <- as.data.frame(cbind(yt,nyt2))
  colnames(da_difft2) <- c("dXt", "dYt")
  
  ll <- list(del_AIC, da_difft, da_difft2)
  
  return(ll)
}

#cut_off truncates the plot at the RHS if needed, Can be be a number, which
#then cuts off all dXt values below this, or it can be "fifty" which works
#out the half way point along dXt (i.e 50%) and uses this
#curve_cross = whether or not to run the deriv_cross function, if true
#this is run and result returned without plotting curves. If generalist returned
#it means this is always above specialist, and vice versa (cross means curves cross)
#der_method means to use our new method (3) which takes derivatives directly from
#functions in sars, or use our original method (2)
deriv_plot <- function(S, G, GS = "partial", ti = "test", cut_off = NULL,
                       NT = "none", HT = "none", HC = NULL,
                       curve_cross = FALSE, der_method = 3){
  
  #get derivative for S and G using mmi_fun
  #3rd element is the deriv (using new method)
  Sder <- mmi_fun(S[,c("a", "s")], GS, NT = NT, HT = HT, 
                  HC = HC, curve_cross = curve_cross)[[der_method]]
  Gder <- mmi_fun(G[,c("a", "s")], GS, NT = NT, HT = HT,
                  HC = HC, curve_cross = curve_cross)[[der_method]]
  
  ##create df for ggplot2
  #first add factor code to each
  Sder$Type <- "Specialists"
  Gder$Type <- "Generalists"
  #combine both into one
  Ader <- as.data.frame(rbind(Sder, Gder))
  
  ##if curve_cross argument = TRUE, then do this calculation and
  #return results
  if (curve_cross){
    cc <- deriv_cross(Ader)
    return(cc)
  }
  
  #truncate graph if cut_off is supplied
  if (!is.null(cut_off)){
    if (cut_off == "fifty"){
      cut_off <- max(Ader$dXt) / 2
      Ader <- filter(Ader, dXt < cut_off)
    } else {
    Ader <- filter(Ader, dXt < cut_off)
    }
  }
  
  #build plot
  g <- ggplot(data = Ader) + geom_line(aes(x=dXt, y=dYt, colour = Type), 
                                       size = 1) +
    scale_colour_manual("",breaks = c("Generalists", "Specialists"), 
                        values = c("red","black")) + theme_bw() +
    ggtitle(ti) + guides(colour=FALSE) + 
    xlab("Area(Ha)") + ylab("dY") + 
    theme(axis.title.y=element_text(size=18)) +
    theme(axis.title.x=element_text(size=18)) + 
    theme(axis.text=element_text(size=16)) + 
    theme(plot.title=element_text(size=20,hjust = 0,vjust=1.5))
  
  return(g)
}

#function to work out whether mmi deriv lines cross or
#not. Works inside the deriv_plot by setting the argument curve_cross
#to TRUE.
deriv_cross <- function(Ader){
  
  #Ader is the dataframe inside deriv_plot with the dxT and dYt values for specialists
  #and generalists. The dxt values should be identical in the two subsets. So here take
  #just specialists and work out how many dxt values there are
  f <- filter(Ader, Type == "Specialists")
  res <- vector(length = nrow(f))
  
  #iterate across each dxt value, then take out the rows from Ader with this 
  #dxt value (i.e. the one specialist and one generalist row with his dxt value)
  for (i in 1:nrow(f)){
    dum <- as.vector(unlist(f[i,1]))
    f2 <- filter(Ader, dXt == dum)
    if (nrow(f2) != 2) stop("comfortably numb")#should only be two occurences of each dxt value
    #work out if, for a given dxt value, whether the generalist or specialist dyt value is larger
    #i.e. which of the curves is "higher"
    #can't use which.max to check if equal as this is, if they are equal, only returns the
    #first instance
    if (f2$dYt[1] == f2$dYt[2]){
      res[i] = "Equal"
    } else {
      w <- which.max(f2$dYt)
      res[i] = as.vector(unlist(f2[w,][3]))
    }
  }
  
  #check if all dyt values are higher for generalists or specialsits (i.e. curves never cross)
  all_G <- all(stringr::str_detect(res, "Generalist"))
  all_S <- all(stringr::str_detect(res, "Specialist"))
  
  #record results, if neither all_G or allS = TRUE, the curves must cross at some point
  if (all_G){
    result <- "Generalist"
  } else if (all_S){
    result <- "Specialist"
  } else{
    result <- "Cross"
  }
  
  return(result)
  
}

#########################################
##plot SG Ratio plots
###############################################

SG_ratio <- function(S, G, ti = "(a)"){
  
  #calculate ratio
  rat <- S$s / G$s
  
  par(adj=0.5)
  plot(S$a,rat, pch=20, cex=1.6, col=col[3], xlab="Area (Ha)", 
       ylab="S:G Ratio", cex.lab=1.88, cex.axis = 1.88)
  par(adj=0)#allows title to be on left margin
  title(main=ti, col.main="black",cex.main=2.05, line=1)
  
}