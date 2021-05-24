
##############################################################
########SG Analyses###################################
#########################################################

source("Source_code.R")

#read data (list named where each element is a dataset)
#note, does not include dataset 21 (dos Anjos 2004) as this was
#provided by the source paper author for revision
ldf <- readRDS(file = "datasets.rds")

##create lists and matrices to hold results

R1 <- matrix(nrow =length(ldf), ncol = 20)
colnames(R1) <- c("zA", "pA", "zS", "pS", "zG", "pG", "F", "Pf", 
                  "D_1", "Pd_1", "varD_1", "D_2", "Pd_2", "varD_2",
                  "D_4", "Pd_4", "varD_4", "D_5", "Pd_5", "varD_5")

R2 <- matrix(ncol = 24, nrow = length(ldf))
colnames(R2) <- c("zA", "pA", "CIL_A", "CIU_A", "nlsCIL_A", "nlsCIU_A", "cA", "R2A",
                  "zG", "pG", "CIL_G", "CIU_G", "nlsCIL_G", "nlsCIU_G", "cG", "R2G",
                  "zS", "pS", "CIL_S", "CIU_S", "nlsCIL_S", "nlsCIU_S", "cS", "R2S")

R3G <- matrix(ncol = 7, nrow = length(ldf))
colnames(R3G) <- c("power", "loga", "negexpo", "monod", "ratio", "heleg", "weibull3")

R3S <- matrix(ncol = 7, nrow = length(ldf))
colnames(R3S) <- c("power", "loga", "negexpo", "monod", "ratio", "heleg", "weibull3")

##Iterate across all datasets to get main results

for (i in 1:length(ldf)){

TI <- paste0("(", letters[i], ")")

##load in a dataset

dat <- ldf[[i]]

##First create the three subset datasets (all, sist and gist) from
#each file.

A <- filter(dat, factor == "a")
S <- filter(dat, factor == "s")
G <- filter(dat, factor == "g")

#some quick sense checks
t1 <- all(A$a == S$a)
t2 <- all(A$a == G$a)
t3 <- all((G$s + S$s) == A$s)
if (!all(c(t1, t2, t3))) stop("sense checks failed")

##Run ANCOVA and Cohen's D

R1[i,] <- ANC(dat = dat, A = A, S = S, G = G, plot = FALSE)

##Non-linear Power results

R2[i, c("zA", "pA", "CIL_A", "CIU_A", 
        "nlsCIL_A", "nlsCIU_A", "cA", "R2A")] <- non_lin(A[,1:2])

R2[i, c("zG", "pG", "CIL_G", "CIU_G", 
        "nlsCIL_G", "nlsCIU_G", "cG", "R2G")] <- non_lin(G[,1:2])

R2[i, c("zS", "pS", "CIL_S", "CIU_S", 
        "nlsCIL_S", "nlsCIU_S", "cS", "R2S")] <- non_lin(S[,1:2])

## MMI fitting and derivative plots

xxG <- mmi_fun(G[,1:2], GS = "partial", verb = T, display = FALSE)
if (!identical(names(xxG[[1]]), colnames(R3G))){
  stop("Hey ohhh")
}

R3G[i,] <- xxG[[1]]


xxS <- mmi_fun(S[,1:2], GS = "partial", verb = T, display = FALSE)
if (!identical(names(xxS[[1]]), colnames(R3S))){
  stop("Hey ohhh")
}

R3S[i,] <- xxS[[1]]

}#eo main for


write.csv(R1, file = "ANCOVA_results.csv")
write.csv(R2, file = "Non_lin_results.csv")
write.csv(round(R3G, 2), file = "MMI_GEN_results.csv")
write.csv(round(R3S, 2), file = "MMI_SPE_results.csv")

############################################################################
#########POST HOC ANALYSES########################################
######################################################################

#turn R1 and R2 into dataframes
#only first 20 rows as these are the main datasets we want
#for overall effect size and boxplots
R1a <- as.data.frame(R1[1:20,])
R2a <- as.data.frame(R2[1:20,])

##Get DPLUS

#select which d value (and variance) you want to use, from our three choices
vv <- select(R1a, "d" = D_1, "var" = varD_1)
dP <- dPlus(d = vv$d, var = vv$var)

# calculate d+ for the other outcome-covariate methods
vv2 <- select(R1a, "d" = D_2, "var" = varD_2)
dP2 <- dPlus(d = vv2$d, var = vv2$var)
vv4 <- select(R1a, "d" = D_4, "var" = varD_4)
dP4 <- dPlus(d = vv4$d, var = vv4$var)
vv5 <- select(R1a, "d" = D_5, "var" = varD_5)
dP5 <- dPlus(d = vv5$d, var = vv5$var)
c(dP, dP2, dP4, dP5)

##t-test for log-log and non-linear z values
#doesn't include 2 non-island datasets, so used version a of R1/R2
t.test(R2a$zS, R2a$zG)
t.test(R1a$zS, R1a$zG)

####FIGURE 1
nams2 <- c(1, 3, 5, 6, 8, 12)

ldf2 <- ldf[nams2]

jpeg(file = "Figure1.jpeg", width = 35, height = 35, units = "cm", res = 300)

par(mfrow = c(2,3))

for (i in 1:length(ldf2)){
  
  ti2 <- paste0("(", letters[i], ")")
  dat2 <- ldf2[[i]]
  A2 <- filter(dat2, factor == "a")
  S2 <- filter(dat2, factor == "s")
  G2 <- filter(dat2, factor == "g")
  par(mar=c(5, 4.5, 4, 2)+0.1)
  if (i %in% c(1, 4)){
    ANC_plot(dat2, A = A2, S = S2, G = G2, ti = ti2, yl = "Log species richness (+0.1)")
  } else {
    ANC_plot(dat2, A = A2, S = S2, G = G2, ti = ti2, yl = "")
  }
  
}

dev.off()

####FIGURE 2

jpeg(file = "Figure2.jpeg", width = 25, height = 13, units = "cm", res = 300)

par(mfrow=c(1,2))

par(adj=0.5)
boxplot(R1a$zA, R1a$zG, R1a$zS, xlab="Species subset", ylab="z value", 
        cex.lab=1.3, names=c("All","Gen","Spec"), cex.axis=1.2, medcol="red", 
        outline=FALSE, col = "white")
par(adj=0)#allows title to be on left margin
title(main="(a)", col.main="black",cex.main=1.4, line=2)


par(adj=0.5)
boxplot(R2a$zA, R2a$zG, R2a$zS, xlab="Species subset", 
        ylab="z value", cex.lab=1.3, names=c("All","Gen","Spec"), 
        cex.axis=1.2,medcol="red",outline=FALSE, col = "white")

par(adj=0)#allows title to be on left margin
title(main="(b)", col.main="black",cex.main=1.4, line=2)

dev.off()

####FIGURE 3

nams3 <- c(1, 5, 12)

ldf3 <- ldf[nams3]

#co <- c() #cut_off values

mb <- vector("list", length = length(ldf3))

for (i in 1:length(ldf3)){
  
  ti3 <- paste0("(", letters[i], ")")
  dat3 <- ldf3[[i]]
  S3 <- filter(dat3, factor == "s")
  G3 <- filter(dat3, factor == "g")
  ##For anciaes and marini, have focused on the first 100ha of the gradient
  if (i == 1){
    mb[[i]] <- deriv_plot(S = S3, G = G3, GS = "partial", ti = ti3, 
                          cut_off = 100, der_method = 3)
    #for the other two, focus is on the first 50% of area gradient
  } else{
    mb[[i]] <- deriv_plot(S = S3, G = G3, GS = "partial", ti = ti3, 
                          cut_off = "fifty", der_method = 3)
  }
}

jpeg(file = "Figure3.jpeg", width = 35, height = 15, units = "cm", res = 300)

grid.arrange(mb[[1]], mb[[2]], mb[[3]], nrow = 1)

dev.off()

#supp info version: split into three (7, 7, 6)
#and also run analysis of seeing whether the deriv curves cross

nams3_SI <- vector("list", length = 3)
nams3_SI[[1]] <- 1:9
nams3_SI[[2]] <- 10:18
nams3_SI[[3]] <- 19:22


pn3_SI <- c("FigS1a.jpeg", "FigS1b.jpeg","FigS1c.jpeg")
hi3 <- c(35, 35, 25)
rw3 <- c(3,3,2)

k <- 1
ccR <- vector(length = length(ldf))#vector to store curve cross information

for (j in 1:3){
  
  print(j)
  
ldf3 <- ldf[nams3_SI[[j]]]

let3 <- letters[nams3_SI[[j]]]

mb <- vector("list", length = length(ldf3))

for (i in 1:length(ldf3)){
  
  print(i)
  
  ti3 <- paste0("(", let3[i], ")")
  dat3 <- ldf3[[i]]
  S3 <- filter(dat3, factor == "s")
  G3 <- filter(dat3, factor == "g")
  mb[[i]] <- deriv_plot(S = S3, G = G3, GS = "partial", 
                        ti = ti3, cut_off = "fifty",
                        curve_cross = FALSE, der_method = 3)
  ccR[k] <- deriv_plot(S = S3, G = G3, GS = "partial", 
                       curve_cross = TRUE, der_method = 3)
  k <- k + 1
  
}#eo i

jpeg(file = pn3_SI[j], width = 35, height = hi3[j], units = "cm", res = 300)

grid.arrange(grobs = mb, nrow = rw3[j], ncol = 3)

dev.off()

}#eo j

ccR #deriv curve cross results
table(ccR[1:20])


#####FIGURE 4

nams4 <- c(1, 12, 14, 17)

ldf4 <- ldf[nams4]

jpeg(file = "Figure4.jpeg", width = 35, height = 35, units = "cm", res = 300)

par(mfrow = c(2,2))

for (i in 1:length(ldf4)){
  
ti4 <- paste0("(", letters[i], ")")
dat4 <- ldf4[[i]]
S4 <- filter(dat4, factor == "s")
G4 <- filter(dat4, factor == "g")
par(mar=c(5, 4.5, 4, 2)+0.1)
SG_ratio(S = S4, G = G4, ti = ti4)

}

dev.off()


#supp info version: split into three (7, 7, 6)

nams4_SI <- vector("list", length = 3)
nams4_SI[[1]] <- 1:9
nams4_SI[[2]] <- 10:18
nams4_SI[[3]] <- 19:22


pn4_SI <- c("FigS2a.jpeg", "FigS2b.jpeg","FigS2c.jpeg")
hi4 <- c(35, 35, 25)
rw4 <- c(3,3,2)


for (j in 1:3){
  
ldf4 <- ldf[nams4_SI[[j]]]

jpeg(file = pn4_SI[j], width = 35, height = hi4[j], units = "cm", res = 300)

par(mfrow = c(rw4[j],3))

let4 <- letters[nams4_SI[[j]]]

for (i in 1:length(ldf4)){
  
  ti4 <- paste0("(", let4[i], ")")
  dat4 <- ldf4[[i]]
  S4 <- filter(dat4, factor == "s")
  G4 <- filter(dat4, factor == "g")
  par(mar=c(5, 4.5, 4, 2)+0.1)
  SG_ratio(S = S4, G = G4, ti = ti4)
  
}#eo i

dev.off()

}#eo j
