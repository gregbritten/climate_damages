# Figure 2
rm(list=ls())

library(R.matlab)

files <- list.files(pattern="equals1")

# GLB: two different discounting rates
rho1 <- 0.04255 # pure rate of time preference (1.5% in DICE-2016)
rho2 <- 0.02955 

#GLB: original temperature timee series started at 1 not 0
ds_ssp370      <- data.frame(readMat(files[1])$SSP370)
ds_ssp370$year <- 2020:(2020+287)
ds_ssp370$rcp  <- "SSP370"
ds_ssp460      <- data.frame(readMat(files[2])$SSP460)
ds_ssp460$year <- 2020:(2020+287)
ds_ssp460$rcp  <- "SSP460"

d <- rbind(ds_ssp370,ds_ssp460)
d <- d[,c(10001,1:10000,10002)]
d[,c(2:10001)] <- d[,c(2:10001)] + 1  #adding one degree here

# Ensemble size
nn = ncol(d) - 2

# Length of time series
n=length(unique(d$year))

# Economic assumptions 
gdp <- 80*10^12         # initial aggregate consumption (80 trillion US dollars)
g   <- 0.019            # annual growth rate of undisturbed aggregate consumption
eta <- c(0,1.45)        # marginal utility of consumption
s   <- 1                # Index for discounting assumptions. When s=1, r = 4.255 + 0 * 0.019 = 4.255.

# Demographic assumptions
pop        <- 7.5     # billion people
pop.asympt <- 11.5    # Asymptotic population, in billion
param      <- 0.134/5 # DICE-2016 "Population 2050 parameter"
pop.ts     <- rep(NA,n)
pop.ts[1]  <- pop
for (i in 2:n) { pop.ts[i] <- pop.ts[i-1]*(pop.asympt/pop.ts[i-1])^param}
pop.ts <- matrix(rep(pop.ts,nn+1),ncol=nn+1,nrow=n)
pop.ts <- pop.ts*10^9 # Population

# Undisturbed output per capita
gdp.ts <- gdp*(1+g)^seq(0,n-1,1)
gdp.ts <- matrix(rep(gdp.ts,nn+1),ncol=nn+1,nrow=n)
gdp.ts <- (gdp.ts/pop.ts)/1000 # undisturbed consumption per capita expressed in thousands of dollars.

######################################################
## GLB: DISCOUNT RATE ################################
######################################################
# Discount factors
dfrho1 <- 1/((1+rho1)^seq(0,n-1,1))
dfrho1 <- matrix(rep(dfrho1,nn+1),ncol=nn+1,nrow=n)
dfrho2 <- 1/((1+rho2)^seq(0,n-1,1))
dfrho2 <- matrix(rep(dfrho2,nn+1),ncol=nn+1,nrow=n)

xwrho1  =xwrho2   <- data.frame() #empty data frame for Weitzman function
xhs1rho1=xhs1rho2 <- data.frame() #howard and sterner D = 1.1450*T^2
xhs2rho1=xhs2rho2 <- data.frame() #howard and sterner D = 0.07438T^2
xkwrho1 =xkwrho2  <- data.frame() 

for (j in c("SSP370","SSP460")) {
  t <- d[d$rcp==j,]
  t <- t[,!(names(t) %in% c("year","rcp","r","g","b"))]
  
  ############################################################################
  ## GLB: DAMAGE FUNCTION CALCULATION FROM Weitzman ##########################
  ############################################################################
  # Isoelastic utility
  xwwrho1 <- colSums((((gdp.ts*(1/(1+((t)/20.46)^2 + ((t)/6.081)^6)))^(1-eta[s]))/(1-eta[s])) * dfrho1 * pop.ts)
  xwwrho2 <- colSums((((gdp.ts*(1/(1+((t)/20.46)^2 + ((t)/6.081)^6)))^(1-eta[s]))/(1-eta[s])) * dfrho2 * pop.ts)
  
  xxhs1rho1 <- colSums(1.1450*t^2 * dfrho1 * pop.ts)
  xxhs1rho2 <- colSums(1.1450*t^2 * dfrho2 * pop.ts)
  
  xxhs2rho1 <- colSums(0.7438*t^2 * dfrho1 * pop.ts)
  xxhs2rho2 <- colSums(0.7438*t^2 * dfrho2 * pop.ts)
  
  xxkwrho1 <- colSums(-(-0.0373*t + (-0.0018/2)*t^2) * dfrho1 * pop.ts)
  xxkwrho2 <- colSums(-(-0.0373*t + (-0.0018/2)*t^2) * dfrho2 * pop.ts)
  
  xwrho1 <- rbind(xwrho1, xwwrho1)
  xwrho2 <- rbind(xwrho2, xwwrho2)
  
  xhs1rho1 <- rbind(xhs1rho1,xxhs1rho1)
  xhs1rho2 <- rbind(xhs1rho2,xxhs1rho2)
  
  xhs2rho1 <- rbind(xhs2rho1,xxhs2rho1)
  xhs2rho2 <- rbind(xhs2rho2,xxhs2rho2)
  
  xkwrho1 <- rbind(xkwrho1,xxkwrho1)
  xkwrho2 <- rbind(xkwrho2,xxkwrho2)
}
xwrho1 <- t(xwrho1)
xwrho2 <- t(xwrho2)
colnames(xwrho1) <- files
colnames(xwrho2) <- files

xhs1rho1 <- t(xhs1rho1)
xhs1rho2 <- t(xhs1rho2)
colnames(xhs1rho1) <- files
colnames(xhs1rho2) <- files

xhs2rho1 <- t(xhs2rho1)
xhs2rho2 <- t(xhs2rho2)
colnames(xhs2rho1) <- files
colnames(xhs2rho2) <- files

xkwrho1 <- t(xkwrho1)
xkwrho2 <- t(xkwrho2)
colnames(xkwrho1) <- files
colnames(xkwrho2) <- files

## Scaling factors for plotting
UofC   <- ((gdp.ts[1,1])^(1-eta[s]))/(1-eta[s]) # Utility of undisturbed per capita consumption today (in thousands of dollars)
scalar <- UofC-(((gdp.ts[1,1]-0.001)^(1-eta[s]))/(1-eta[s])) # Utility value of the marginal dollar

UofCrho1   <- sum((((gdp.ts[,1])^(1-eta[s]))/(1-eta[s]))*pop.ts[,1]*dfrho1[,1]) # NPV of utility of undisturbed aggregate consumption
UofCrho2   <- sum((((gdp.ts[,1])^(1-eta[s]))/(1-eta[s]))*pop.ts[,1]*dfrho2[,1])

dam_370_w_rho1 <- (UofCrho1 - xwrho1[2:nrow(xwrho1),grep("SSP370",colnames(xwrho1))])/scalar
dam_460_w_rho1 <- (UofCrho1 - xwrho1[2:nrow(xwrho1),grep("SSP460",colnames(xwrho1))])/scalar
dam_370_w_rho2 <- (UofCrho2 - xwrho2[2:nrow(xwrho2),grep("SSP370",colnames(xwrho2))])/scalar
dam_460_w_rho2 <- (UofCrho2 - xwrho2[2:nrow(xwrho2),grep("SSP460",colnames(xwrho2))])/scalar

dam_370_hs1_rho1 <- (UofCrho1 - xhs1rho1[2:nrow(xhs1rho1),grep("SSP370",colnames(xhs1rho1))])/scalar
dam_460_hs1_rho1 <- (UofCrho1 - xhs1rho1[2:nrow(xhs1rho1),grep("SSP460",colnames(xhs1rho1))])/scalar
dam_370_hs1_rho2 <- (UofCrho2 - xhs1rho2[2:nrow(xhs1rho2),grep("SSP370",colnames(xhs1rho2))])/scalar
dam_460_hs1_rho2 <- (UofCrho2 - xhs1rho2[2:nrow(xhs1rho2),grep("SSP460",colnames(xhs1rho2))])/scalar

dam_370_hs2_rho1 <- (UofCrho1 - xhs2rho1[2:nrow(xhs2rho1),grep("SSP370",colnames(xhs2rho1))])/scalar
dam_460_hs2_rho1 <- (UofCrho1 - xhs2rho1[2:nrow(xhs2rho1),grep("SSP460",colnames(xhs2rho1))])/scalar
dam_370_hs2_rho2 <- (UofCrho2 - xhs2rho2[2:nrow(xhs2rho2),grep("SSP370",colnames(xhs2rho2))])/scalar
dam_460_hs2_rho2 <- (UofCrho2 - xhs2rho2[2:nrow(xhs2rho2),grep("SSP460",colnames(xhs2rho2))])/scalar

dam_370_kw_rho1 <- (UofCrho1 - xkwrho1[2:nrow(xkwrho1),grep("SSP370",colnames(xkwrho1))])/scalar
dam_460_kw_rho1 <- (UofCrho1 - xkwrho1[2:nrow(xkwrho1),grep("SSP460",colnames(xkwrho1))])/scalar
dam_370_kw_rho2 <- (UofCrho2 - xkwrho2[2:nrow(xkwrho2),grep("SSP370",colnames(xkwrho2))])/scalar
dam_460_kw_rho2 <- (UofCrho2 - xkwrho2[2:nrow(xkwrho2),grep("SSP460",colnames(xkwrho2))])/scalar

write.csv(file='damages_05_04_2022.csv',
          data.frame(SSP370_w_rho1=dam_370_w_rho1,
                     SSP460_w_rho1=dam_460_w_rho1,
                     SSP370_w_rho2=dam_370_w_rho2,
                     SSP460_w_rho2=dam_460_w_rho2,
                     SSP370_hs1_rho1=dam_370_hs1_rho1,
                     SSP460_hs1_rho1=dam_460_hs1_rho1,
                     SSP370_hs1_rho2=dam_370_hs1_rho2,
                     SSP460_hs1_rho2=dam_460_hs1_rho2,
                     SSP370_hs2_rho1=dam_370_hs2_rho1,
                     SSP460_hs2_rho1=dam_460_hs2_rho1,
                     SSP370_hs2_rho2=dam_370_hs2_rho2,
                     SSP460_hs2_rho2=dam_460_hs2_rho2,
                     SSP370_kw_rho1=dam_370_kw_rho1,
                     SSP460_kw_rho1=dam_460_kw_rho1,
                     SSP370_kw_rho2=dam_370_kw_rho2,
                     SSP460_kw_rho2=dam_460_kw_rho2),row.names=FALSE)




