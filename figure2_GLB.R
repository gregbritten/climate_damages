# Figure 2
rm(list=ls())

library(R.matlab)

files <- list.files(pattern="equals1")

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
rho <- c(0.04255,0.015) # pure rate of time preference (1.5% in DICE-2016)
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

# Discount factors
df.2 <- 1/((1+rho[s])^seq(0,n-1,1))
df.2 <- matrix(rep(df.2,nn+1),ncol=nn+1,nrow=n)

x <- data.frame()
for (j in c("SSP370","SSP460")) {
  t <- d[d$rcp==j,]
  t <- t[,!(names(t) %in% c("year","rcp","r","g","b"))]
  
  # Isoelastic utility
  x.Weitzman <- colSums((((gdp.ts*(1/(1+((t)/20.46)^2 + ((t)/6.081)^6)))^(1-eta[s]))/(1-eta[s])) * df.2 * pop.ts)
  x <- rbind(x, x.Weitzman)
}
x <- t(x)
colnames(x) <- files

## Scaling factors for plotting
UofC   <- ((gdp.ts[1,1])^(1-eta[s]))/(1-eta[s]) # Utility of undisturbed per capita consumption today (in thousands of dollars)
scalar <- UofC-(((gdp.ts[1,1]-0.001)^(1-eta[s]))/(1-eta[s])) # Utility value of the marginal dollar
UofC   <- sum((((gdp.ts[,1])^(1-eta[s]))/(1-eta[s]))*pop.ts[,1]*df.2[,1]) # NPV of utility of undisturbed aggregate consumption



dam_370 <- (UofC - x[2:nrow(x),grep("SSP370",colnames(x))])/scalar
dam_460 <- (UofC - x[2:nrow(x),grep("SSP460",colnames(x))])/scalar

pdf('damage_distributions_11_19_2021.pdf',height=4.5,width=5)
par(mfrow=c(1,1))
hist(dam_370/1E12,breaks=500,xlim=c(0,400),ylim=c(0,0.018),col=rgb(1,0,0,0.5),xlab='',ylab='',main='',freq=FALSE)
hist(dam_460/1E12,breaks=500,add=TRUE,col=rgb(0,0,1,0.5),freq=FALSE)
  legend("topright", legend=c("SSP370","SSP460"), col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pt.cex=2, pch=15,bty='n',cex=1.2)
  mtext(side=1,'Trillion USD',line=2.5)
  mtext(side=2,'Density',line=2.5)
dev.off()  

write.csv(file='damages_11_19_2021.csv',data.frame(SSP370=dam_370,SSP460=dam_460),row.names=FALSE)


