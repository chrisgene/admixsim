# 07/15/18
# simulate admixed population allele frequencies from marginal betas

# original function by Noah Zaitlen, modified by Gillian Belbin and Christopher Gignoux



fhudson = function(p1,p2,n1,n2) {
# calculate fst using hudson's formula, from Bhatia et al Gen Res 2013
numerator = (p1-p2) ^2  - (p1*(1-p1)/(2*n1-1)) - (p2*(1-p2)/(2*n2-1))
denominator = p1*(1-p2) + p2*(1-p1)
# across all values, return ratio of averages
return(mean(numerator) / mean(denominator))
}



N = 500 #number of individuals
M = 1000 #number of SNPs



# AFR-AMR, AFR-EUR , AMR-EUR Fsts from Moreno-Estrada plos gen
fsts = c(0.27 , 0.16 , 0.16 )
# use lsbls to calculate drift from ancestral pop, a la Shriver et al 2004, AFR/AMR/EUR
# note: used modified formula here (2xlsbl) to get beta distributions to be parameterized correctly

plsbls = c( (fsts[1] + fsts[2] - fsts[3]), (fsts[1] + fsts[3] - fsts[2]) , (fsts[2] + fsts[3] - fsts[1]) )
mlsbls = c( (fsts[1] + fsts[2] - fsts[3]), (fsts[1] + fsts[3] - fsts[2]) , (fsts[2] + fsts[3] - fsts[1]) )


# 3-way population admixture using values derived from real-world data...

# thetas parameterized as AFR, AMR, EUR, calculated from GALA II (Galanter et al 2012)
prthetas = c(0.20763 , 0.10332 , 0.6890)
prsds = c( 0.1067561 , 0.04132107 , 0.1051037)

mxthetas = c(0.03807 , 0.5781 , 0.38383)
mxsds = c(0.1638042 , 0.1736818 , 0.02851813)

# simple truncated normal distribution for now...
ptheta1 = rnorm(N,prthetas[1],prsds[1])
ptheta1[ptheta1 < 0] = 0
ptheta2 = rnorm(N,prthetas[2],prsds[2])
ptheta2[ptheta2 < 0] = 0
ptheta2[ptheta2+ptheta1>1] = 1-ptheta1[ptheta2+ptheta1>1]
ptheta3 = 1-(ptheta1+ptheta2)
ptheta=cbind(ptheta1,ptheta2,ptheta3)

# simple truncated normal distribution for now...
mtheta1 = rnorm(N,mxthetas[1],mxsds[1])
mtheta1[mtheta1 < 0] = 0
mtheta2 = rnorm(N,mxthetas[2],mxsds[2])
mtheta2[mtheta2 < 0] = 0
mtheta2[mtheta2+mtheta1>1] = 1-mtheta1[mtheta2+mtheta1>1]
mtheta3 = 1-(mtheta1+mtheta2)
mtheta=cbind(mtheta1,mtheta2,mtheta3)



ancFreqs = runif(M,0.05,0.95)

pfreqs1=pfreqs2=pfreqs3=0
for (i in 1:M) {
    pfreqs1[i] = rbeta(1,ancFreqs[i]*(1-plsbls[1])/plsbls[1],(1-ancFreqs[i])*(1-plsbls[1])/plsbls[1])
    pfreqs2[i] = rbeta(1,ancFreqs[i]*(1-plsbls[2])/plsbls[2],(1-ancFreqs[i])*(1-plsbls[2])/plsbls[2])
    pfreqs3[i] = rbeta(1,ancFreqs[i]*(1-plsbls[3])/plsbls[3],(1-ancFreqs[i])*(1-plsbls[3])/plsbls[3])
}

mfreqs1=mfreqs2=mfreqs3=0
for (i in 1:M) {
    mfreqs1[i] = rbeta(1,ancFreqs[i]*(1-mlsbls[1])/mlsbls[1],(1-ancFreqs[i])*(1-mlsbls[1])/mlsbls[1])
    mfreqs2[i] = rbeta(1,ancFreqs[i]*(1-mlsbls[2])/mlsbls[2],(1-ancFreqs[i])*(1-mlsbls[2])/mlsbls[2])
    mfreqs3[i] = rbeta(1,ancFreqs[i]*(1-mlsbls[3])/mlsbls[3],(1-ancFreqs[i])*(1-mlsbls[3])/mlsbls[3])
}



## decided to do the gamma matrix as 2x matrices of alleles instead of 1x matrix of genotypes

pgamma1 = matrix(0,nrow=N,ncol=M)
for(i in 1:N) {
       for(j in 1:M) {    
    x <- rmultinom(1, size=1, prob=ptheta[i,])
    pgamma1[i,j]=row(x)[which(x!=0,arr.ind = T)]
    }
}
pgamma2 = matrix(0,nrow=N,ncol=M)
for(i in 1:N) {
    for(j in 1:M) {
    x <- rmultinom(1, size=1, prob=ptheta[i,])
    pgamma2[i,j]=row(x)[which(x!=0,arr.ind = T)]
    }
}


pgenos = matrix(0,nrow=N,ncol=M)
for(i in 1:N) {
    for(j in 1:M) {
        if (pgamma1[i,j]==1 & pgamma2[i,j]==1) {
            pgenos[i,j] = rbinom(1,2,pfreqs1[j])
        }
        if (pgamma1[i,j]==2 & pgamma2[i,j]==2) {
            pgenos[i,j] = rbinom(1,2,pfreqs2[j])
        }
         if (pgamma1[i,j]==3 & pgamma2[i,j]==3) {
            pgenos[i,j] = rbinom(1,2,pfreqs3[j])
        }
        if ((pgamma1[i,j]== 1 & pgamma2[i,j]==2)  | (pgamma1[i,j]== 2 & pgamma2[i,j]== 1)) {
            pgenos[i,j] = rbinom(1,1,pfreqs1[j])+rbinom(1,1,pfreqs2[j])
        }
        if ((pgamma1[i,j]== 1 & pgamma2[i,j]== 3)  | (pgamma1[i,j]== 3 & pgamma2[i,j]== 1)) {
            pgenos[i,j] = rbinom(1,1,pfreqs1[j])+rbinom(1,1,pfreqs3[j])
        }
         if ((pgamma1[i,j]== 2 & pgamma2[i,j]== 3)  | (pgamma1[i,j]== 3 & pgamma2[i,j]== 2)) {
            pgenos[i,j] = rbinom(1,1,pfreqs2[j])+rbinom(1,1,pfreqs3[j])
        }
    }        
}


mgamma1 = matrix(0,nrow=N,ncol=M)
for(i in 1:N) {
       for(j in 1:M) {    
    x <- rmultinom(1, size=1, prob=mtheta[i,])
    mgamma1[i,j]=row(x)[which(x!=0,arr.ind = T)]
    }
}
mgamma2 = matrix(0,nrow=N,ncol=M)
for(i in 1:N) {
    for(j in 1:M) {
    x <- rmultinom(1, size=1, prob=mtheta[i,])
    mgamma2[i,j]=row(x)[which(x!=0,arr.ind = T)]
    }
}


mgenos = matrix(0,nrow=N,ncol=M)
for(i in 1:N) {
    for(j in 1:M) {
        if (mgamma1[i,j]==1 & mgamma2[i,j]==1) {
            mgenos[i,j] = rbinom(1,2,mfreqs1[j])
        }
        if (mgamma1[i,j]==2 & mgamma2[i,j]==2) {
            mgenos[i,j] = rbinom(1,2,mfreqs2[j])
        }
         if (mgamma1[i,j]==3 & mgamma2[i,j]==3) {
            mgenos[i,j] = rbinom(1,2,mfreqs3[j])
        }
        if ((mgamma1[i,j]== 1 & mgamma2[i,j]==2)  | (mgamma1[i,j]== 2 & mgamma2[i,j]== 1)) {
            mgenos[i,j] = rbinom(1,1,mfreqs1[j])+rbinom(1,1,mfreqs2[j])
        }
        if ((mgamma1[i,j]== 1 & mgamma2[i,j]== 3)  | (mgamma1[i,j]== 3 & mgamma2[i,j]== 1)) {
            mgenos[i,j] = rbinom(1,1,mfreqs1[j])+rbinom(1,1,mfreqs3[j])
        }
         if ((mgamma1[i,j]== 2 & mgamma2[i,j]== 3)  | (mgamma1[i,j]== 3 & mgamma2[i,j]== 2)) {
            mgenos[i,j] = rbinom(1,1,mfreqs2[j])+rbinom(1,1,mfreqs3[j])
        }
    }        
}


pmafs = apply(pgenos,2,mean)/2
mmafs = apply(mgenos,2,mean)/2
mafs = (pmafs+mmafs)/2
FstEst = mean((pmafs-mmafs)^2)/mean(2*mafs*(1-mafs)) ## Noah's fstEstimator

h_fst_est <- fhudson(pmafs,mmafs,N,N)

m <- cbind(mmafs,"Mexican")
p <- cbind(pmafs,"Puerto Rican")
x <- as.data.frame(rbind(m,p))
colnames(x) <- c("MAF","Population")
ggplot(data=x) + geom_density(aes(x=as.numeric(MAF),fill=Population),alpha=0.2)
ggplot(data=as.data.frame(cbind(mmafs,pmafs))) + geom_point(aes(x=pmafs,y=mmafs)) + theme_classic() + xlab("Puerto Rican (Allele Frequency)") + ylab("Mexican (Allele Frequency)")
hist(pmafs-mmafs)

## Compare allele freq differences between upper and lower quartiles of dominant ancestral component WITHIN populations

pquartiles <- as.data.frame(as.vector(quantile(ptheta[,3], probs = seq(0, 1, 0.25))))

ptheta3 <- as.data.frame(ptheta3)
upper_p_genos = matrix(0,nrow=N,ncol=M)
for(i in 1:N) {
    if (ptheta3[i,] >= pquartiles[4,]){
        upper_p_genos[i,]=pgenos[i,]
    }
    else {
        upper_p_genos[i,]=NA
    }
}
lower_p_genos = matrix(0,nrow=N,ncol=M)
for(i in 1:N) {
    if (ptheta3[i,] <= pquartiles[2,]){
        lower_p_genos[i,]=pgenos[i,]
    }
    else {
        lower_p_genos[i,]=NA
    }
}
lpmafs = apply(lower_p_genos,2,mean,na.rm=T)/2
upmafs = apply(upper_p_genos,2,mean,na.rm=T)/2

withinpr_fst_est <- fhudson(lpmafs,upmafs,length(ptheta3[which(ptheta3 >= pquartiles[4,]),]),length(ptheta3[which(ptheta3 >= pquartiles[4,]),]))

mquartiles <- as.data.frame(as.vector(quantile(mtheta[,2], probs = seq(0, 1, 0.25))))
mtheta2 <- as.data.frame(mtheta2)
upper_m_genos = matrix(0,nrow=N,ncol=M)
for(i in 1:N) {
    if (mtheta2[i,] >= mquartiles[4,]){
        upper_m_genos[i,]=mgenos[i,]
    }
    else {
        upper_m_genos[i,]=NA
    }
}
lower_m_genos = matrix(0,nrow=N,ncol=M)
for(i in 1:N) {
    if (mtheta2[i,] <= mquartiles[2,]){
        lower_m_genos[i,]=mgenos[i,]
    }
    else {
        lower_m_genos[i,]=NA
    }
}

lmmafs = apply(lower_m_genos,2,mean,na.rm=T)/2
ummafs = apply(upper_m_genos,2,mean,na.rm=T)/2
withinmex_fst_est <- fhudson(lmmafs,ummafs,length(mtheta2[which(mtheta2 >= mquartiles[4,]),]),length(mtheta2[which(mtheta3 >= mquartiles[4,]),]))

