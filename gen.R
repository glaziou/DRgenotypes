#########################################
# genotyping - metaanalysis
# PG, init 14/02/2017
#########################################
library(metafor)

setwd('../data/genotypes_2017/')

#----------------------------------------
# run this block only if dataset updated
# read data and drop rows with no useful data

# gen <- fread('gen.csv')
# gen <- fread('gen070317b.csv')
gen <- fread('gen140317.csv')
sel <- gen$res.mut + gen$res.nomut > 0
gen <- gen[sel]
sel <- gen$country != 'Global'
gen <- gen[sel]
gen$res <- rowSums(cbind(gen$res.mut, gen$res.nomut))
gen$sens <- rowSums(cbind(gen$sens.mut, gen$sens.nomut))
table(gen$patientGroup)
gen$patientGroup[gen$patientGroup=='New Patients'] <- 'New patients'

# aggregate ZAF
zaf <- gen[country %in% c('ZAF_GP', 'ZAF_KZN')]
zaf2 <- zaf[,.(country='ZAF',
              res.mut=sum(res.mut),
              sens.mut=sum(sens.mut),
              res.nomut=sum(res.nomut),
              sens.nomut=sum(sens.nomut),
              se=NA, sp=NA, ppv=NA, npv=NA,
              sampleSize=sum(sampleSize),
              mutations=sum(mutations),
              res=sum(res),
              sens=sum(sens)
              ), by=list(patientGroup, drug)]
zaf2 <- within(zaf2, {
  se <- res.mut/res
  sp <- sens.nomut/sens
  ppv <- res.mut/(res.mut + sens.mut)
  npv <- sens.nomut/(sens.mut + sens.nomut)
})

gen2 <- rbind(gen[country %ni% c('ZAF_GP', 'ZAF_KZN')], zaf2)
gen <- gen2

# calculate SDs and 95% CI assuming binomial errors
out <- cii(gen$res, gen$res.mut)
gen$se.sd <- out$se
gen$se.lo <- out$lower95ci
gen$se.hi <- out$upper95ci

sel <- gen$sens>0 & !is.na(gen$sens)
out2 <- cii(gen$sens[sel], gen$sens.nomut[sel])
gen$sp.sd[sel] <- out2$se
gen$sp.lo[sel] <- out2$lower95ci
gen$sp.hi[sel] <- out2$upper95ci

sel <- gen$se.sd==0
gen$se.sd[sel] <- (gen$se.hi[sel] - gen$se.lo[sel]) / 4

# save data
save(gen, file = 'gen.Rdata')
#----------------------------------------




# basic forest plots
load('gen.Rdata')

drugs <- unique(gen$drug)
gen$formerSoviet <- gen$country %in% c('AZE', 'BLR', 'UKR')


for (i in drugs){
  sel <- gen$drug==i & gen$patientGroup == 'All patients' # & gen$se>0 & gen$se<1
  fit <- rma(xi=res.mut, mi=res.nomut, measure='PLO', data=gen[sel], method='REML')
  pdf(file=paste('forestAll_',i,'.pdf', sep=''), width=10, height=8)
  forest(fit, atransf=transf.ilogit,
         slab=paste(gen$country[sel], '/', gen$patientGroup[sel]),
         xlab=paste('Sensitivity (', gen$drug[sel][1], ')', sep=''),
         subset=order(gen$country[sel]))
  dev.off()
}


# effect modifiers
mods <- data.table(drugs=drugs, testHistory=NA, testSoviet=NA, testRifRes=NA)
for (i in drugs){
  sel.all <- gen$drug==i & gen$patientGroup %in% c('All patients')
  sel.hx <- gen$drug==i & gen$patientGroup %in% c('New patients','Previously treated')

  fit.hx <- rma(xi=res.mut, mi=res.nomut, mods=~patientGroup, measure='PLO', data=gen[sel.hx], method='REML')
  fit.sov <- rma(xi=res.mut, mi=res.nomut, mods=~formerSoviet, measure='PLO', data=gen[sel.all], method='REML')
#  fit.hx <- rma(yi=se, sei=se.sd, mods=~patientGroup, measure='GEN', data=gen[sel.hx], method='REML')
#  fit.sov <- rma(yi=se, sei=se.sd, mods=~formerSoviet, measure='GEN', data=gen[sel.all], method='REML')

#  if (i %ni% c('rif_h','rif_m')) {
  if (i %ni% c('rif')) {
    sel.rif <- gen$drug==i & gen$patientGroup %in% c('Rif resistant','Rif susceptible')
    fit.rif <- rma(xi=res.mut, mi=res.nomut, mods=~factor(patientGroup=='Rif resistant'), measure='PLO', data=gen[sel.rif], method='REML')
#    fit.rif <- rma(yi=se, sei=se.sd, mods=~factor(patientGroup=='Rif resistant'), measure='GEN', data=gen[sel.rif], method='REML')
    mods$testRifRes[mods$drugs==i] <- fit.rif$QMp
  }

  mods$testHistory[mods$drugs==i] <- fit.hx$QMp
  mods$testSoviet[mods$drugs==i] <- fit.sov$QMp
  mods <- within(mods, {
    testRifRes.signif <- testRifRes < 0.05
    testSoviet.signif <- testSoviet < 0.05
    testHistory.signif <- testHistory < 0.05
  })
}
(mods)



# pooled separately in former soviet countries
# for (i in c('inh_h', 'inh2_h', 'rif_h', 'inh_m', 'inh2_m', 'rif_m')){
#   sel <- gen$drug==i & gen$patientGroup == 'All patients'
#   fit <- rma(yi=se, sei=se.sd, measure='GEN', method='REML', data=gen[sel])
#   fit.ns <- rma(yi=se, sei=se.sd, measure='GEN', method='REML', data=gen[sel & formerSoviet==FALSE])
#   fit.s <- rma(yi=se, sei=se.sd, measure='GEN', method='REML', data=gen[sel & formerSoviet==TRUE])
#   pdf(file=paste('forestAll_',i,'.pdf', sep=''), width=10, height=8)
#   forest(fit, atransf=transf.ilogit, refline=NA, addfit=FALSE,
#              slab=gen$country[sel],
#              rows=c(13:11, 6:3), ylim=c(1, 16),
#              xlab=paste('Pooled Sensitivity (', gen$drug[sel][1], ')', sep=''),
#              order=order(-gen$formerSoviet[sel]))
#   addpoly(fit.s, row= 9.5, atransf=transf.ilogit, mlab="Pooled, former Soviet Union")
#   addpoly(fit.ns, row= 1.5, atransf=transf.ilogit, mlab="Pooled, all other countries")
#   dev.off()
# }
#
# for (i in c('ofx', 'ofx2')){
#   sel <- gen$drug==i & gen$patientGroup == 'All patients' & gen$country %ni% c('PHL')
#   fit <- rma(yi=se, sei=se.sd, measure='GEN', method='REML', data=gen[sel])
#   fit.ns <- rma(yi=se, sei=se.sd, measure='GEN', method='REML', data=gen[sel & formerSoviet==FALSE])
#   fit.s <- rma(yi=se, sei=se.sd, measure='GEN', method='REML', data=gen[sel & formerSoviet==TRUE])
# #  pdf(file=paste('forestAll_',i,'.pdf', sep=''), width=10, height=8)
#   forest(fit, atransf=transf.ilogit, refline=NA, addfit=FALSE,
#          slab=gen$country[sel],
#          rows=c(12:10, 5:3), ylim=c(1, 16),
#          xlab=paste('Pooled Sensitivity (', gen$drug[sel][1], ')', sep=''),
#          order=order(-gen$formerSoviet[sel]))
#   addpoly(fit.s, row= 8.5, atransf=transf.ilogit, mlab="Pooled, former Soviet Union")
#   addpoly(fit.ns, row= 1.5, atransf=transf.ilogit, mlab="Pooled, all other countries")
# #  dev.off()
# }


# pooled separately by rif res
for (i in c('inh', 'inh2', 'inh', 'inh2', 'ofx','ofx2')){
  selrr <- gen$drug==i & gen$patientGroup == 'Rif resistant'
  selrs <- gen$drug==i & gen$patientGroup == 'Rif susceptible'
  fit.rr <- rma(xi=res.mut, mi=res.nomut, measure='PLO', method='REML', data=gen[selrr])
  fit.rs <- rma(xi=res.mut, mi=res.nomut, measure='PLO', method='REML', data=gen[selrs])

  pdf(file=paste('forestRifRes_',i,'.pdf', sep=''), width=10, height=8)
  forest(fit.rr, atransf=transf.ilogit, refline=NA, addfit=TRUE,
         slab=gen$country[selrr],
         xlab=paste('Pooled Sensitivity, Rif Resistant (', gen$drug[selrr][1], ')', sep=''),
         order=order(-gen$se[selrr]))
  dev.off()

  pdf(file=paste('forestRifSus_',i,'.pdf', sep=''), width=10, height=8)
  forest(fit.rs, atransf=transf.ilogit, refline=NA, addfit=TRUE,
         slab=gen$country[selrs],
         xlab=paste('Pooled Sensitivity, Rif Susceptible (', gen$drug[selrs][1], ')', sep=''),
         order=order(-gen$se[selrs]))
  dev.off()
}

for (i in setdiff(drugs, c('rif'))){
  selrr <- gen$drug==i & gen$patientGroup == 'Rif resistant'
  selrs <- gen$drug==i & gen$patientGroup == 'Rif susceptible'
  fit.rr <- rma(xi=res.mut, mi=res.nomut, measure='PLO', method='REML', data=gen[selrr])
  fit.rs <- rma(xi=res.mut, mi=res.nomut, measure='PLO', method='REML', data=gen[selrs])

  pdf(file=paste('forestRifRes_',i,'.pdf', sep=''), width=10, height=8)
  forest(fit.rr, atransf=transf.ilogit, refline=NA, addfit=TRUE,
         slab=gen$country[selrr],
         xlab=paste('Pooled Sensitivity, Rif Resistant (', gen$drug[selrr][1], ')', sep=''),
         order=order(-gen$se[selrr]))
  dev.off()

  pdf(file=paste('forestRifSus_',i,'.pdf', sep=''), width=10, height=8)
  forest(fit.rs, atransf=transf.ilogit, refline=NA, addfit=TRUE,
         slab=gen$country[selrs],
         xlab=paste('Pooled Sensitivity, Rif Susceptible (', gen$drug[selrs][1], ')', sep=''),
         order=order(-gen$se[selrs]))
  dev.off()
}



# store pooled results for non-injectables
idrugs <- c('kan','kan2','amk','cap')
nidrugs <- setdiff(drugs, idrugs)
res <- data.table(drug=drugs, se=NA, se.lo=NA, se.hi=NA, se.sd=NA, sp=NA, sp.lo=NA, sp.hi=NA, sp.sd=NA)

for (i in nidrugs){
  sel <- gen$drug==i & gen$patientGroup == 'All patients' # & gen$se>0 & gen$se<1
  sel2 <- res$drug==i

  sefit <- rma(xi=res.mut, mi=res.nomut, measure='PLO', data=gen[sel], method='REML')
  spfit <- rma(xi=sens.nomut, mi=sens.mut, measure='PLO', data=gen[sel], method='REML')

  res$se[sel2] <- invlogit(predict(sefit)$pred)
  res$se.lo[sel2] <- invlogit(predict(sefit)$ci.lb)
  res$se.hi[sel2] <- invlogit(predict(sefit)$ci.ub)

  res$sp[sel2] <- invlogit(predict(spfit)$pred)
  res$sp.lo[sel2] <- invlogit(predict(spfit)$ci.lb)
  res$sp.hi[sel2] <- invlogit(predict(spfit)$ci.ub)
}

# store pooled results for injectables
for (i in idrugs){
  sel <- gen$drug==i & gen$patientGroup == 'Rif resistant' # & gen$se>0 & gen$se<1
  sel2 <- res$drug==i

  sefit <- rma(xi=res.mut, mi=res.nomut, measure='PLO', data=gen[sel3], method='REML')
  spfit <- rma(xi=sens.nomut, mi=sens.mut, measure='PLO', data=gen[sel], method='REML')

  res$se[sel2] <- invlogit(predict(sefit)$pred)
  res$se.lo[sel2] <- invlogit(predict(sefit)$ci.lb)
  res$se.hi[sel2] <- invlogit(predict(sefit)$ci.ub)

  res$sp[sel2] <- invlogit(predict(spfit)$pred)
  res$sp.lo[sel2] <- invlogit(predict(spfit)$ci.lb)
  res$sp.hi[sel2] <- invlogit(predict(spfit)$ci.ub)
}

res$se.sd <- (res$se.hi - res$se.lo)/4
res$sp.sd <- (res$sp.hi - res$sp.lo)/4


# actual prevalence
gen$prevRes <- gen$res/gen$sampleSize
out <- cii(gen$sampleSize, gen$res)
gen$prevRes.lo <- out$lower95ci
gen$prevRes.hi <- out$upper95ci
gen$prevRes.sd <- out$se

vget.beta <- Vectorize(get.beta, c('ev','sd'))

out1 <- vget.beta(res$se, res$se.sd)
out2 <- vget.beta(res$sp, res$sp.sd)
res$se1 <- out1[1, ]
res$se2 <- out1[2, ]
res$sp1 <- out2[1, ]
res$sp2 <- out2[2, ]

sel <- (gen$patientGroup=="All patients" & gen$drug %in% nidrugs) | (gen$patientGroup=="Rif resistant" & gen$drug %in% idrugs)
gen2 <- merge(gen[sel,.(country, drug, tpos=res.mut+sens.mut, n=sampleSize, prevRes, prevRes.lo, prevRes.hi, prevRes.sd)],
              res, by='drug')


# posterior prevalence
library(rjags)

bcprev <- "model {
tpos ~ dbin(theta, n)
theta <- se*phi + (1-sp)*(1-phi)
se ~ dbeta(se1, se2)
sp ~ dbeta(sp1, sp2)
phi ~ dbeta(1, 1)
}"

set.seed(123)
for (i in 1:dim(gen2)[1]){
  print(i)

  out <- jags.model(textConnection(bcprev),
                    data = list(tpos = gen2$tpos[i],
                                n = gen2$n[i],
                                se1 = gen2$se1[i],
                                se2 = gen2$se2[i],
                                sp1 = gen2$sp1[i],
                                sp2 = gen2$sp2[i]
                    ),
                    n.chains = 3)

  update(out, 10000)

  mcmc <- coda.samples(out, variable.names=c("phi"), n.iter=20000)
  omcmc <- summary(mcmc)
  gen2$post[i] <- omcmc$statistics[1]
  gen2$post.lo[i] <- omcmc$quantiles[1]
  gen2$post.hi[i] <- omcmc$quantiles[5]
  gen2$post.sd[i] <- omcmc$statistics[2]
}


# plot posterior vs actual prevalence
druglab <- c('OFX','OFX (2)', 'MFX', 'MFX (2)', 'PZA', 'PZA (2)', 'RIF', 'INH', 'INH (2)', 'KAN', 'KAN (2)', 'AMK', 'CAP')

for (i in 1:length(drugs)){
  sel <- gen2$drug == drugs[i]
  pdf(file=paste('post_', drugs[i],'.pdf', sep=''), width=10, height=8)
  p <-  ggplot(data=gen2[sel]) +
    geom_point(aes(prevRes, country), shape=I('|'), size=I(5), colour=I('lightblue')) +
    geom_point(aes(post, country), shape=I('|'), size=I(5), colour=I('red')) +
    geom_segment(aes(x=prevRes.lo, xend=prevRes.hi, y=country, yend=country), colour=I('lightblue'), size=I(3)) +
    geom_segment(aes(x=post.lo, xend=post.hi, y=country, yend=country), colour=I('red'), size=I(1)) +
    geom_point(aes(tpos/n, country), shape=I(4)) +
    xlab(paste('Prevalence of ', druglab[i], ' resistance (blue=true, red=test positive, adjusted)', sep='')) +
    ylab('') + theme_bw(base_size = 18)
  print(p)
  dev.off()
}


save(gen2, file='gen2.Rdata')




# assume all genotypically resistant are also phenotypically resistant (sp = 1)
library(rjags)
load('gen2.Rdata')
load('gen.Rdata')

drugs <- unique(gen$drug)
idrugs <- c('kan','kan2','amk','cap')
nidrugs <- setdiff(drugs, idrugs)

gen$res2 <- gen$mutations + gen$res.nomut
gen$prevRes2 <- gen$res2/gen$sampleSize
out <- cii(gen$sampleSize, gen$res2)

gen$prevRes2.lo <- out$lower95ci
gen$prevRes2.hi <- out$upper95ci
gen$prevRes2.sd <- out$se


sel <- (gen$patientGroup=="All patients" & gen$drug %in% nidrugs) | (gen$patientGroup=="Rif resistant" & gen$drug %in% idrugs)
gen3 <- merge(gen2, gen[sel,.(country, drug, prevRes2, prevRes2.lo, prevRes2.hi, prevRes2.sd)], by=c('country', 'drug'))

bcprev2 <- "model {
tpos ~ dbin(theta, n)
theta <- se*phi
se ~ dbeta(se1, se2)
phi ~ dbeta(1, 1)
}"

set.seed(1234)
for (i in 1:dim(gen3)[1]){
  print(i)

  out <- jags.model(textConnection(bcprev2),
                    data = list(tpos = gen3$tpos[i],
                                n = gen3$n[i],
                                se1 = gen3$se1[i],
                                se2 = gen3$se2[i]
                    ),
                    n.chains = 3)

  update(out, 10000)

  mcmc <- coda.samples(out, variable.names=c("phi"), n.iter=20000)
  omcmc <- summary(mcmc)
  gen3$post2[i] <- omcmc$statistics[1]
  gen3$post2.lo[i] <- omcmc$quantiles[1]
  gen3$post2.hi[i] <- omcmc$quantiles[5]
  gen3$post2.sd[i] <- omcmc$statistics[2]
}

gen3[,.(country,drug,prevRes,prevRes2,se,sp,post,post2)]

# plot posterior vs actual prevalence
druglab <- c('OFX','OFX (2)', 'MFX', 'MFX (2)', 'PZA', 'PZA (2)', 'RIF', 'INH', 'INH (2)', 'KAN', 'KAN (2)', 'AMK', 'CAP')

for (i in 1:length(drugs)){
  sel <- gen3$drug == drugs[i]
  pdf(file=paste('post_sp1_', drugs[i],'.pdf', sep=''), width=10, height=8)
  p <-  ggplot(data=gen3[sel]) +
    geom_point(aes(prevRes2, country), shape=I('|'), size=I(5), colour=I('lightblue')) +
    geom_point(aes(post2, country), shape=I('|'), size=I(5), colour=I('red')) +
    geom_segment(aes(x=prevRes2.lo, xend=prevRes2.hi, y=country, yend=country), colour=I('lightblue'), size=I(3)) +
    geom_segment(aes(x=post2.lo, xend=post2.hi, y=country, yend=country), colour=I('red'), size=I(1)) +
    geom_point(aes(tpos/n, country), shape=I(4)) +
    xlab(paste('Prevalence of ', druglab[i], ' resistance (blue=true (including pheno- geno+), red=test positive, adjusted)', sep='')) +
    ylab('') + theme_bw(base_size = 18)
  print(p)
  dev.off()
}

save(gen3, file='gen3.Rdata')



# by rif status in selected drugs
rdrugs <- c('inh2','ofx2', 'mfx2', 'pza')

resRS <- data.table(drug=rdrugs, RR=0, se=NA, se.lo=NA, se.hi=NA, se.sd=NA)
resRR <- data.table(drug=rdrugs, RR=1, se=NA, se.lo=NA, se.hi=NA, se.sd=NA)

for (i in rdrugs){
  sel <- gen$drug==i & gen$patientGroup %in% c('Rif susceptible')
  sel2 <- resRS$drug==i

  sefit <- rma(xi=res.mut, mi=res.nomut, measure='PLO', data=gen[sel], method='REML')

  resRS$se[sel2] <- invlogit(predict(sefit)$pred)
  resRS$se.lo[sel2] <- invlogit(predict(sefit)$ci.lb)
  resRS$se.hi[sel2] <- invlogit(predict(sefit)$ci.ub)
}
for (i in rdrugs){
  sel <- gen$drug==i & gen$patientGroup %in% c('Rif resistant')
  sel2 <- resRR$drug==i

  sefit <- rma(xi=res.mut, mi=res.nomut, measure='PLO', data=gen[sel], method='REML')

  resRR$se[sel2] <- invlogit(predict(sefit)$pred)
  resRR$se.lo[sel2] <- invlogit(predict(sefit)$ci.lb)
  resRR$se.hi[sel2] <- invlogit(predict(sefit)$ci.ub)
}

res <- rbind(resRS, resRR)

res$se.sd <- (res$se.hi - res$se.lo)/4

vget.beta <- Vectorize(get.beta, c('ev','sd'))
out1 <- vget.beta(res$se, res$se.sd)
res$se1 <- out1[1, ]
res$se2 <- out1[2, ]


gen$res2 <- gen$mutations + gen$res.nomut
gen$prevRes <- gen$res/gen$sampleSize
gen$prevRes2 <- gen$res2/gen$sampleSize
out <- cii(gen$sampleSize, gen$res2)

gen$prevRes2.lo <- out$lower95ci
gen$prevRes2.hi <- out$upper95ci
gen$prevRes2.sd <- out$se

sel <- ((gen$patientGroup=="Rif susceptible" & gen$drug %in% nidrugs) | (gen$patientGroup=="Rif resistant")) & gen$drug %in% rdrugs
gen$RR[gen$patientGroup=='Rif susceptible'] <- 0
gen$RR[gen$patientGroup=='Rif resistant'] <- 1
gen4 <- merge(gen[sel,.(country, RR, patientGroup, drug, tpos=res.mut+sens.mut, n=sampleSize,
                        prevRes, prevRes2, prevRes2.lo, prevRes2.hi, prevRes2.sd)],
              res, by=c('drug', 'RR'))


bcprev2 <- "model {
tpos ~ dbin(theta, n)
theta <- se*phi
se ~ dbeta(se1, se2)
phi ~ dbeta(1, 1)
}"

set.seed(12345)
for (i in 1:dim(gen4)[1]){
  print(i)

  out <- jags.model(textConnection(bcprev2),
                    data = list(tpos = gen4$tpos[i],
                                n = gen4$n[i],
                                se1 = gen4$se1[i],
                                se2 = gen4$se2[i]
                    ),
                    n.chains = 3)

  update(out, 10000)

  mcmc <- coda.samples(out, variable.names=c("phi"), n.iter=20000)
  omcmc <- summary(mcmc)
  gen4$post2[i] <- omcmc$statistics[1]
  gen4$post2.lo[i] <- omcmc$quantiles[1]
  gen4$post2.hi[i] <- omcmc$quantiles[5]
  gen4$post2.sd[i] <- omcmc$statistics[2]
}

gen4[,.(country,drug,prevRes,prevRes2,se,post2)]

druglab2 <- c('INH (2)', 'OFX (2)', 'MFX (2)', 'PZA')

for (i in 1:length(rdrugs)){
  sel <- gen4$drug == rdrugs[i]
  pdf(file=paste('post_sp1_byRif_', rdrugs[i],'.pdf', sep=''), width=10, height=8)
  p <-  ggplot(data=gen4[sel]) +
    geom_point(aes(prevRes2, country), shape=I('|'), size=I(5), colour=I('lightblue')) +
    geom_point(aes(post2, country), shape=I('|'), size=I(5), colour=I('red')) +
    geom_segment(aes(x=prevRes2.lo, xend=prevRes2.hi, y=country, yend=country), colour=I('lightblue'), size=I(3)) +
    geom_segment(aes(x=post2.lo, xend=post2.hi, y=country, yend=country), colour=I('red'), size=I(1)) +
    geom_point(aes(tpos/n, country), shape=I(4)) +
    xlab(paste('Prevalence of ', druglab2[i], ' resistance (blue=true, red=test positive, adjusted)', sep='')) +
    ylab('') + theme_bw(base_size = 18) +
    facet_wrap(~patientGroup, scales='free_x')
  print(p)
  dev.off()
}

save(gen4, file='gen4.Rdata')





