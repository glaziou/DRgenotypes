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
gen <- fread('gen070317b.csv')
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


# save data
save(gen, file = 'gen.Rdata')
#----------------------------------------


# basic forest plots
load('gen.Rdata')

drugs <- unique(gen$drug)
gen$formerSoviet <- gen$country %in% c('AZE', 'BLR', 'UKR')

for (i in drugs){
  sel <- gen$drug==i & gen$patientGroup != 'All patients'

  pdf(file=paste('forestGroups_',i,'.pdf', sep=''), width=10, height=8)
  forest(x=gen$se[sel], sei=gen$se.sd[sel], ci.lb=gen$se.lo[sel], ci.ub=gen$se.hi[sel],
         slab=paste(gen$country[sel], '/', gen$patientGroup[sel]),
         xlab=paste('Sensitivity (', gen$drug[sel][1], ')', sep=''), subset=order(gen$country[sel]))
  dev.off()
}



# effect modifiers
mods <- data.table(drugs=drugs, testHistory=NA, testSoviet=NA, testRifRes=NA)
for (i in drugs){
  sel.all <- gen$drug==i & gen$patientGroup %in% c('All patients')
  sel.hx <- gen$drug==i & gen$patientGroup %in% c('New patients','Previously treated')

  fit.hx <- rma(xi=res.mut, ni=res, mods=~patientGroup, measure='PLO', data=gen[sel.hx], method='REML')
  fit.sov <- rma(xi=res.mut, ni=res, mods=~formerSoviet, measure='PLO', data=gen[sel.all], method='REML')

  if (i %ni% c('rif_h','rif_m')) {
    sel.rif <- gen$drug==i & gen$patientGroup %in% c('Rif resistant','Rif susceptible')
    fit.rif <- rma(xi=res.mut, ni=res, mods=~factor(patientGroup=='Rif resistant'), measure='PLO', data=gen[sel.rif], method='REML')
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


# for (i %in% c('inh','inh2','rif')){
#   fit.ns <- rma(yi=se, sei=se.sd, measure='PLO', method='REML', data=gen[drug==i & patientGroup=='New patients'])
#   fit.s <- rma(yi=se, sei=se.sd, measure='PLO', method='REML', data=gen[drug==i & patientGroup=='Previously treated'])
#   pdf(file=paste('forestByHistory_',i,'.pdf', sep=''), width=10, height=8)
#   forest(fit.hx, atransf=transf.ilogit, refline=NA, addfit=FALSE,
#          slab=gen$country[sel],
#          rows=c(18:13, 8:3), ylim=c(1, 21),
#          xlab=paste('Pooled Sensitivity (', gen$drug[sel][1], ')', sep=''),
#          order=order(-gen$formerSoviet[sel]))
#   addpoly(fit.ns, row= 1.5, atransf=transf.ilogit, mlab="Pooled, New patients")
#   addpoly(fit.s, row= 9.5, atransf=transf.ilogit, mlab="Pooled, Previously treated")
#   dev.off()
# }



# sel <- gen$drug==i
# forest(fit.hx, atransf=transf.ilogit, refline=NA, addfit=FALSE,
#        slab=gen$country[sel],
#        xlab=paste('Pooled Sensitivity (', gen$drug[sel][1], ')', sep=''),
#        order=order(-gen$formerSoviet[sel]))
# addpoly(fit.ns, row= 1.5, atransf=transf.ilogit, mlab="Pooled, New patients")
# addpoly(fit.s, row= 9.5, atransf=transf.ilogit, mlab="Pooled, Previously treated")



# pooled separately in former soviet countries
for (i in c('inh_h', 'inh2_h', 'rif_h', 'inh_m', 'inh2_m', 'rif_m')){
  sel <- gen$drug==i & gen$patientGroup == 'All patients'
  fit <- rma(yi=se, sei=se.sd, measure='PLO', method='REML', data=gen[sel])
  fit.ns <- rma(yi=se, sei=se.sd, measure='PLO', method='REML', data=gen[sel & formerSoviet==FALSE])
  fit.s <- rma(yi=se, sei=se.sd, measure='PLO', method='REML', data=gen[sel & formerSoviet==TRUE])
  pdf(file=paste('forestAll_',i,'.pdf', sep=''), width=10, height=8)
  forest(fit, atransf=transf.ilogit, refline=NA, addfit=FALSE,
             slab=gen$country[sel],
             rows=c(13:11, 6:3), ylim=c(1, 16),
             xlab=paste('Pooled Sensitivity (', gen$drug[sel][1], ')', sep=''),
             order=order(-gen$formerSoviet[sel]))
  addpoly(fit.s, row= 9.5, atransf=transf.ilogit, mlab="Pooled, former Soviet Union")
  addpoly(fit.ns, row= 1.5, atransf=transf.ilogit, mlab="Pooled, all other countries")
  dev.off()
}

for (i in c('ofx', 'ofx2')){
  sel <- gen$drug==i & gen$patientGroup == 'All patients' & gen$country %ni% c('PHL')
  fit <- rma(yi=se, sei=se.sd, measure='PLO', method='REML', data=gen[sel])
  fit.ns <- rma(yi=se, sei=se.sd, measure='PLO', method='REML', data=gen[sel & formerSoviet==FALSE])
  fit.s <- rma(yi=se, sei=se.sd, measure='PLO', method='REML', data=gen[sel & formerSoviet==TRUE])
#  pdf(file=paste('forestAll_',i,'.pdf', sep=''), width=10, height=8)
  forest(fit, atransf=transf.ilogit, refline=NA, addfit=FALSE,
         slab=gen$country[sel],
         rows=c(12:10, 5:3), ylim=c(1, 16),
         xlab=paste('Pooled Sensitivity (', gen$drug[sel][1], ')', sep=''),
         order=order(-gen$formerSoviet[sel]))
  addpoly(fit.s, row= 8.5, atransf=transf.ilogit, mlab="Pooled, former Soviet Union")
  addpoly(fit.ns, row= 1.5, atransf=transf.ilogit, mlab="Pooled, all other countries")
#  dev.off()
}


# pooled separately by rif res
for (i in c('inh_h', 'inh2_h', 'inh_m', 'inh2_m', 'ofx','ofx2')){
  selrr <- gen$drug==i & gen$patientGroup == 'Rif resistant'
  selrs <- gen$drug==i & gen$patientGroup == 'Rif susceptible'
  fit.rr <- rma(yi=se, sei=se.sd, measure='PLO', method='REML', data=gen[selrr])
  fit.rs <- rma(yi=se, sei=se.sd, measure='PLO', method='REML', data=gen[selrs])
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

for (i in setdiff(drugs, c('rif_h','rif_m','kan','kan2'))){
  selrr <- gen$drug==i & gen$patientGroup == 'Rif resistant'
  selrs <- gen$drug==i & gen$patientGroup == 'Rif susceptible'
  fit.rr <- rma(yi=se, sei=se.sd, measure='PLO', method='REML', data=gen[selrr])
  fit.rs <- rma(yi=se, sei=se.sd, measure='PLO', method='REML', data=gen[selrs])
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

for (i in c('kan','kan2')){
  selrr <- gen$drug==i & gen$patientGroup == 'Rif resistant' & gen$se>0
  selrs <- gen$drug==i & gen$patientGroup == 'Rif susceptible' & gen$se>0
  fit.rr <- rma(yi=se, sei=se.sd, measure='PLO', method='REML', data=gen[selrr])
  fit.rs <- rma(yi=se, sei=se.sd, measure='PLO', method='REML', data=gen[selrs])
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








