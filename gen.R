#########################################
# genotyping - metaanalysis
# PG, init 14/02/2017
#########################################
library(metafor)

setwd('../data/genotypes_2017/')

#----------------------------------------
# run this block only if dataset updated
# read data and drop rows with no useful data
gen <- fread('gen.csv')
sel <- gen$res.mut + gen$res.nomut > 0
gen <- gen[sel]
sel <- gen$country != 'Global'
gen <- gen[sel]

gen$res <- rowSums(cbind(gen$res.mut, gen$res.nomut))
gen$sens <- rowSums(cbind(gen$sens.mut, gen$sens.nomut))

# calculate SDs and 95% CI assuming binomial errors
out <- cii(gen$res, gen$res.mut)
gen$se.sd <- out$se
gen$se.lo <- out$lower95ci
gen$se.hi <- out$upper95ci

out2 <- cii(gen$sens, gen$sens.nomut)
gen$sp.sd <- out2$se
gen$sp.lo <- out2$lower95ci
gen$sp.hi <- out2$upper95ci


# save data
save(gen, file = 'gen.Rdata')
#----------------------------------------


# forest plots
load('gen.Rdata')

drugs <- unique(gen$drug)
mods <- data.table(drugs=drugs, testGroup=NA, testSoviet=NA)
gen$formerSoviet <- gen$country %in% c('AZE', 'BLR', 'UKR')

for (i in drugs){
  sel <- gen$drug==i & gen$patientGroup != 'All patients'

  pdf(file=paste('forestGroups_',i,'.pdf', sep=''), width=10, height=8)
  forest(x=gen$se[sel], sei=gen$se.sd[sel], ci.lb=gen$se.lo[sel], ci.ub=gen$se.hi[sel],
         slab=paste(gen$country[sel], '/', gen$patientGroup[sel]),
         xlab=paste('Sensitivity (', gen$drug[sel][1], ')', sep=''), subset=order(gen$country[sel]))
  dev.off()
}

for (i in drugs){
  sel <- gen$drug==i & gen$patientGroup != 'All patients'
  fit.g <- rma(xi=res.mut, ni=res, mods=~patientGroup, measure='PLO', data=gen[sel], method='REML')
  fit.s <- rma(xi=res.mut, ni=res, mods=~formerSoviet, measure='PLO', data=gen[sel], method='REML')
  mods$testGroup[mods$drugs==i] <- fit.g$QMp
  mods$testSoviet[mods$drugs==i] <- fit.s$QMp
}
(mods) # patient group effect significant for kan only

for (i in c('inh', 'inh2', 'rif')){
  sel <- gen$drug==i & gen$patientGroup == 'All patients'
  fit <- rma(yi=se, sei=se.sd, measure='PLO', method='REML', data=gen[sel])
  fit.ns <- rma(yi=se, sei=se.sd, measure='PLO', method='REML', data=gen[sel & formerSoviet==FALSE])
  fit.s <- rma(yi=se, sei=se.sd, measure='PLO', method='REML', data=gen[sel & formerSoviet==TRUE])
  pdf(file=paste('forestAll_',i,'.pdf', sep=''), width=10, height=8)
  forest(fit, atransf=transf.ilogit, refline=NA,
             slab=gen$country[sel],
             rows=c(14:12, 7:3), ylim=c(-1, 17),
             xlab=paste('Pooled Sensitivity (', gen$drug[sel][1], ')', sep=''),
             order=order(-gen$formerSoviet[sel]))
  addpoly(fit.s, row= 10.5, atransf=transf.ilogit, mlab="Model for countries of the former Soviet Union")
  addpoly(fit.ns, row= 1.5, atransf=transf.ilogit, mlab="Model for all other countries")
  dev.off()
}







