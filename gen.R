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


# calculate SDs and 95% CI assuming binomial errors
out <- cii(gen$res.mut + gen$res.nomut, gen$res.mut)
gen$se.sd <- out$se
gen$se.lo <- out$lower95ci
gen$se.hi <- out$upper95ci

out2 <- cii(gen$sens.mut + gen$sens.nomut, gen$sens.nomut)
gen$sp.sd <- out2$se
gen$sp.lo <- out2$lower95ci
gen$sp.hi <- out2$upper95ci

# save data
save(gen, file = 'gen.Rdata')
#----------------------------------------


# forest plots
load('gen.Rdata')
sel <- gen$drug=='rif' & gen$patientGroup != 'All patients'
forest(x=gen$se[sel], sei=gen$se.sd[sel], ci.lb=gen$se.lo[sel], ci.ub=gen$se.hi[sel],
       slab=paste(gen$country[sel], '/', gen$patientGroup[sel]),
       xlab=paste('Sensitivity (', gen$drug[sel][1], ')', sep=''))


