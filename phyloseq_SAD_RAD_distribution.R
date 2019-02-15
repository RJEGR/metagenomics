# Motivation
# The distribution of abundances of species is one of the basic patterns of ecological communities. The empirical distributions of abundances (SADs) or their ranks (RADs) are modelled through probability distributions. Hence, the maximum likelihood method can be used to fit and compare such mode.m.


# Species abundance distributions (SADs) are one of the basic patterns of ecological communities
# (McGill et al., 2007). The empirical distributions are traditionally modeled through
# probability distributions. Hence, the maximum likelihood method can be used to fit and
# compare competing models for SADs. The package sads provides functions to fit the most
# used models to empirical SADs. Species-abundance models assign a probability for each abundance value. 
# Thus, these models are probability density functions (PDFs) of abundances of species. Rank-abundance
# models assign a probability for each abundance rank

# install.packages('sads')
#library(devtools)
#install_github(repo = 'piLaboratory/sads', ref= 'dev', build_vignettes = TRUE)

library(ggplot2)

phyloseq = prune_samples(sample_sums(phyloseq) > 0, phyloseq)
# What about the total reads per sample, and what does the distribution look like?

df = data.frame(nreads = sort(taxa_sums(phyloseq), TRUE), sorted = 1:ntaxa(phyloseq), 
                type = "OTUs")
df = rbind(df, data.frame(nreads = sort(sample_sums(phyloseq), 
                                        TRUE), sorted = 1:nsamples(phyloseq), type = "Samples"))
title = "Total number of reads"
p = ggplot(df, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity") + scale_x_continuous(expand = c(-0,0))
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")



library(sads)

A <- sort(taxa_sums(phyloseq::subset_samples(phyloseq, Transecto == "A")))
Y <- sort(taxa_sums(phyloseq::subset_samples(phyloseq, Transecto == "Y")))
H <- sort(taxa_sums(phyloseq::subset_samples(phyloseq, Transecto == "H")))


# Function octav tabulates the number of species in classes of logarithm of abundances at
# base 2 (Preston’s octaves) and returns a data frame
# A logical argument preston allows smoothing the numbers as proposed by Preston (1948).
# The octave number is the upper limit of the class in log2 scale. Hence, for abundance values
# smaller than one (e.g. biomass data) the octave numbers are negative. A Preston plot is a
# histogram of this table, obtainable by applying the function plot to the data frame
# 

sad <- df[,1]
(sad.oc <- octav(sad))
plot(sad.oc, ylab="Number of Taxa")

# The plot method for objects of class octav has a logical argument prop that rescales the yaxis
# to relative frequencies of species (taxas) in each octave, which can be used to compare different
# data sets:

plot(sad.oc, prop = TRUE, border=NA, col=NA)
lines(octav(Y), mid = FALSE, prop = TRUE, col="red")
lines(octav(A), mid = FALSE, prop = TRUE, col="blue")
lines(octav(H), mid = FALSE, prop = TRUE, col="green")
legend("topright", c("Transect Y", "Transect A", "Transect H" ), col=c("red", "blue", "green"), lty=1, bty="n")


# Function rad returns a data frame of sorted abundances and their ranks 2

head(sad.rad <- rad(sad))

# To get the rank-abundance or Whitaker’s plot apply the function plot on the data frame:
plot(sad.rad, ylab="Number of Reads", xlab = "Otus")

# and plot all samples
head(A.rad <- rad(A))
plot(A.rad, prop = TRUE, type="n", ylab="Taxa Relative Abundance", xlab = "Otus Rank")
lines(rad(Y), prop = TRUE, col="red")
lines(rad(A), prop = TRUE, col="blue")
lines(rad(H), prop = TRUE, col="green")
legend("topright", c("Transect Y", "Transect A", "Transect H" ), col=c("red", "blue", "green"), lty=1, bty="n")

###############
# model fitting
###############

(sad.m <- fitsad(sad,'lnorm'))
summary(sad.m)


summary(sad.m)
coef(sad.m)
logLik(sad.m)
AIC(sad.m)

sad.m.prf <- profile(sad.m)
likelregions(sad.m.prf) #likelihood intervals


confint(sad.m.prf)

# Then use plotprofmle to plot likelihood profiles at the original scale (relative negative
# log-likelihood) and function plot to get plots at chi-square scale (square-root of twice the relative log-likelihood)

par(mfrow=c(2,2))
plotprofmle(sad.m.prf)# log-likelihood profile
plot(sad.m.prf)# z-transformed profile
par(mfrow=c(1,1))


# When applied on a sad model object, the function plot returns four diagnostic plots:
par(mfrow=c(2,2))
plot(sad.m)
par(mfrow=c(1,1))

# It’s possible to fit other models to the same data set, such as the Poisson-lognormal and a truncated lognormal:

(sad.pl <- fitsad(x=sad, sad="poilog")) #default is zero-truncated
(sad.ln <- fitsad(x=sad, sad="lnorm", trunc=0.5)) # lognormal truncated at 0.5
# (sad.gam <- fitsad(x=sad, sad="gamma"))
(sad.geo <- fitsad(x=sad, sad="geom"))

# Null, Preemption, Lognormal, Zipf & ZipfMandelbrot
AICtab(sad.pl, sad.ln, sad.geo, base=TRUE)



head(sad.pl.oc <- octavpred(sad.pl))
head(sad.ln.oc <- octavpred(sad.ln))
head(sad.geo.oc <- octavpred(sad.geo))


# (moths.oc <- octav(moths))

plot(sad.oc)
lines(sad.pl.oc, col="blue")
lines(sad.ln.oc, col="red")
lines(sad.geo.oc, col="green")

legend("topright",
       c("Poisson-lognormal", "Truncated lognormal", "Geometric"),
       lty=1, col=c("blue","red", "green"))

#
head(sad.pl.rad <- radpred(sad.pl))

# Note: extreme values generated by radpred. Calculations will take a while...
# Error in qfinder(dpoilog, p[1], list(mu = mu, sig = sig), 0) : quantile function did not converge!
head(sad.ln.rad <- radpred(sad.ln))
head(sad.geo.rad <- radpred(sad.geo))


plot(sad.rad)
#lines(sad.pl.rad, col="blue")
lines(sad.ln.rad, col="red")
lines(sad.geo.rad, col="green")
legend("topright",
         c("Truncated lognormal", "Geometric"),
         lty=1, col=c("red", "green"))
