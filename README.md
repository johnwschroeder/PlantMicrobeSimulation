README
================
John Schroeder
1/16/2020

# Introduction

The following R project includes code necessary to run the simulation
for the manuscript entitled “Mutualist and pathogen traits interact to
affect plant community structure in a spatially explicit model.” All
functions necessary to run the simulation are included in
Simulation\_functions.R.

The following R packages are required to run the simulation

``` r
library(doSNOW) #For parallelizing simulation runs
library(randomForest) #For conducting random forest analyses between batches of runs (used for optimization)
library(abind) #For binding lists of results
library(poweRlaw) #For generating dispersal kernels
library(vegan)
```

Source code for the simulation functions:

``` r
source("./Simulation_functions.R")
```

The following code chunks conduct an example simulation run using the
parameter values presented in Table 1 of Schroeder et al. (2020, in
review)

``` r
numcores <- 4 #Define the number of cores to use to simultaneously run multiple simulations
cl <- parallel::makeCluster(numcores) #Create cluster
registerDoSNOW(cl) #Register cluster
```

Define parameter values and settings for simulation runs (this example
runs simulations with identical parameter values four
times):

``` r
position <- as.matrix(cbind('g' = 8.939777,'h' = 0.4086455, 'b.t' = 2.124761, 'f' = 0.8,
               's.m' = 0.7367272, 'b.m' = 2.318665, 'alpha.m' = 1.164921,
               'gamma.m' = 0.2254993, 'r.m' = 0.6809609, 'q.m' = -0.07401046, 'c.m' = 1.677422,
               's.p' = 0.04722473, 'b.p' = 2.840781, 'alpha.p' = 0.8030246,
               'gamma.p' = 0.9312558, 'r.p' = 0.5455649, 'q.p' = 1.857254, 'c.p' = 1.442391))
particle.positions <- as.data.frame(position[rep(1,4),])
number.of.runs <- 4
tree.species <- 5
trees.in.forest <- 499
mort = 0.1
mortality.replacement.steps <- 3000
eliminate.feedback <- c(rep(FALSE,4))
fix.feedback <- c(rep(FALSE,4))
remove.mutualists <- c(rep(FALSE,4))
remove.pathogens <- c(rep(FALSE,4))
f.vals <- t(sapply(particle.positions[,4],function(x) {(c(1:tree.species-1)/(tree.species-1))*(1-x)+x}))
pb <- txtProgressBar(max = number.of.runs, style = 3)
```

    ## 
      |                                                                       
      |                                                                 |   0%

``` r
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
```

Run multiple simulations using the ‘dopar’ function, and store the
output in a list called simulation.results:

``` r
simulation.results <- foreach(g=particle.positions$g,
                                  h=particle.positions$h,
                                  b.t=particle.positions$b.t,
                                  s.m = particle.positions$s.m,
                                  b.m = particle.positions$b.m,
                                  alpha.m=particle.positions$alpha.m,
                                  gamma.m=particle.positions$gamma.m,
                                  r.m=particle.positions$r.m,
                                  q.m=particle.positions$q.m,
                                  c.m=particle.positions$c.m,
                                  s.p = particle.positions$s.p,
                                  b.p = particle.positions$b.p,
                                  alpha.p=particle.positions$alpha.p,
                                  gamma.p =particle.positions$gamma.p,
                                  r.p = particle.positions$r.p,
                                  q.p = particle.positions$q.p,
                                  c.p = particle.positions$c.p,
                                  fix.feedback = fix.feedback,
                                  eliminate.feedback = eliminate.feedback,
                                  remove.mutualists = remove.mutualists,
                                  remove.pathogens = remove.pathogens,
                                  index=c(1:number.of.runs),
                                  .packages=c("poweRlaw","vegan"),.options.snow = opts) %dopar%
        (
        psf.simulation(m=trees.in.forest,
                            mort=mort,
                            tree.species = tree.species,
                            mutualist.species.per.tree = 1,
                            pathogen.species.per.tree = 1,
                            time.steps=mortality.replacement.steps,
                            mutualist.effect.function.consp=mutualist.effect.function.consp, 
                            mutualist.effect.function.heterosp=mutualist.effect.function.heterosp,
                            pathogen.effect.function.consp=pathogen.effect.function.consp,
                            pathogen.effect.function.heterosp=pathogen.effect.function.heterosp,
                            g=g,
                            h=h,
                            s.p=s.p,
                            s.m=s.m,
                            b.p=b.p,
                            b.m=b.m,
                            b.t=b.t,
                            f.vals= rev(f.vals[index,]),
                            index=index,
                            gamma.m=gamma.m,
                            gamma.p=gamma.p,
                            r.m=r.m,
                            r.p=r.p,
                            q.m=q.m,
                            q.p=q.p,
                            c.m=c.m,
                            c.p=c.p,
                            alpha.m=alpha.m,
                            alpha.p=alpha.p,
                            fix.feedback = fix.feedback,
                            eliminate.feedback = eliminate.feedback,
                            remove.mutualists = remove.mutualists,
                            remove.pathogens = remove.pathogens,
                            initiate.forest.matrix=initiate.forest.matrix,
                            initiate.fungal.matrix=initiate.fungal.matrix,
                            trial.function=trial.function,
                            dispersal.function=dpldis,
                            microbe.dispersal.function=dpldis,
                            track.over.time=TRUE))
```

    ## 
      |                                                                       
      |================                                                 |  25%
      |                                                                       
      |================================                                 |  50%
      |                                                                       
      |=================================================                |  75%
      |                                                                       
      |=================================================================| 100%
