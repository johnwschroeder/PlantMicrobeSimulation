#John Schroeder
#16-Jan-2019
#Spatially explicit plant-microbe interactions simulation

#######################################################################################################################################################################################################################################################
#Parts of the simulation
#######################################################################################################################################################################################################################################################

#######################
#initiate.forest.matrix
#returns an m by n matrix with a tree community in which each tree starts with an equal abundance
#######################
initiate.forest.matrix <- function(m, #length of forest vector
                                   tree.species, #total number of tree species
                                   mutualist.species.per.tree, #total number of mutualist species per tree
                                   pathogen.species.per.tree #total number of pathogen species per tree
) { #Create data frame with equal numbers of individuals per species.
  X = c(1:m) #Position along x axis... randomized
  rankAbundance <- round(m/tree.species) # #Randomly assign abundances
  while (sum(rankAbundance)!=m) { #Just to make sure the total number of trees equals the number of cells
    if (sum(rankAbundance)>m) {
      rankAbundance[1] <- rankAbundance[1]-1
    }
    if (sum(rankAbundance)<m) {
      rankAbundance[1] <- rankAbundance[1]+1
    }
  }
  TreeSpecies <- rep(c(1:tree.species),rankAbundance) #Create species ID integer
  ID.num <- c(1:m) #Create individual ID number for each tree
  community.matrix <- data.frame("TreeSpecies"=sample(TreeSpecies,m,replace=FALSE),"X"=X) #Create data frame without fungal matrix
  mutualist.matrix <- mat.or.vec(m,mutualist.species.per.tree*tree.species) #Create empty matrix to be filled with mutualist abundance values
  colnames(mutualist.matrix) <- paste(
    rep(paste("Mu",c(1:mutualist.species.per.tree),sep="."),tree.species),
    paste("TSp",rep(c(1:tree.species),each=mutualist.species.per.tree),sep="."),sep="_") #Create labels for each mutualist species
  pathogen.matrix <- mat.or.vec(m,pathogen.species.per.tree*tree.species) #Create empty matrix to be filled with pathogen abundance values
  colnames(pathogen.matrix) <- paste(
    rep(paste("Pa",c(1:pathogen.species.per.tree),sep="."),tree.species),
    paste("TSp",rep(c(1:tree.species),each=pathogen.species.per.tree),sep="."),sep="_") #Create labels for each pathogen species
  for (i in c(1:length(community.matrix[,1]))) { #The following for loop assigns a random mutualist/pathogen to each tree
    mu <- sample(c(1:mutualist.species.per.tree),size=1)
    mutualist.matrix[i,which(colnames(mutualist.matrix)==paste(paste("Mu.",mu,sep=""),paste("TSp",community.matrix[i,1],sep="."),sep="_"))] <- 1
    pa <- sample(c(1:pathogen.species.per.tree),size=1)
    pathogen.matrix[i,which(colnames(pathogen.matrix)==paste(paste("Pa.",pa,sep=""),paste("TSp",community.matrix[i,1],sep="."),sep="_"))] <- 1
  }
  
  output <- list()
  output$mutualists.adult <- mutualist.matrix
  output$pathogens.adult <- pathogen.matrix
  output$forest.matrix <- community.matrix
  output$rankAbundance <- rankAbundance
  return(output)
}

#######################
#initiate.fungal.matrix
#returns a matrix of the affinity of each fungal taxon for each host
#######################
initiate.fungal.matrix <- function(tree.species, #total number of plant species
                                   fungal.species.per.tree, #number of fungal species per plant (either mutualist of pathogen)
                                   Mu.or.p, #character either "Mu" or "Pa" for "mutualist" or "pathogen"
                                   conspecific.effect.function, #Function to determine the relationship between the affinity value and the host effect value on conspecific hosts
                                   heterospecific.effect.function, #Function to determine the relationship between the affinity value and the host effect value on heterospecific hosts
                                   heterospecific.affinity #The affinity that a given microbe has on a non-preferred host
) {
  affinities <- c(1:(fungal.species.per.tree))/(fungal.species.per.tree) #Allows for multiple fungal taxa per tree. Equal to 1 for all simulations presented in the manuscript
  conspecific.effects <- conspecific.effect.function(affinities,heterospecific.affinity) #Uses conspecific effect function to determine the effect values for microbes on their preferred hosts
  heterospecific.effects <- heterospecific.effect.function(affinities,heterospecific.affinity) #Uses heterospecific effect function to determine the effect values for microbes on their non-preferred hosts
  
  effect.matrix.combined <- mat.or.vec(fungal.species.per.tree*tree.species, tree.species) #Empty matrix of effect values
  rownames(effect.matrix.combined) <- paste(
    rep(paste(Mu.or.p,c(1:fungal.species.per.tree),sep="."),tree.species),
    paste("TSp",rep(c(1:tree.species),each=fungal.species.per.tree),sep="."),sep="_") #Add names of fungal species to rows
  colnames(effect.matrix.combined) <- paste("TSp",c(1:tree.species),sep=".") #Add names of tree species to columns
  
  for (i in 1:tree.species) { 
    effect.matrix <- mat.or.vec(fungal.species.per.tree, tree.species) #Empty matrix
    effect.matrix[,i] <- conspecific.effects #Assign the maximum effect to the row pointing to the focal species
    effect.matrix[,-i] <- heterospecific.effects #Indicates heterospecific host
    effect.matrix.combined[((i-1)*fungal.species.per.tree+1):(i*fungal.species.per.tree),] <- effect.matrix
  } #Create distribution of effect values
  return(effect.matrix.combined)
}

#######################
#Effect functions
#These relate host affinity to the effect value on a given host
#######################
mutualist.effect.function.consp <- function(mutualist.affinities, tradeoff.intercept){mutualist.affinities}
mutualist.effect.function.heterosp <- function(mutualist.affinities, tradeoff.intercept){tradeoff.intercept*mutualist.affinities}
pathogen.effect.function.consp <- function(pathogen.affinities, tradeoff.intercept){pathogen.affinities}
pathogen.effect.function.heterosp <- function(pathogen.affinities, tradeoff.intercept){tradeoff.intercept*pathogen.affinities}

#######################
#trial.function
#Converts total mutualist and pathogen effects to a seedling survival probability
#######################
trial.function <- function(mutualist.effect,
                           pathogen.effect,g,h) 
{ 
  1/(1+exp(-g*(h*mutualist.effect-pathogen.effect)))
}


#######################################################################################################################################################################################################################################################
#Functions to run the simulation
#######################################################################################################################################################################################################################################################

#######################
#psf.simulation
#returns the results of a single simulation run
#######################
psf.simulation <- function(m, #length of forest vector
                           mort, #percentage of trees that die at each time step
                           tree.species, #Number of tree species in the simulation
                           mutualist.species.per.tree, #Each tree species has an equal number of mutualists
                           pathogen.species.per.tree, #Each tree species has an equal number of pathogens
                           time.steps, #Number of mortality/replacement events
                           g, #Constant scaling the impact of microbes on host survival
                           h, #Constant scaling the relative effect of mutualists on host survival
                           s.m, #Affinity of mutualists for their non-preferred hosts
                           s.p, #Affinity of pathogens for their non-preferred hosts
                           b.m, #Mutualist dispersal parameter
                           b.p, #Pathogen dispersal parameter
                           b.t, #Plant dispersal parameter
                           f.vals, #Vector of survival probability multipliers per species
                           gamma.m, #Fecundity of mutualists
                           gamma.p, #Fecundity of pathogens
                           r.m, #Intrinsic growth rate of mutualists
                           r.p, #Intrinsic growht rate of pathogens
                           c.m, #Exponent that scales mutualist competitive ability with host affinity
                           c.p, #Exponent that scales pathogen competitive ability with host affinity
                           q.m, #Exponent that scales mutualist intrinsic growth rate with host affinity
                           q.p, #Exponent that scales pathogen intrinsic growth rate with host affinity
                           alpha.m, #Mutualist competition value
                           alpha.p, #Pathogen competition value
                           track.over.time, #True or False - whether to track populations of mutualists, pathogens, and plants through time, and save results to output
                           index,
                           fix.feedback = FALSE, #True or False - whether to fix the species-specific microbial communities half-way through the simulation
                           eliminate.feedback = FALSE, #TRUE or FALSE - whether to eliminate feedback halfway through the simulation
                           remove.mutualists = FALSE, #TRUE or FALSE - whether to eliminate effect of mutualists halfway through the simulation
                           remove.pathogens = FALSE, #TRUE or FALSE - whether to eliminate effect of pathogens halfway through the simulation
                           subsample = FALSE, #TRUE or FALSE - whether to only take data from every 10 time steps
                           ################################ functions
                           initiate.forest.matrix,
                           initiate.fungal.matrix,
                           trial.function, #Function determining relationship between microbial abundances and seedling survival probability
                           dispersal.function, #Dispersal kernel function for plants
                           microbe.dispersal.function, #Dispersal kernel function for pathogens and mutualists
                           mutualist.effect.function.consp, #Relationship between host affinity and benefit to the preferred host
                           mutualist.effect.function.heterosp, #Relationship between host affinity and benefit to the non-preferred host
                           pathogen.effect.function.consp, #Relationship between host affinity and impact on the preferred host
                           pathogen.effect.function.heterosp #Relationship between host affinity and impact on the non-preferred host
) {
  exit.status <- NA #Create variable to indicate whether the simulation stopped because of an error
  init <- initiate.forest.matrix(m,tree.species,mutualist.species.per.tree,pathogen.species.per.tree) #Create forest matrix
  forest.matrix <- init$forest.matrix #Plants IDs and positions
  mutualists.adult <- init$mutualists.adult #Mutualist IDs, positions, and abundancnes
  pathogens.adult <- init$pathogens.adult #Pathogen IDs, positions, and abundancnes
  
  distance.vec <- c(0:(floor(m/2)),(floor(m/2)):1) #Distance values between one point and all other points (with wrap-around)
  distance.pdf.mutualist <- microbe.dispersal.function(c(1:(floor(m/2)+1)),1,b.m) #PDF of mutualist abundances (dispersal kernel)
  distance.pdf.pathogen <- microbe.dispersal.function(c(1:(floor(m/2)+1)),1,b.p) #PDF of pathogen abundances (dispersal kernel)
  distance.pdf.plants <- dispersal.function(c(1:(floor(m/2)+1)),1,b.t) #PDF of seed distributions (dispersal kernel)
  distance.matrix <- mat.or.vec(m,m) #Initialize distance matrix between all cells in the vector
  matrix.pdf.mutualist <- mat.or.vec(m,m) #Initialize the relative abundances of mutualist spores from each plant in each cell
  matrix.pdf.pathogen <- mat.or.vec(m,m) #Initialize the relative abundances of pathogen spores from each plant in each cell
  matrix.pdf.plants <- mat.or.vec(m,m) #Initialize the relative abundances of plant seeds from each plant in each cell
  
  distance.matrix[1,] <- distance.vec #First row of distance matrix
  for (i in c(1:(m-1))) { 
    distance.matrix[i+1,] <- distance.vec[c((m-i+1):m,1:(m-i))]
  } #Fill in pairwise distance matrix
  
  for (i in c(1:nrow(matrix.pdf.mutualist))) {
    matrix.pdf.mutualist[i,] <- distance.pdf.mutualist[distance.matrix[i,]+1]/sum(distance.pdf.mutualist[distance.matrix[i,]+1])
    matrix.pdf.pathogen[i,] <- distance.pdf.pathogen[distance.matrix[i,]+1]/sum(distance.pdf.pathogen[distance.matrix[i,]+1])
    matrix.pdf.plants[i,] <- distance.pdf.plants[distance.matrix[i,]+1]/sum(distance.pdf.plants[distance.matrix[i,]+1])
  } #Use dispersal kernels to calculate density of spores/seeds of each microbe/plant at each location
  
  tree.matrix <- mat.or.vec(m,tree.species) #Initialize tree matrix
  fact.matrix <- diag(1,tree.species,tree.species) #Initialize diagonal matrix of 1s
  
  for (i in c(1:length(forest.matrix[,1]))) { #Create matrix version of tree species vector
    tree.matrix[i,] <- fact.matrix[forest.matrix[i,1],]
  }
  
  mutualist.propagule.matrix <- gamma.m*matrix.pdf.mutualist%*%mutualists.adult #Multiply fecundity by the dispersal kernel and the abundance of mutualists on adults
  pathogen.propagule.matrix <- gamma.p*matrix.pdf.pathogen%*%pathogens.adult #Multiply fecundity by the dispersal kernel and the abundance of pathogens on adults
  
  seed.matrix <- matrix.pdf.plants%*%tree.matrix #Abundance of each tree's seeds across the plot
  mutualist.effect.matrix <- initiate.fungal.matrix(tree.species,
                                                    fungal.species.per.tree=mutualist.species.per.tree, 
                                                    Mu.or.p="Mu", 
                                                    conspecific.effect.function=mutualist.effect.function.consp,
                                                    heterospecific.effect.function=mutualist.effect.function.heterosp,
                                                    heterospecific.affinity=s.m)
  mutualist.growth.matrix <- r.m*mutualist.effect.matrix^q.m #Create growth constant for each mutualist on each host
  mutualist.competition.matrix <- alpha.m*mutualist.effect.matrix^c.m #Create competition coefficients for each mutualist on each host
  pathogen.effect.matrix <- initiate.fungal.matrix(tree.species, 
                                                   fungal.species.per.tree=pathogen.species.per.tree, 
                                                   Mu.or.p="Pa", 
                                                   conspecific.effect.function=pathogen.effect.function.consp,
                                                   heterospecific.effect.function=pathogen.effect.function.heterosp,
                                                   heterospecific.affinity=s.p)
  pathogen.growth.matrix <- r.p*pathogen.effect.matrix^q.p #Create growth constant for each pathogen on each host
  pathogen.competition.matrix <- alpha.p*pathogen.effect.matrix^c.p #Create competition coefficients for each pathogen on each host
  if (track.over.time==TRUE) { #Initialize variables to track over time
    trees.over.time <- array(numeric(0),c(m,tree.species,time.steps)) #Array to store tree community over time
    trees.over.time[,,1] <- tree.matrix #Initial tree community
    colnames(trees.over.time) <- paste("TSp",c(1:tree.species),sep=".") #Tree species names
    mutualists.over.time <- array(numeric(0),c(m,mutualist.species.per.tree*tree.species,time.steps)) #Array to store mutualist composition through time
    pathogens.over.time <- array(numeric(0),c(m,pathogen.species.per.tree*tree.species,time.steps)) #Array to store pathogen composition through time
    survival.over.time <- array(numeric(0),c(m,tree.species,time.steps)) #Array to store pathogen composition through time
    colnames(mutualists.over.time) <- row.names(mutualist.effect.matrix) #Mutualist species names
    colnames(pathogens.over.time) <- row.names(pathogen.effect.matrix) #Pathogen species names
  }
  feedback.over.time <- mat.or.vec(tree.species,time.steps) #Mean feedback per species over time
  feedback.per.species <- numeric(tree.species) #Initialize feedback per species vector
  
  for (i in c(1:time.steps)) { #Keep looping through mortality/replacement/microbe population change for the predefined number of time steps
    if ((remove.mutualists==TRUE)&(i==(time.steps/2))) {
      mutualist.effect.matrix <- mat.or.vec(tree.species,tree.species)
    } #If the simulation is set to turn off the effect of mutualists, turn them off halfway through the simulation
    if ((remove.pathogens==TRUE)&(i==(time.steps/2))) {
      pathogen.effect.matrix <- mat.or.vec(tree.species,tree.species)
    } #If the simulation is set to turn off the effect of mutualists, turn them off halfway through the simulation
    Xm <- mutualists.adult + mutualist.propagule.matrix #All mutualists (including spores) across the simulation
    Xp <- pathogens.adult + pathogen.propagule.matrix #All pathogens (including spores) across the simulation
    total.mutualist.effects <- Xm%*%mutualist.effect.matrix #Multiply abundances by effect values per species to get total mutualist effects
    total.pathogen.effects <- Xp%*%pathogen.effect.matrix #Multiply abundances by effect values per species to get total pathogen effects
    surv.prob <- trial.function(mutualist.effect=total.mutualist.effects,pathogen.effect=total.pathogen.effects,g,h) #Survival probabilities of each species beneath each adult
    feedback.per.species <- numeric(tree.species) #Initialize the vector for average feedback per species
    tree.vec <- c(forest.matrix[,1],c(1:tree.species),c(1:tree.species)) #Make a temporary vec that includes all species (even if they've gone extinct in the simulation)
    surv.prob.temp <- rbind(surv.prob,matrix(rep(NA,2*tree.species^2),nrow=2*tree.species)) 
    surv.mat <- mat.or.vec(tree.species,tree.species) #Empty matrix to fill with reciprocal survival probabilities
    for (k in c(1:tree.species)) {
      surv.mat[k,] <- colMeans(surv.prob.temp[which(tree.vec==k),],na.rm=TRUE) 
    } #Calculate mean survival probability of each species beneath each other species
    feedback.mat <- mat.or.vec(tree.species,tree.species) #Empty matrix to fill with pairwise feedback values
    for (y in c(1:tree.species)) { 
      for (z in c(1:tree.species)) {
        feedback.mat[y,z] <- ((surv.mat[y,y]-surv.mat[y,z])/(surv.mat[y,y]+surv.mat[y,z])+(surv.mat[z,z]-surv.mat[z,y])/(surv.mat[z,z]+surv.mat[z,y]))/2
      }
    } #Calculate feedback per species pair
    feedback.mat[which(feedback.mat==0)] <- NA
    feedback.per.species <- rowSums(feedback.mat,na.rm=TRUE)/(tree.species-1) #Take average feedback score per species
    feedback.per.species[which(feedback.per.species==0)] <- NA
    feedback.over.time[,i] <- feedback.per.species #Add feedback values to the matrix tracking feedback over time
    
    dead.tree.IDs <- sample(forest.matrix[,2],round(m*mort),replace = FALSE) #Select random trees
    dead.trees.comm <- forest.matrix[dead.tree.IDs,] #Subset out data for dead trees
    
    survival.probabilities <- surv.prob[dead.tree.IDs,] #Survival probabilities of each species beneath the dead trees
    if ((eliminate.feedback==TRUE)&(i>=(time.steps/2))) {
      survival.probabilities <- (survival.probabilities+1)/(survival.probabilities+1)
    } #If the effect of feedback is set to stop, stop feedback halfway through the simulation
    recruitment.probabilities <- t(f.vals*t(survival.probabilities*seed.matrix[dead.tree.IDs,])) #Recruitment probabilities of each species beneath each of the dead trees
    
    if ((fix.feedback==FALSE)&(eliminate.feedback==FALSE)&(remove.pathogens==FALSE)) {
      if (length(which(is.na(recruitment.probabilities)))>0) { 
        exit.status <- 1 #Na.recruitment.probabilities
        break
      } #If there are NA recruitment probabilities, stop the simulation
      if (length(which(recruitment.probabilities<0))>0) { 
        exit.status <- 2 #Negative.recruitment.probabilities
        break
      } #If there are negative recruitment probabilities, stop the simulation
      if (length(which(recruitment.probabilities==0))==tree.species) { 
        exit.status <- 3 #No.nonzero.recruitment.probabilities
        break
      } #If there are no non-zero recruitment probabilities, stop the simulation
    }
    
    new.tree.IDs <- apply(recruitment.probabilities,1,function(x){ 
      which(rmultinom(1,1,x)==1)
    }) #Draw from the recruitment probabilities to select recruits at each position
    
    forest.matrix[dead.tree.IDs,1] <- new.tree.IDs #Update the species identities of the new recuits
    tree.matrix[dead.tree.IDs,] <- fact.matrix[forest.matrix[dead.tree.IDs,1],] #Update the tree community matrix
    if ((fix.feedback==FALSE)|(i<time.steps/2)) {
      dXm1.dt <- t(mutualist.growth.matrix[,forest.matrix[,1]])*Xm #Update mutualist abundances acording to competition equations part 1
      dXm2.dt <- dXm1.dt*(1-((1-t(mutualist.competition.matrix[,forest.matrix[,1]]))*Xm+rowSums(t(mutualist.competition.matrix[,forest.matrix[,1]])*Xm))) #Update mutualist abundances acording to competition equations part 2
      dXp1.dt <- t(pathogen.growth.matrix[,forest.matrix[,1]])*Xp #Update pathogen abundances acording to competition equations part 1
      dXp2.dt <- dXp1.dt*(1-((1-t(pathogen.competition.matrix[,forest.matrix[,1]]))*Xp+rowSums(t(pathogen.competition.matrix[,forest.matrix[,1]])*Xp))) #Update pathogen abundances acording to competition equations part 2
      mutualists.adult <- mutualists.adult + dXm2.dt
      mutualists.adult[which(mutualists.adult<0)] <- 0 #Prohibit negative abundances
      pathogens.adult <- pathogens.adult + dXp2.dt
      pathogens.adult[which(pathogens.adult<0)] <- 0 #Prohibit negative abundances
    }
    if ((fix.feedback==TRUE)&(i==time.steps/2)) {
      mutualists.fixed <- mat.or.vec(tree.species,tree.species)
      pathogens.fixed <- mat.or.vec(tree.species,tree.species)
      for (k in c(1:tree.species)) {
        mutualists.fixed[k,] <- colMeans(mutualists.adult[which(forest.matrix[,1]==k),])
        pathogens.fixed[k,] <- colMeans(pathogens.adult[which(forest.matrix[,1]==k),])
      }
    }
    if ((fix.feedback==TRUE)&(i>=time.steps/2)) {
      for (j in dead.tree.IDs) {
        mutualists.adult[j,] <- mutualists.fixed[forest.matrix[j,1],]
        pathogens.adult[j,] <- pathogens.fixed[forest.matrix[j,1],]
      }
    }
    
    mutualist.propagule.matrix <- gamma.m*matrix.pdf.mutualist%*%mutualists.adult #Create spore abundances for next round
    pathogen.propagule.matrix <- gamma.p*matrix.pdf.pathogen%*%pathogens.adult #Create spore abundances for next round
    seed.matrix <- matrix.pdf.plants%*%tree.matrix #Create seed dispersal kernels for each tree
    treeTab <- colSums(tree.matrix) #Tree abundance vector
    if (track.over.time==TRUE) {
      trees.over.time[,,i] <- tree.matrix
      mutualists.over.time[,,i] <- Xm
      pathogens.over.time[,,i] <- Xp
      survival.over.time[,,i] <- surv.prob/rowSums(surv.prob)
    } #Record tree community, mutualists, and pathogens for current time step
    
    if ((fix.feedback==FALSE)&(eliminate.feedback==FALSE)&(remove.pathogens==FALSE)) {
      if (length(which(treeTab>0))<tree.species) {
        exit.status <- 4 #Tree_species_loss
        break
      } #If tree species go extinct, exit simulation
    }
  } #Iterate through simulation
  returns <- list()
  returns$exit.status <- exit.status
  returns$forest.matrix <- forest.matrix
  returns$mutualist.matrix <- mutualist.effect.matrix
  returns$pathogen.matrix <- pathogen.effect.matrix
  returns$gamma.p <- gamma.p
  returns$gamma.m <- gamma.m
  returns$mutualists.adult <- mutualists.adult
  returns$pathogens.adult <- pathogens.adult
  returns$tree.community <- colSums(tree.matrix)
  returns$tree.species <- tree.species
  returns$g <- g
  returns$h <- h
  time.vec <- c(1:time.steps)
  if (track.over.time==TRUE) {
    if (subsample==TRUE) {
      returns$trees.over.time <- trees.over.time[,,time.vec%%10==0]
      returns$mutualists.over.time <- mutualists.over.time[,,time.vec%%10==0]
      returns$pathogens.over.time <- pathogens.over.time[,,time.vec%%10==0]
      returns$feedback.over.time <- feedback.over.time[,time.vec%%10==0]
      returns$survival.over.time <- survival.over.time[,,time.vec%%10==0]
    }
    else {
      returns$trees.over.time <- trees.over.time
      returns$mutualists.over.time <- mutualists.over.time
      returns$pathogens.over.time <- pathogens.over.time
      returns$feedback.over.time <- feedback.over.time
      returns$survival.over.time <- survival.over.time
    }
  } #Returns if abundances were tracked over time
  feedback.over.time[which(feedback.over.time==0)] <- -1
  returns$max.feedback.over.time <- max(feedback.over.time,na.rm=TRUE)
  returns$b.m <- b.m
  returns$b.p <- b.p
  returns$b.t <- b.t
  returns$f.vals <- f.vals
  returns$h <- h
  returns$position <- c(g,
                        h,
                        b.t,
                        min(f.vals),
                        s.m,
                        b.m,
                        alpha.m,
                        gamma.m,
                        r.m,
                        q.m,
                        c.m,
                        s.p,
                        b.p,
                        alpha.p,
                        gamma.p,
                        r.p,
                        q.p,
                        c.p)
  names(returns$position) <- c("g",
                               "h",
                               "b.t",
                               "f",
                               "s.m",
                               "b.m",
                               "alpha.m",
                               "gamma.m",
                               "r.m",
                               "q.m",
                               "c.m",
                               "s.p",
                               "b.p",
                               "alpha.p",
                               "gamma.p",
                               "r.p",
                               "q.p",
                               "c.p")
  return(returns)
}

#######################
#measure.diversity.abundance.correlation
#Returns results of correlation between mutualist/pathogen diversity and host abundance
#######################
measure.diversity.abundance.correlation <- function(modelOutput) {
  diversity.cor.mat <- mat.or.vec(length(modelOutput),2) #Empty matrix to fill with diversity-abundance correlation stats
  for (i in c(1:length(modelOutput))) {
    output <- modelOutput[[i]] #Take output from one simulation
    diversity.per.species.mutualists <- numeric(output$tree.species) #Initialize vector to fill with diversity of mutualists on each tree species
    diversity.per.species.pathogens <- numeric(output$tree.species) #Initialize vector to fill with diversity of pathogens on each tree species
    abundance <- numeric(output$tree.species) #Initialize vector to fill with tree abundances
    for (j in c(1:output$tree.species)) {
      indices <- which(output$forest.matrix[,1]==j) #Find which trees belong to a given species
      subs.mu <- output$mutualists.adult[indices,] #Take the subset of mutualists that occur on all individuals of a given tree species
      subs.pa <- output$pathogens.adult[indices,] #Take the subset of pathogens that occur on all individuals of a given tree species
      diversity.per.species.mutualists[j] <- mean(diversity(subs.mu),na.rm=TRUE) #Mean diversity of mutualists on a given species
      diversity.per.species.pathogens[j] <- mean(diversity(subs.pa),na.rm=TRUE) #Mean diversity of pathogens on a given species
      abundance[j] <- length(indices) #Count the number of trees for species j
    }
    diversity.cor.mat[i,1] <- cor.test(abundance,diversity.per.species.mutualists,na.action=na.rm)$estimate #Correlation coefficient for mutualist diversity-abunadnce correlation
    diversity.cor.mat[i,2] <- cor.test(abundance,diversity.per.species.pathogens,na.action=na.rm)$estimate #Correlation coefficient for pathogen diversity-abundance correlation
  }
  return(data.frame("Div.cor.m"=diversity.cor.mat[,1],"Div.cor.p"=diversity.cor.mat[,2]))
}

#######################
#measure.psf.strength
#Returns plant-soil feedback values of each tree
#######################
measure.PSF.strength <- function(modelOutput,
                                 microbe.dispersal.function) {
  forest.matrix <- modelOutput$forest.matrix #Take the forest matrix from model output
  m <- length(forest.matrix[,1]) #Number of trees in the simulation
  
  distance.vec <- c(0:(floor(m/2)),(floor(m/2)):1) #Distance values
  
  distance.pdf.mutualist <- microbe.dispersal.function(c(1:(floor(m/2)+1)),1,modelOutput$b.m) #Relative dispersal kernel of mutualists
  distance.pdf.pathogen <- microbe.dispersal.function(c(1:(floor(m/2)+1)),1,modelOutput$b.p) #Relative dispersal kernel of pathogens
  distance.matrix <- mat.or.vec(m,m) #Initialize distance matrix
  matrix.pdf.mutualist <- mat.or.vec(m,m) #Initialize tree-wise dispersal kernels of mutualists
  matrix.pdf.pathogen <- mat.or.vec(m,m) #Initialize tree-wise dispersal kernels of pathogens
  distance.matrix[1,] <- distance.vec #Assign first row of pairwise distance matrix
  
  for (i in c(1:(m-1))) {
    distance.matrix[i+1,] <- distance.vec[c((m-i+1):m,1:(m-i))]
  } #Fill pairwise distance matrix
  for (i in c(1:nrow(matrix.pdf.mutualist))) {
    matrix.pdf.mutualist[i,] <- distance.pdf.mutualist[distance.matrix[i,]+1]/sum(distance.pdf.mutualist[distance.matrix[i,]+1])
    matrix.pdf.pathogen[i,] <- distance.pdf.pathogen[distance.matrix[i,]+1]/sum(distance.pdf.pathogen[distance.matrix[i,]+1])
  } #Fill tree-wise dispersal kernels for pathogens and mutualists
  
  mutualist.propagule.matrix <- matrix.pdf.mutualist%*%modelOutput$mutualists.adult #Mutualist propagule abundances
  pathogen.propagule.matrix <- matrix.pdf.pathogen%*%modelOutput$pathogens.adult #Pathogen propagule abundances
  
  mutualist.effects <- (modelOutput$mutualists.adult + modelOutput$gamma.m*mutualist.propagule.matrix)%*%modelOutput$mutualist.matrix #Total effect of mutualists on each species at each position
  pathogen.effects <- (modelOutput$pathogens.adult + modelOutput$gamma.p*pathogen.propagule.matrix)%*%modelOutput$pathogen.matrix #Total effect of pathogens on each species at each position
  
  surv.prob <- trial.function(mutualist.effect=mutualist.effects,
                              pathogen.effect=pathogen.effects,
                              g=modelOutput$g,
                              h=modelOutput$h) #Calculate survival probability of each species at each location
  feedback.vec <- numeric(length(surv.prob[,1])) #Initialize vector of feedback strengths
  conspecific.vec <- numeric(length(trees))  #Initialize vector of survival beneath conspecifics
  
  trees <- modelOutput$tree.species #Number of tree species
  feedback.per.species <- numeric(trees) #Initialize vector to fill with feedback per species
  
  tree.vec <- c(forest.matrix[,1],c(1:trees),c(1:trees)) #Make a temporary vec that includes all species (even if they've gone extinct in the simulation)
  surv.prob.temp <- rbind(surv.prob,matrix(rep(NA,2*trees^2),nrow=2*trees)) 
  mutualist.effects.temp <- rbind(mutualist.effects,matrix(rep(NA,2*trees^2),nrow=2*trees)) 
  pathogen.effects.temp <- rbind(pathogen.effects,matrix(rep(NA,2*trees^2),nrow=2*trees))
  surv.mat <- mat.or.vec(trees,trees) #Initialize matrix to fill with reciprocal survival probabilities
  mutualist.mat <- mat.or.vec(trees,trees) #Initialize matrix to fill with mutualist abundances beneath each tree
  pathogen.mat <- mat.or.vec(trees,trees) #Initialize matrix to fill with pathogen abundances beneath each tree
  
  for (i in c(1:trees)) {
    surv.mat[i,] <- colMeans(surv.prob.temp[which(tree.vec==i),],na.rm=TRUE)
    mutualist.mat[i,] <- colMeans(mutualist.effects.temp[which(tree.vec==i),],na.rm=TRUE)
    pathogen.mat[i,] <- colMeans(pathogen.effects.temp[which(tree.vec==i),],na.rm=TRUE)
  } #Calculate mean survival probability, mutualist and pathogen abundances beneath each tree species
  print(surv.mat)
  feedback.mat <- mat.or.vec(trees,trees) #Initialize matrix to fill with feedback strengths between each species
  feedback.mat.mutualist <- mat.or.vec(trees,trees) #Initialize matrix to fill with mutualist feedbacks
  feedback.mat.pathogen <- mat.or.vec(trees,trees) #Initialize matrix to fill with pathogen feedbacks
  heterospecific.vec <- numeric(length(trees)) #Initialize vector to fill with survival probability of each tree species beneath heterospecifics
  conspecific.vec <- numeric(length(trees)) #Initialize vector to fill with survival probability of each tree species beneath conspecifics
  mutualist.consp.vec <- numeric(length(trees)) #Initialize vector to fill with abundances of mutualists beneath conspecifics
  mutualist.heterosp.vec <- numeric(length(trees)) #Initialize vector to fill with abundances of mutualists beneath heterospecifics
  pathogen.consp.vec <- numeric(length(trees)) #Initialize vector to fill with abundances of pathogens beneath conspecifics 
  pathogen.heterosp.vec <- numeric(length(trees)) #Initialize vector to fill with abundances of pathogens beneath heterospecifics
  for (i in c(1:trees)) {
    for (j in c(1:trees)) {
      feedback.mat[i,j] <- ((surv.mat[i,i]-surv.mat[i,j])/(surv.mat[i,i]+surv.mat[i,j])+(surv.mat[j,j]-surv.mat[j,i])/(surv.mat[j,j]+surv.mat[j,i]))/2
      feedback.mat.mutualist[i,j] <- ((mutualist.mat[i,i]-mutualist.mat[i,j])/(mutualist.mat[i,i]+mutualist.mat[i,j])+(mutualist.mat[j,j]-mutualist.mat[j,i])/(mutualist.mat[j,j]+mutualist.mat[j,i]))/2
      feedback.mat.pathogen[i,j] <- ((pathogen.mat[i,i]-pathogen.mat[i,j])/(pathogen.mat[i,i]+pathogen.mat[i,j])+(pathogen.mat[j,j]-pathogen.mat[j,i])/(pathogen.mat[j,j]+pathogen.mat[j,i]))/2
    } #Calculate pairwise feedbacks
    conspecific.vec[i] <- surv.mat[i,i]/sum(surv.mat[i,],na.rm=TRUE) #Relative survival probability of species i beneath conspecifics
    heterospecific.vec[i] <- mean(surv.mat[-i,i]/rowSums(surv.mat[-i,],na.rm=TRUE),na.rm=TRUE) #Relative survival probability of species i beneath heterospecifics
    mutualist.consp.vec[i] <- mutualist.mat[i,i]/sum(mutualist.mat[i,],na.rm=TRUE) #Relative abundances of mutualists of species i beneath conspecifics
    mutualist.heterosp.vec[i] <- mean(mutualist.mat[-i,i]/rowSums(mutualist.mat[-i,],na.rm=TRUE),na.rm=TRUE) #Relative abundances of mutualists of species i beneath heterospecifics
    pathogen.consp.vec[i] <- pathogen.mat[i,i]/sum(pathogen.mat[i,],na.rm=TRUE) #Relative abundances of pathogens of species i beneath conspecifics
    pathogen.heterosp.vec[i] <- mean(pathogen.mat[-i,i]/rowSums(pathogen.mat[-i,],na.rm=TRUE),na.rm=TRUE) #Relative abundances of pathogens of species i beneath heterospecifics
  } #Calculate pairwise feedbacks and average survival probabilities/mutalist/pathogen abundances beneath conspecifics and heterospecifics
  
  feedback.mat[which(feedback.mat==0)] <- NA 
  feedback.per.species <- rowSums(feedback.mat,na.rm=TRUE)/(trees-1) #Average feedback per species
  feedback.per.species[which(feedback.per.species==0)] <- NA
  mutualist.feedback.per.species <- rowSums(feedback.mat.mutualist,na.rm=TRUE)/(trees-1) #Average mutualist feedback per species
  pathogen.feedback.per.species <- rowSums(feedback.mat.pathogen,na.rm=TRUE)/(trees-1) #Average pathogen feedback per species
  feedbacks <- data.frame("TreeSpecies"=c(1:trees),
                          "Feedback"=feedback.per.species,
                          "Mu.feed"=mutualist.feedback.per.species,
                          "Pa.feed"=pathogen.feedback.per.species,
                          "Conspecific"=conspecific.vec,
                          "Heterospecific"=heterospecific.vec,
                          "Mu.Consp"=mutualist.consp.vec,
                          "Mu.heterosp"=mutualist.heterosp.vec,
                          "Pa.Consp"=pathogen.consp.vec,
                          "Pa.heterosp"=pathogen.heterosp.vec,
                          "Abundance"=table(tree.vec)-2)
  return(feedbacks)
}

#######################
#feedback.abundance.correlation
#Measures correlation between host abundance and negative feedback
#######################
feedback.abundance.correlation <- function(simulationOutput) {#feedback.consp.hetero refers to choice between feedback score, conspecific survival, or heterospecific survival probabilities
  nrows <- nrow(simulationOutput)
  variance.vec <- numeric(nrows)
  coefficient.feedback.vec <- numeric(nrows)
  coefficient.feedback.mutualist.vec <- numeric(nrows)
  coefficient.feedback.pathogen.vec <- numeric(nrows)
  coefficient.conspecific.vec <- numeric(nrows)
  coefficient.heterospecific.vec <- numeric(nrows)
  coefficient.consp.mutualist.vec <- numeric(nrows)
  coefficient.heterosp.mutualist.vec <- numeric(nrows)
  coefficient.consp.pathogen.vec <- numeric(nrows)
  coefficient.heterosp.pathogen.vec <- numeric(nrows)
  max.strength.vec <- numeric(nrows)
  min.strength.vec <- numeric(nrows)
  median.strength.vec <- numeric(nrows)
  mean.strength.vec <- numeric(nrows)
  
  feedback.matrix <- simulationOutput[,grepl("Feedback",colnames(simulationOutput))]
  feedback.matrix.mutualist <- simulationOutput[,grepl("Mu.feed",colnames(simulationOutput))]
  feedback.matrix.pathogen <- simulationOutput[,grepl("Pa.feed",colnames(simulationOutput))]
  conspecific.matrix <- simulationOutput[,grepl("Conspecific",colnames(simulationOutput))]
  heterospecific.matrix <- simulationOutput[,grepl("Heterospecific",colnames(simulationOutput))]
  mutualist.consp.matrix <- simulationOutput[,grepl("Mu.Consp",colnames(simulationOutput))]
  mutualist.heterosp.matrix <- simulationOutput[,grepl("Mu.heterosp",colnames(simulationOutput))]
  pathogen.consp.matrix <- simulationOutput[,grepl("Pa.Consp",colnames(simulationOutput))]
  pathogen.heterosp.matrix <- simulationOutput[,grepl("Pa.heterosp",colnames(simulationOutput))]
  abundance.matrix <- simulationOutput[,grepl("Abundance",colnames(simulationOutput))]
  for (i in c(1:nrow(feedback.matrix))) {
    abundance.vec <- as.numeric(abundance.matrix[i,])
    feedback.vec <- as.numeric(feedback.matrix[i,])
    feedback.vec.mutualist <- as.numeric(feedback.matrix.mutualist[i,])
    feedback.vec.pathogen <- as.numeric(feedback.matrix.pathogen[i,])
    conspecific.vec <- as.numeric(conspecific.matrix[i,])
    heterospecific.vec <- as.numeric(heterospecific.matrix[i,])
    mutualist.consp.vec <- as.numeric(mutualist.consp.matrix[i,])
    mutualist.heterosp.vec <- as.numeric(mutualist.heterosp.matrix[i,])
    pathogen.consp.vec <- as.numeric(pathogen.consp.matrix[i,])
    pathogen.heterosp.vec <- as.numeric(pathogen.heterosp.matrix[i,])
    max.strength.vec[i] <- max(feedback.vec,na.rm=TRUE)
    min.strength.vec[i] <- min(feedback.vec,na.rm=TRUE)
    median.strength.vec[i] <- median(feedback.vec,na.rm=TRUE)
    mean.strength.vec[i] <- mean(feedback.vec,na.rm=TRUE)
    variance.vec[i] <- sd(feedback.vec,na.rm=TRUE)
    
    if (length(which((feedback.vec!=0)&(!is.na(feedback.vec)))>=3)) {
      model.feedback <- cor.test(abundance.vec,feedback.vec)
      coefficient.feedback.vec[i] <- model.feedback$estimate
    }
    if (length(which((conspecific.vec!=0)&(!is.na(conspecific.vec)))>=3)) {
      model.conspecific <- cor.test(abundance.vec,conspecific.vec)
      coefficient.conspecific.vec[i] <- model.conspecific$estimate
    }
    if (length(which((heterospecific.vec!=0)&(!is.na(heterospecific.vec)))>=3)) {
      model.heterospecific <- cor.test(abundance.vec,heterospecific.vec)
      coefficient.heterospecific.vec[i] <- model.heterospecific$estimate
    }
    if (length(which((feedback.vec.mutualist!=0)&(!is.na(feedback.vec.mutualist)))>=3)) {
      model.feedback.mutualist <- cor.test(abundance.vec,feedback.vec.mutualist)
      coefficient.feedback.mutualist.vec[i] <- model.feedback.mutualist$estimate
      
    }
    if (length(which((mutualist.consp.vec!=0)&(!is.na(mutualist.consp.vec)))>=3)) {
      model.conspecific.mutualist <- cor.test(abundance.vec,mutualist.consp.vec)
      coefficient.consp.mutualist.vec[i] <- model.conspecific.mutualist$estimate
      
    }
    if (length(which((mutualist.heterosp.vec!=0)&(!is.na(mutualist.heterosp.vec)))>=3)) {
      model.heterospecific.mutualist <- cor.test(abundance.vec,mutualist.heterosp.vec)
      coefficient.heterosp.mutualist.vec[i] <- model.heterospecific.mutualist$estimate
      
    }
    if (length(which((feedback.vec.pathogen!=0)&(!is.na(feedback.vec.pathogen)))>=3)) {
      model.feedback.pathogen <- cor.test(abundance.vec,feedback.vec.pathogen)
      coefficient.feedback.pathogen.vec[i] <- model.feedback.pathogen$estimate
    }
    if (length(which((pathogen.consp.vec!=0)&(!is.na(pathogen.consp.vec)))>=3)) {
      model.conspecific.pathogen <- cor.test(abundance.vec,pathogen.consp.vec)
      coefficient.consp.pathogen.vec[i] <- model.conspecific.pathogen$estimate
    }
    if (length(which((pathogen.heterosp.vec!=0)&(!is.na(pathogen.heterosp.vec)))>=3)) {
      model.heterospecific.pathogen <- cor.test(abundance.vec,pathogen.heterosp.vec)
      coefficient.heterosp.pathogen.vec[i] <- model.heterospecific.pathogen$estimate
    }
  }
  output.data.frame <- data.frame("coefficient.feedback" = coefficient.feedback.vec,
                                  "coefficient.feedback.mutualist" = coefficient.feedback.mutualist.vec,
                                  "coefficient.feedback.pathogen" = coefficient.feedback.pathogen.vec,
                                  "coefficient.conspecific" = coefficient.conspecific.vec,
                                  "coefficient.heterospecific" = coefficient.heterospecific.vec,
                                  "coefficient.consp.mutualist" = coefficient.consp.mutualist.vec,
                                  "coefficient.heterosp.mutualist" = coefficient.heterosp.mutualist.vec,
                                  "coefficient.consp.pathogen" = coefficient.consp.pathogen.vec,
                                  "coefficient.heterosp.pathogen" = coefficient.heterosp.pathogen.vec,
                                  "max.strength" = max.strength.vec,
                                  "min.strength" = min.strength.vec,
                                  "mean.strength" = mean.strength.vec,
                                  "median.strength" = median.strength.vec,
                                  "sd.strength" = variance.vec)
  print(output.data.frame)
  output.data.frame
}

#######################
#nearest.neighbor.heterospecific
#Calculate the proportion of nearest neighbors that are conspecific vs heterospecific
#######################
nearest.neighbor.heterospecific <- function(forestMatrix) {
  dist.mat <- as.matrix(dist(forestMatrix$X,diag=TRUE,upper=TRUE))
  nearest.neighbor.indices <- unlist(apply(dist.mat,1,FUN=function(x) {which(x==sort(x)[2])[1]}))
  is.heterospecific <- forestMatrix$TreeSpecies!=forestMatrix$TreeSpecies[nearest.neighbor.indices]
  return(length(which(is.heterospecific==TRUE))/(length(is.heterospecific)))
}

#########################################################################################################################################################
optimize_swarm_azure <- function(continue=FALSE, #Wethere this run is a continuation of a previous set
                                 previous.output=NULL, #Output from previous set of runs
                                 subset.vector=NULL, #Indices of previous positions to use from previous set of runs
                                 swarm.size, #Number of independent simulations to run in each batch
                                 bounds.mat, #Value ranges for each of the parameters
                                 trees.in.forest, #Number of trees (cells) in the simulation
                                 mort, #Number of trees to replace in each time step in the simulation
                                 tree.species, #Number of tree species
                                 mortality.replacement.steps, #Number of time steps per simulation
                                 w, #Velocity weighting for particle swarm optimization
                                 c, #Distance weighting for particle swarm optimization
                                 nreplace, #Number of particles to reassign to a random position in each iteration of PSO
                                 maxRuns, #Maximum number of runs
                                 filePath, #Where to save run output
                                 azure, #TRUE or FALSE - whether to run the simulation using doAzureParallel
                                 #######################################################
                                 #Functions
                                 #See definitions
                                 #######################################################
                                 simulation.function,
                                 initiate.forest.matrix,
                                 initiate.fungal.matrix,
                                 mutualist.effect.function.consp,
                                 mutualist.effect.function.heterosp,
                                 pathogen.effect.function.consp,
                                 pathogen.effect.function.heterosp,
                                 trial.function,
                                 measure.PSF.strength,
                                 feedback.abundance.correlation,
                                 dispersal.function,
                                 microbe.dispersal.function
) {
  prop.soln.achieved <- numeric(0) #Initialize proportion of simulations achieving desired outcome
  colnames <- c("g",
                "h",
                "b.t",
                "f",
                "s.m",
                "b.m",
                "alpha.m",
                "gamma.m",
                "r.m",
                "q.m",
                "c.m",
                "s.p",
                "b.p",
                "alpha.p",
                "gamma.p",
                "r.p",
                "q.p",
                "c.p")
  if (continue==FALSE) {
    position.array <- array(numeric(),c(swarm.size,length(bounds.mat[,1]),1)) #Initialize array for particle positions
    random.particle.positions <- mat.or.vec(swarm.size,length(bounds.mat[,1])) #Initizalie matrix of particle positions
    for (d in c(1:length(bounds.mat[,1]))) {
      random.particle.positions[,d] <- runif(swarm.size,min=bounds.mat[d,1],max=bounds.mat[d,2])
    } #Select particle positions from random uniform distributions
    position.array[,,1] <- random.particle.positions #Assign random positions to first slot of positions array
    colnames(position.array) <- colnames
    explored.positions <- position.array[,,1] #Matrix of positions for which simulations have been run (add to this after each round of simulations)
    pareto.response.repository <- mat.or.vec(2,2) #Initialize matrix to store responses of best solutions (can be pareto front if optimizing for multiple conditions)
    pareto.position.repository <- mat.or.vec(2,nrow(bounds.mat)) #Initialize matrix to store positions of best solutions
    velocity.mat <- mat.or.vec(swarm.size,nrow(bounds.mat)) #Initizalize velocity matrix to store particle velocities for PSO
    for (d in c(1:length(bounds.mat[,1]))) {
      jitter <- runif(length(velocity.mat[,1]),min=-1,max=1)
      velocity.mat[,d] <- jitter*(bounds.mat[d,2]-bounds.mat[d,1])
    } #Generate random velocities for each particle
    t=1 #Set time (number of rounds of simulation batches) to 1
  } #If starting a new set of runs
  else {
    position.array <- array(numeric(0),c(swarm.size,length(bounds.mat[,1]),(previous.output$number.of.steps+1))) #Positions from last batch of previous runs
    colnames(position.array) <- colnames
    explored.positions <- as.matrix(previous.output$explored.positions[subset.vector,]) #All positions used for runs in previous set
    explored.responses <- as.matrix(previous.output$explored.responses[subset.vector,]) #Responses from previous set of runs
    pareto.response.repository <- previous.output$pareto.best.responses #Pareto responses from previous runs
    pareto.position.repository <- previous.output$pareto.best.positions #Pareto positions from previous runs
    velocity.mat <- mat.or.vec(length(position.array[,1,1]),nrow(bounds.mat)) #Reassign velocities to random velocities
    for (d in c(1:length(bounds.mat[,1]))) {
      jitter <- runif(length(velocity.mat[,1]),min=-1,max=1)
      velocity.mat[,d] <- jitter*(bounds.mat[d,2]-bounds.mat[d,1])
    } #Randomly assign particle velocities
    t <- previous.output$number.of.steps+1 
    new.positions <- mat.or.vec(swarm.size,length(bounds.mat[,1])) #Create new matrix of random position values (note: continuing runs does not work mid-optimization)
    for (d in c(1:length(bounds.mat[,1]))) {
      new.positions[,d] <- runif(swarm.size,min=bounds.mat[d,1],max=bounds.mat[d,2])
    } #Assign random position values
    for (i in c(1:swarm.size)) {
      position.array[i,,t] <- new.positions[i,]
    } #Assign to position to position array
  } #If continuing a batch of runs
  while(t<maxRuns) {
    if ((t>1)|(continue==TRUE)) {
      explored.positions <- rbind(explored.positions,position.array[,,t])
    } #Add positions from last batch onto list of explored positions
    f.vals <- sapply(position.array[,4,t],function(x) {(c(1:tree.species-1)/(tree.species-1))*(1-x)+x}) #Fitness values for all of the plant species
    if (azure==TRUE) {
      opts <- list(enableCloudCombine = TRUE)
    } #Define azure options
    if (azure==FALSE) {
      pb <- txtProgressBar(max = swarm.size, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
    } #Define options to create progress bar for runs if running simulations locally
    ptm <- proc.time() #Start the clock
    simulation.results <- foreach(g=position.array[,1,t],
                                  h=position.array[,2,t],
                                  b.t=position.array[,3,t],
                                  s.m = position.array[,5,t],
                                  b.m = position.array[,6,t],
                                  alpha.m=position.array[,7,t],
                                  gamma.m=position.array[,8,t],
                                  r.m=position.array[,9,t],
                                  q.m=position.array[,10,t],
                                  c.m=position.array[,11,t],
                                  s.p = position.array[,12,t],
                                  b.p = position.array[,13,t],
                                  alpha.p=position.array[,14,t],
                                  gamma.p =position.array[,15,t],
                                  r.p = position.array[,16,t],
                                  q.p = position.array[,17,t],
                                  c.p = position.array[,18,t],
                                  index=c(1:swarm.size),
                                  .options.azure = opts,.options.snow = opts) %dopar% #.packages=c("poweRlaw","vegan"), 
      (
        simulation.function(m=trees.in.forest,
                            mort=mort,
                            tree.species = tree.species,
                            mutualist.species.per.tree = 1,
                            pathogen.species.per.tree = 1,
                            time.steps=mortality.replacement.steps, #Number of mortality/replacement events
                            mutualist.effect.function.consp=mutualist.effect.function.consp, #Relationship between host affinity and benefit to the host
                            mutualist.effect.function.heterosp=mutualist.effect.function.heterosp,
                            pathogen.effect.function.consp=pathogen.effect.function.consp, #Relationship between host affinity and impact on the host
                            pathogen.effect.function.heterosp=pathogen.effect.function.heterosp,
                            g=g,
                            h=h,
                            s.p=s.p,
                            s.m=s.m,
                            b.p=b.p,
                            b.m=b.m,
                            b.t=b.t,
                            f.vals= rev(f.vals[,index]), #Probability that a seed will undergo trial if selected (per species)
                            index=index,
                            gamma.m=gamma.m, #Fecundity of mutualists
                            gamma.p=gamma.p, #Fecundity of pathogens
                            r.m=r.m, #Intrinsic growth rate of mutualists
                            r.p=r.p,
                            q.m=q.m,
                            q.p=q.p,
                            c.m=c.m,
                            c.p=c.p,
                            alpha.m=alpha.m,
                            alpha.p=alpha.p,
                            ################################ functions
                            initiate.forest.matrix=initiate.forest.matrix,
                            initiate.fungal.matrix=initiate.fungal.matrix,
                            trial.function=trial.function,
                            dispersal.function=dispersal.function,
                            microbe.dispersal.function=microbe.dispersal.function,
                            track.over.time=FALSE))
    
    cat("Simulation run time: ", (proc.time()-ptm)[3]/60, " min\n") #Print simulation run time
    ############################################################################
    #Evaluate output
    ############################################################################
    PSF.strength <- lapply(simulation.results, FUN = measure.PSF.strength, microbe.dispersal.function=microbe.dispersal.function) #Measure PSF strength for each tree species in each simulation run
    
    PSF.strength <- data.frame("Feedback"=matrix(unlist(lapply(PSF.strength,function(x) {x$Feedback})),ncol=tree.species,byrow = TRUE),
                               "Mu.feed"=matrix(unlist(lapply(PSF.strength,function(x) {x$Mu.feed})),ncol=tree.species,byrow = TRUE),
                               "Pa.feed"=matrix(unlist(lapply(PSF.strength,function(x) {x$Pa.feed})),ncol=tree.species,byrow = TRUE),
                               "Conspecific"=matrix(unlist(lapply(PSF.strength,function(x) {x$Conspecific})),ncol=tree.species,byrow=TRUE),
                               "Heterospecific"=matrix(unlist(lapply(PSF.strength,function(x) {x$Heterospecific})),ncol=tree.species,byrow=TRUE),
                               "Mu.Consp"=matrix(unlist(lapply(PSF.strength,function(x) {x$Mu.Consp})),ncol=tree.species,byrow=TRUE),
                               "Mu.heterosp"=matrix(unlist(lapply(PSF.strength,function(x) {x$Mu.heterosp})),ncol=tree.species,byrow=TRUE),
                               "Pa.Consp"=matrix(unlist(lapply(PSF.strength,function(x) {x$Pa.Consp})),ncol=tree.species,byrow=TRUE),
                               "Pa.heterosp"=matrix(unlist(lapply(PSF.strength,function(x) {x$Pa.heterosp})),ncol=tree.species,byrow=TRUE),
                               "Abundance"=matrix(unlist(lapply(PSF.strength,function(x) {x$Abundance.Freq})),ncol=tree.species,byrow=TRUE)) #Convert PSF output to data frame
    feedback.correlation <- feedback.abundance.correlation(PSF.strength) #Calculate feedback-abundance correlations
    correl.feedback <- feedback.correlation$coefficient.feedback 
    correl.feedback[which(is.na(correl.feedback))] <- 0 #Assign value of 0 to NA values (careful to not use stats in which 0 is a meaningful value)
    correl.feedback.mutualist <- feedback.correlation$coefficient.feedback.mutualist
    correl.feedback.mutualist[which(is.na(correl.feedback.mutualist))] <- 0
    correl.feedback.pathogen <- feedback.correlation$coefficient.feedback.pathogen
    correl.feedback.pathogen[which(is.na(correl.feedback.pathogen))] <- 0
    correl.conspecific <- feedback.correlation$coefficient.conspecific
    correl.conspecific[which(is.na(correl.conspecific))] <- 0
    correl.heterospecific <- feedback.correlation$coefficient.heterospecific
    correl.heterospecific[which(is.na(correl.heterospecific))] <- 0
    max.feedback.strength <- feedback.correlation$max.strength
    max.feedback.strength[which(is.na(max.feedback.strength))] <- 0
    mean.feedback.strength <- feedback.correlation$mean.strength
    mean.feedback.strength[which(is.na(mean.feedback.strength))] <- 0
    sd.plant.abundance <- unlist(lapply(simulation.results,function(x) {sd(x$tree.community)}))
    heterospecific.neighbor.proportion <- unlist(lapply(simulation.results,function(x) {nearest.neighbor.heterospecific(x$forest.matrix)}))
    plant.richness.vec <- unlist(lapply(simulation.results,function(x) {specnumber(x$tree.community)}))
    plant.diversity.vec <- unlist(lapply(simulation.results,function(x) {diversity(x$tree.community)}))
    plant.abundance.variation <- unlist(lapply(simulation.results,function(x) {sd(x$tree.community)}))
    mut.richness.vec <- unlist(lapply(simulation.results,function(x) {specnumber(colSums(x$mutualists.adult))}))
    mut.diversity.vec <- unlist(lapply(simulation.results,function(x) {diversity(colSums(x$mutualists.adult))}))
    pat.richness.vec <- unlist(lapply(simulation.results,function(x) {specnumber(colSums(x$pathogens.adult))}))
    pat.diversity.vec <- unlist(lapply(simulation.results,function(x) {diversity(colSums(x$pathogens.adult))}))
    exit.status.vec <- unlist(lapply(simulation.results,function(x) {x$exit.status}))
    sd.feedback.strength <- feedback.correlation$sd.strength
    sd.feedback.strength[which(is.na(sd.feedback.strength))] <- 0
    diversity.abundance.cor <- measure.diversity.abundance.correlation(simulation.results)
    max.feedback.over.time <- unlist(lapply(simulation.results,function(x) {x$max.feedback.over.time}))
    #random.seed <- unlist(lapply(simulation.results,function(x) {x$random.seed}))
    
    response.matrix <- cbind(plant.diversity.vec,
                             plant.richness.vec,
                             mut.richness.vec,
                             pat.richness.vec,
                             mut.diversity.vec,
                             pat.diversity.vec,
                             max.feedback.strength,
                             mean.feedback.strength,
                             correl.feedback,
                             correl.feedback.mutualist,
                             correl.feedback.pathogen,
                             correl.conspecific,
                             correl.heterospecific,
                             heterospecific.neighbor.proportion,
                             plant.abundance.variation,
                             sd.plant.abundance,
                             sd.feedback.strength,
                             diversity.abundance.cor$Div.cor.m,
                             diversity.abundance.cor$Div.cor.p,
                             max.feedback.over.time,
                             exit.status.vec)
    colnames(response.matrix) <- c("plant.shannon.diversity",
                                   "plant.richness",
                                   "mutualist.richness",
                                   "pathogen.richness",
                                   "mutualist.diversity",
                                   "pathogen.diversity",
                                   "max.feedback.strength",
                                   "mean.feedback.strength",
                                   "correl.feedback",
                                   "correl.feedback.mutualist",
                                   "correl.feedback.pathogen",
                                   "correl.conspecific",
                                   "correl.heterospecific",
                                   "heterospecific.neighbor.proportion",
                                   "plant.abundance.varition",
                                   "sd.plant.abundance",
                                   "sd.strength",
                                   "Div.cor.m",
                                   "Div.cor.p",
                                   "Max.feedback.over.time",
                                   "Exit.status")
    #print(response.matrix)
    if ((t==1)&(continue==FALSE)) {
      explored.responses <- response.matrix
    } #Assign response matrix to explored responses
    else {
      explored.responses <- rbind(explored.responses,response.matrix)
    } #Bind responses to explored responses
    ##############################################################################################################
    #Determine accuracy of random forest classifier
    ##############################################################################################################
    responseCorr <- as.factor((explored.responses[,9] > 0.8)*(explored.responses[,17] > 0.02)*(explored.responses[,2]==5)*(explored.responses[,3]==5)*(explored.responses[,7]<0)*(is.na(explored.responses[,21]))) #Categorize runs into those that did/did not maintain diversity and generate a positive feedback-abunndace correlation
    responseCorr[which(is.na(responseCorr))] <- 0 #Set NA values to 0 (i.e. FALSE)
    sampleSizesCorr <- rep(min(table(responseCorr)),length(unique(responseCorr))) #Maximum number of samples giving an equal number of 1s and 0s
    print(sampleSizesCorr)
    if (length(sampleSizesCorr)==2) {
      randFor <- randomForest(x=explored.positions, y=responseCorr, ntree = 1000,na.action=na.omit,sampsize = sampleSizesCorr,norm.votes=TRUE)
      votesCorr <- randFor$votes[,2]
    } #Fit random forest
    else {
      votesCorr <- 0
      print("Still building")
    } #If there aren't any runs that maintain diversity and create a positive correlation between abundance and diversity, don't fit a random forest
    
    if (length(sampleSizesCorr)==2) {
      cat("OOB error rate corr: ",randFor$err.rate[1000,1],"\n")
    } #Print error rate of random forest
    
    ############################################################################
    #Find pareto best responses... take worst values for each response
    ############################################################################
    pareto.indices <- which(votesCorr%in%sort(votesCorr,decreasing=TRUE)[c(1:10)]) #Take 10 best solutions based on votes from random forest
    
    pareto.response.repository <- explored.responses[pareto.indices,]
    pareto.position.repository <- explored.positions[pareto.indices,]
    if (length(sampleSizesCorr)==2) {
      print(cbind(pareto.response.repository[,c(2,3,9,10,12,17,21)],votesCorr[pareto.indices]))
      print(pareto.position.repository[,c(5,12)])
    } #If the random forest was fitted, print key results and positions
    
    ####################################################################################################################################################################################
    position.array <- abind(position.array,mat.or.vec(length(position.array[,1,1]),length(position.array[1,,1])),along=3) #Initialize new positions in position array
    if (length(pareto.indices)>1) {
      random.pareto.positions <- pareto.position.repository[sample(x=c(1:length(pareto.position.repository[,1])),size = length(position.array[,1,t+1]),replace=TRUE),]
    } #Select random particle for each particle to move toward
    else {
      random.pareto.positions <- matrix(pareto.position.repository,nrow=1) 
    }
    
    for (i in c(1:length(position.array[,1,t+1]))) {
      if (length(pareto.indices)>1) {
        pareto.position <- random.pareto.positions[i,]
      } 
      else {
        pareto.position <- random.pareto.positions
      }
      for (d in c(1:length(position.array[1,,t+1]))) {
        velocity.mat[i,d] <- w*velocity.mat[i,d]+c*runif(1,0,1)*(pareto.position[d]-position.array[i,d,t]) #Particle swarm equation        print(velocity.mat[i,d])
        position.array[i,d,t+1] <- position.array[i,d,t]+velocity.mat[i,d] #Update positions acording to new velocities
        if (position.array[i,d,t+1] < bounds.mat[d,1]) {
          position.array[i,d,t+1] <- bounds.mat[d,1]+abs(position.array[i,d,t+1]-bounds.mat[d,1])%%(bounds.mat[d,2]-bounds.mat[d,1])
        } #Enforce bounds
        if (position.array[i,d,t+1] > bounds.mat[d,2]) {
          position.array[i,d,t+1] <- bounds.mat[d,2]-abs(position.array[i,d,t+1]-bounds.mat[d,2])%%(bounds.mat[d,2]-bounds.mat[d,1])
        } #Enforce bounds
      }
    } #Particle swarm optimization algorithm
    
    ####################################################################################################################################################################################
    replace.indices <- sample(c(1:swarm.size),nreplace) #Randomly selected particles to replace
    number.replace <- nreplace
    if (length(sampleSizesCorr)==1) {
      replace.indices <- c(1:swarm.size)
      number.replace <- swarm.size
    }
    
    if (number.replace>0) {
      new.positions <- mat.or.vec(number.replace,length(bounds.mat[,1])) #Initialize new positions
      for (d in c(1:length(bounds.mat[,1]))) {
        new.positions[,d] <- runif(number.replace,min=bounds.mat[d,1],max=bounds.mat[d,2])
      } #Generate random positions
      cat("Number of particles replaced: ",number.replace,"\n")
      for (i in c(1:number.replace)) {
        position.array[replace.indices[i],,t+1] <- new.positions[i,]
      } #Assign random positions
    } #Move some (or all) particles to random positions if specified
    
    ##############################################################################################################
    #Create returns
    ##############################################################################################################
    
    returns <- list()
    position.colnames <- c("g",
                           "h",
                           "b.t",
                           "f",
                           "s.m",
                           "b.m",
                           "alpha.m",
                           "gamma.m",
                           "r.m",
                           "q.m",
                           "c.m",
                           "s.p",
                           "b.p",
                           "alpha.p",
                           "gamma.p",
                           "r.p",
                           "q.p",
                           "c.p")
    
    returns$pareto.best.positions <- pareto.position.repository
    if (is.matrix(returns$pareto.best.positions)) {
      colnames(returns$pareto.best.positions) <- position.colnames
    }
    returns$pareto.best.responses <- pareto.response.repository
    returns$explored.positions <- as.data.frame(explored.positions)
    colnames(returns$explored.positions) <- position.colnames
    returns$explored.responses <- as.data.frame(explored.responses)
    returns$number.of.steps <- t
    returns$position.array <- position.array
    returns$velocity.mat <- velocity.mat
    save(returns,file=filePath)
    cat("Step ",t," completed at ")
    cat(as.character(Sys.time()),"\n","==================================\n")
    t <- t+1
  }
  return(returns)
}
