#John Schroeder
#5-April-2020
#Spatially explicit plant-microbe interactions simulation plotting functions

#######################################################################################################################################################################################################################################################
#Functions to plot simulaton results
#######################################################################################################################################################################################################################################################

#######################
#calculate.microbes.through.space
#Returns data frame of average abundance values of microbes with increasing distance from focal hosts
#######################
calculate.microbes.through.space <- function(modelOutput, #List of simulation results
                                          mu.or.pa, #Mutualists or pathogens
                                          indices, #Vector of indices that correspond to the subset of modelOutput you would like to plot
                                          ncells #Number x axis positions defining spatial window for plotting
                                          ) {
  for (i in indices) {
    forest.matrix <- modelOutput[[i]]$forest.matrix
    m <- length(forest.matrix[,1])
    
    distance.vec <- c(0:(floor(m/2)),(floor(m/2)):1)
    
    distance.pdf.mutualist <- dpldis(c(1:(floor(m/2)+1)),1,modelOutput[[i]]$b.m)
    distance.pdf.pathogen <- dpldis(c(1:(floor(m/2)+1)),1,modelOutput[[i]]$b.p)
    distance.matrix <- mat.or.vec(m,m)
    matrix.pdf.mutualist <- mat.or.vec(m,m)
    matrix.pdf.pathogen <- mat.or.vec(m,m)
    distance.matrix[1,] <- distance.vec
    
    for (j in c(1:(m-1))) {
      distance.matrix[j+1,] <- distance.vec[c((m-j+1):m,1:(m-j))]
    }
    for (j in c(1:nrow(matrix.pdf.mutualist))) {
      matrix.pdf.mutualist[j,] <- distance.pdf.mutualist[distance.matrix[j,]+1]/sum(distance.pdf.mutualist[distance.matrix[j,]+1])
      matrix.pdf.pathogen[j,] <- distance.pdf.pathogen[distance.matrix[j,]+1]/sum(distance.pdf.pathogen[distance.matrix[j,]+1])
    }
    
    mutualist.propagule.matrix <- matrix.pdf.mutualist%*%modelOutput[[i]]$mutualists.adult
    pathogen.propagule.matrix <- matrix.pdf.pathogen%*%modelOutput[[i]]$pathogens.adult
    
    mutualist <- (modelOutput[[i]]$mutualists.adult + modelOutput[[i]]$gamma.m*mutualist.propagule.matrix)
    pathogen <- (modelOutput[[i]]$pathogens.adult + modelOutput[[i]]$gamma.p*pathogen.propagule.matrix)
    
    tree.community <- modelOutput[[i]]$tree.community
    dummy <- c(1,3,5)
    indicesd <- dummy
    if (mu.or.pa=="mu") {
      microbe <- mutualist
      ylab <- "Mutualist abundance"
    }
    else {
      microbe <- pathogen
      ylab <- "Pathogen abundance"
    }
    position <- forest.matrix$X
    dat <- data.frame("Tree.species"=as.factor(paste("TSp.",forest.matrix$TreeSpecies,sep="")),"Position" = position, "Abundance" = microbe[,indicesd])
    dat$host <- forest.matrix$TreeSpecies
    
    dist.mat <- mat.or.vec(length(forest.matrix$TreeSpecies),3)
    for (k in c(1:length(forest.matrix$TreeSpecies))) {
      dist <- abs(forest.matrix$X - k)
      matTemp <- cbind(forest.matrix$TreeSpecies,dist)
      for (j in c(1:3)) {
        subs <- matTemp[which(matTemp[,1]==dummy[j]),]
        dist.mat[k,j] <- min(subs[,2])
      }
    }
    colnames(dat) <- c("Tree.species","Position","Mic.1","Mic.2","Mic.3","host")
    colnames(dist.mat) <- c("Dist.1","Dist.2","Dist.3")
    
    dat <- cbind(dat,dist.mat)
    dat$host[which((dat$host != dummy[1])&(dat$host != dummy[2])&(dat$host != dummy[3]))] <- 10
    dat$host[which(dat$host==dummy[1])] <- 1
    dat$host[which(dat$host==dummy[2])] <- 2
    dat$host[which(dat$host==dummy[3])] <- 3
    
    df <- data.frame("Focal.sp"=numeric((ncells+1)*6),"Host"=numeric((ncells+1)*6),"Abundance"=numeric((ncells+1)*6),
                     "std.error.Abundance"=numeric((ncells+1)*6),"Distance"=numeric((ncells+1)*6))
    df$Adult <- rep(c(c("High",rep("Heterospecific",ncells)),
                      c("Medium",rep("Heterospecific",ncells)),
                      c("Low",rep("Heterospecific",ncells))),2)
    df$Focal.sp <- c(rep("High",2*(ncells+1)),rep("Medium",2*(ncells+1)),rep("Low",2*(ncells+1)))
    df$Host <- c(rep("High",(ncells+1)),rep("Heterospecific",(ncells+1)),
                 rep("Medium",(ncells+1)),rep("Heterospecific",(ncells+1)),
                 rep("Low",(ncells+1)),rep("Heterospecific",(ncells+1)))
    df$Distance <- rep(c(0:ncells),6)
    for (j in c(1:3)) {
      for (k in c(0:ncells)) {
        datSubs <- dat[which(dat[,(6+j)]==k),]
        df$Abundance[(j-1)*2*(ncells+1)+k+1] <- mean(datSubs[,(2+j)])
        df$std.error.Abundance[(j-1)*2*(ncells+1)+k+1] <- std.error(datSubs[,(2+j)])
        df$Abundance[(j-1)*2*(ncells+1)+(ncells+1)+k+1] <- mean(unlist(datSubs[,(2+c(1:3)[-j])]))
        df$std.error.Abundance[(j-1)*2*(ncells+1)+(ncells+1)+k+1] <- std.error(unlist(datSubs[,(2+c(1:3)[-j])]))
      }
    }
    df$run.number <- i
    
    if (i == min(indices)) {
      dff <- df
    }
    else {
      dff <- bind_rows(dff,df)
    }
    
  }

  df <- aggregate(Abundance ~ Focal.sp+Host+Distance+Adult,data=dff,FUN=mean)
  df$std.error.Abundance <- aggregate(Abundance ~ Focal.sp+Host+Distance+Adult,data=dff,FUN=std.error)$Abundance
  
  if (mu.or.pa=="mu") {
    df$Guild <- "Mutualist"
  }
  else {
    df$Guild <- "Pathogen"
  }
  return(df)
}


#######################
#plot.microbes.through.space
#Returns plot of average microbe abundances with increasing distance from focal host
#######################
plot.microbes.through.space <- function(mutualists, #Output from calculate.microbes.through.space
                                     pathogens, #Output from calculate.microbes.through.space
                                     survival #Output from plot.microbes.through.space
                                     ) {
  names(survival)[which(names(survival)=="Survival")] <- "Abundance"
  names(survival)[which(names(survival)=="std.error.survival")] <- "std.error.Abundance"
  survival$Guild <- "Survival"
  dat <- rbind(mutualists,pathogens,survival)
  
  colourCount <- 3
  
  colourCount <- 3
  getPalette1 = colorRampPalette(brewer.pal(3, "Dark2"))
  getPalette2 = colorRampPalette(brewer.pal(5, "PuOr"))
  colors1 <- c(getPalette1(1))
  colors2 <- c(getPalette2(2))
  colors <- c(colors2[1],colors2[2],colors1)
  
  dat <- subset(dat,Host!="Heterospecific")
  dat$Distance <- dat$Distance*6
  plot <- ggplot(data=dat,aes(x=Distance, y=Abundance, color=Host, linetype=Guild)) +
    theme_bw() +
    geom_line(size=1) +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    scale_linetype_manual(values=c('dotted','dashed','solid'))+
    geom_ribbon(data=subset(dat,Focal.sp=="High"),aes(ymin=Abundance - std.error.Abundance,ymax=Abundance + std.error.Abundance,fill=Host),alpha=0.5,color=NA) +
    geom_ribbon(data=subset(dat,Focal.sp=="Medium"),aes(ymin=Abundance - std.error.Abundance,ymax=Abundance + std.error.Abundance,fill=Host),alpha=0.5,color=NA) +
    geom_ribbon(data=subset(dat,Focal.sp=="Low"),aes(ymin=Abundance - std.error.Abundance,ymax=Abundance + std.error.Abundance,fill=Host),alpha=0.5,color=NA) +
    facet_wrap(~Guild,scales = "free_y") +
    labs(color="Fitness of\ntree species/\npreferred host") +
    guides(fill=FALSE,alpha=FALSE) +
    labs(x = "Distance (m) from focal tree",y="Abundance\nSurvival prob.")
  plot
}

#######################
#plot.survival.through.space
#Returns plot of average seedling survival with increasing distance from focal host
#######################
plot.survival.through.space <- function(modelOutput, #List of simulation results
                                        indices, #Vector of indices of model output
                                        ncells, #Distance (in cells) to show survival distribution
                                        fitness.dif #True or false. Whether to account for fitness differences in calculating survival probability
                                        ) {
  for (k in indices) {
    forest.matrix <- modelOutput[[k]]$forest.matrix
    m <- length(forest.matrix[,1])
    f.vals <- modelOutput[[k]]$f.vals
    distance.vec <- c(0:(floor(m/2)),(floor(m/2)):1)
    
    distance.pdf.mutualist <- dpldis(c(1:(floor(m/2)+1)),1,modelOutput[[k]]$b.m)
    distance.pdf.pathogen <- dpldis(c(1:(floor(m/2)+1)),1,modelOutput[[k]]$b.p)
    distance.matrix <- mat.or.vec(m,m)
    matrix.pdf.mutualist <- mat.or.vec(m,m)
    matrix.pdf.pathogen <- mat.or.vec(m,m)
    distance.matrix[1,] <- distance.vec
    for (i in c(1:(m-1))) {
      distance.matrix[i+1,] <- distance.vec[c((m-i+1):m,1:(m-i))]
    }
    for (i in c(1:nrow(matrix.pdf.mutualist))) {
      matrix.pdf.mutualist[i,] <- distance.pdf.mutualist[distance.matrix[i,]+1]/sum(distance.pdf.mutualist[distance.matrix[i,]+1])
      matrix.pdf.pathogen[i,] <- distance.pdf.pathogen[distance.matrix[i,]+1]/sum(distance.pdf.pathogen[distance.matrix[i,]+1])
    }
    mutualist.propagule.matrix <- matrix.pdf.mutualist%*%modelOutput[[k]]$mutualists.adult
    pathogen.propagule.matrix <- matrix.pdf.pathogen%*%modelOutput[[k]]$pathogens.adult
    
    mutualist.effects <- (modelOutput[[k]]$mutualists.adult + modelOutput[[k]]$gamma.m*mutualist.propagule.matrix)%*%modelOutput[[k]]$mutualist.matrix
    pathogen.effects <- (modelOutput[[k]]$pathogens.adult + modelOutput[[k]]$gamma.p*pathogen.propagule.matrix)%*%modelOutput[[k]]$pathogen.matrix
    
    surv.prob <- trial.function(mutualist.effect=mutualist.effects,
                                pathogen.effect=pathogen.effects,
                                g=modelOutput[[k]]$g,
                                h=modelOutput[[k]]$h)
    if (fitness.dif == TRUE) {
      surv.prob <- t(as.matrix(f.vals))[rep(1,length(surv.prob[,1])),]*surv.prob
    }
    surv.prob <- surv.prob/rowSums(surv.prob)
    surv.vec <- c(surv.prob)
    
    taxon <- rep(colnames(surv.prob),each=length(surv.prob[,1]))
    position <- rep(c(1:length(surv.prob[,1])),ncol(surv.prob))
    dat <- data.frame("Tree.species" = taxon,"Position" = position, "Survival" = surv.vec)
    dat$host <- rep(forest.matrix$TreeSpecies,length(mutualist.propagule.matrix[1,]))
    
    dist.mat <- mat.or.vec(length(forest.matrix$TreeSpecies),3)
    tree.community <- modelOutput[[k]]$tree.community
    dummy <-c(1,ceiling(length(unique(forest.matrix$TreeSpecies))/2),length(unique(forest.matrix$TreeSpecies)))
    dummy <- c(which(tree.community==sort(tree.community,decreasing=TRUE)[1]),
               which(tree.community==sort(tree.community,decreasing=TRUE)[dummy[2]]),
               which(tree.community==sort(tree.community,decreasing=TRUE)[dummy[3]]))
    dummy <- c(1,ceiling(length(unique(forest.matrix$TreeSpecies))/2),length(unique(forest.matrix$TreeSpecies)))#c(1,3,5)
    dummy <- c(1,3,5)
    
    for (i in c(1:length(forest.matrix$TreeSpecies))) {
      dist <- abs(forest.matrix$X - i)
      matTemp <- cbind(forest.matrix$TreeSpecies,dist)
      for (j in c(1:3)) {
        subs <- matTemp[which(matTemp[,1]==dummy[j]),]
        dist.mat[i,j] <- min(subs[,2])
      }
    }
    colnames(dist.mat) <- c("Dist.1","Dist.2","Dist.3")
    
    dat <- cbind(dat,dist.mat)
    
    dat$host[which((dat$host != dummy[1])&(dat$host != dummy[2])&(dat$host != dummy[3]))] <- 8
    dat$host[which(dat$host==dummy[1])] <- 1
    dat$host[which(dat$host==dummy[2])] <- 2
    dat$host[which(dat$host==dummy[3])] <- 3
    
    df <- data.frame("Focal.sp"=numeric((ncells+1)*6),"Host"=numeric((ncells+1)*6),"Survival"=numeric((ncells+1)*6),
                     "std.error.survival"=numeric((ncells+1)*6),"Distance"=numeric((ncells+1)*6))
    df$Focal.sp <- rep(c(rep("High",(ncells+1)),
                         rep("Medium",(ncells+1)),
                         rep("Low",(ncells+1))),2)
    df$Adult <- c(c("High",rep("Heterospecific",ncells)),
                  c("Medium",rep("Heterospecific",ncells)),
                  c("Low",rep("Heterospecific",ncells)))
    df$Host <- c(c(rep("High",(ncells+1)),rep("Medium",(ncells+1)),rep("Low",(ncells+1))),
                 rep("Heterospecific",(ncells+1)*3))
    df$Distance <- rep(c(0:ncells),6)
    
    for (j in c(1:3)) {
      dat1 <- subset(dat,Tree.species%in%paste("TSp.",dummy[j],sep=""))
      datConsp <- dat1[which(dat1[,(4+j)]==0),]
      df$Survival[((j-1)*(ncells+1)+1)] <- mean(datConsp[,3])
      df$std.error.survival[((j-1)*(ncells+1)+1)] <- std.error(datConsp[,3])
      for (i in c(1:ncells)) {
        datHet <- dat1[which(dat1[,(4+j)]==i),]
        df$Survival[((j-1)*(ncells+1)+i+1)] <- mean(datHet[,3])
        df$std.error.survival[((j-1)*(ncells+1)+i+1)] <- std.error(datHet[,3])
      }
    }
    
    for (j in c(1:3)) {
      dat1 <- subset(dat,Tree.species%in%paste("TSp.",c(dummy[c(1:3)[-j]]),sep=""))
      datConsp <- dat1[which(dat1[,(4+j)]==0),]
      df$Survival[3*(ncells+1) + ((j-1)*(ncells+1)+1)] <- mean(datConsp[,3])
      df$std.error.survival[3*(ncells+1) + ((j-1)*(ncells+1)+1)] <- std.error(datConsp[,3])
      for (i in c(1:ncells)) {
        datHet <- dat1[which(dat1[,(4+j)]==i),]
        df$Survival[3*(ncells+1) + ((j-1)*(ncells+1)+i+1)] <- mean(datHet[,3])
        df$std.error.survival[3*(ncells+1) + ((j-1)*(ncells+1)+i+1)] <- std.error(datHet[,3])
      }
    }
    
    df$run <- k
    if (k == min(indices)) {
      dff <- df
    }
    else {
      dff <- bind_rows(dff,df)
    }
  }
  df <- aggregate(Survival ~ Focal.sp+Host+Distance+Adult,data=dff,FUN=mean)
  df$std.error.survival <- aggregate(Survival ~ Focal.sp+Host+Distance+Adult,data=dff,FUN=std.error)$Survival

  colourCount <- 3
  getPalette1 = colorRampPalette(brewer.pal(3, "Dark2"))
  getPalette2 = colorRampPalette(brewer.pal(5, "PuOr"))
  colors1 <- c(getPalette1(1))
  colors2 <- c(getPalette2(2))
  colors <- c(colors2[1],colors2[2],colors1)
  
  df <- subset(df,Host!="Heterospecific")
  
  return(df)
}

#######################
#plot.feedback.per.species
#Returns plot of average feedback strength and abundance of species across simulation runs
#######################
plot.feedback.per.species <- function(cndd.strength,
                                      simulation.results,
                                      index
                                      ) {
  df <- list()
  tree.species <- simulation.results[[1]]$tree.species
  PSF.frame <- data.frame("Feedback"=matrix(unlist(lapply(cndd.strength,function(x) {x$Feedback})),ncol=tree.species,byrow = TRUE),
                          "Mu.feed"=matrix(unlist(lapply(cndd.strength,function(x) {x$Mu.feed})),ncol=tree.species,byrow = TRUE),
                          "Pa.feed"=matrix(unlist(lapply(cndd.strength,function(x) {x$Pa.feed})),ncol=tree.species,byrow = TRUE),
                          "Conspecific"=matrix(unlist(lapply(cndd.strength,function(x) {x$Conspecific})),ncol=tree.species,byrow=TRUE),
                          "Heterospecific"=matrix(unlist(lapply(cndd.strength,function(x) {x$Heterospecific})),ncol=tree.species,byrow=TRUE),
                          "Mu.Consp"=matrix(unlist(lapply(cndd.strength,function(x) {x$Mu.Consp})),ncol=tree.species,byrow=TRUE),
                          "Mu.heterosp"=matrix(unlist(lapply(cndd.strength,function(x) {x$Mu.heterosp})),ncol=tree.species,byrow=TRUE),
                          "Pa.Consp"=matrix(unlist(lapply(cndd.strength,function(x) {x$Pa.Consp})),ncol=tree.species,byrow=TRUE),
                          "Pa.heterosp"=matrix(unlist(lapply(cndd.strength,function(x) {x$Pa.heterosp})),ncol=tree.species,byrow=TRUE),
                          "Abundance"=matrix(unlist(lapply(cndd.strength,function(x) {x$Abundance.Freq/sum(x$Abundance.Freq)})),ncol=tree.species,byrow=TRUE))

  for (j in c(1:length(index))) {
    diversity.per.species.mutualists <- numeric(simulation.results[[j]]$tree.species)
    diversity.per.species.pathogens <- numeric(simulation.results[[j]]$tree.species)
    abundance <- numeric(simulation.results[[j]]$tree.species)
    for (i in c(1:simulation.results[[j]]$tree.species)) {
      indices <- which(simulation.results[[j]]$forest.matrix[,1]==i)
      subs.mu <- simulation.results[[j]]$mutualists.adult[indices,]
      subs.pa <- simulation.results[[j]]$pathogens.adult[indices,]
      diversity.per.species.mutualists[i] <- mean(diversity(subs.mu))
      diversity.per.species.pathogens[i] <- mean(diversity(subs.pa))
      abundance[i] <- length(indices)
    }
    newDf <- cndd.strength[[j]]
    df <- rbind(df,data.frame("Abundance" = rep(newDf$Abundance.Freq/sum(newDf$Abundance.Freq),6),
                              "Value"=c(newDf$Feedback,
                                        (newDf$Conspecific*newDf$Abundance.Freq + newDf$Heterospecific*(sum(newDf$Abundance.Freq)-newDf$Abundance.Freq))/sum(newDf$Abundance.Freq),
                                        newDf$Mu.feed,newDf$Pa.feed,
                                        diversity.per.species.mutualists,
                                        diversity.per.species.pathogens),

                              "Variable"=c(rep("Feedback",length(newDf[,1])),
                                           rep("Survival",length(newDf[,1])),#rep("Conspecific",length(newDf[,1])),rep("Heterospecific",length(newDf[,1])),
                                           rep("Feedback",length(newDf[,1])),rep("Feedback",length(newDf[,1])),
                                           rep("Diversity",length(newDf[,1])),rep("Diversity",length(newDf[,1]))),
                              "Type"=c(rep("Net PSF",length(newDf[,1])),
                                       rep("Survival",length(newDf[,1])),
                                       rep("Mutualist",length(newDf[,1])), rep("Pathogen",length(newDf[,1])),
                                       rep("Mutualist",length(newDf[,1])), rep("Pathogen",length(newDf[,1]))),
                              "Run"=as.factor(j)))
  }

  
  colourCount <- 3
  getPalette1 = colorRampPalette(brewer.pal(3, "Dark2"))
  getPalette2 = colorRampPalette(brewer.pal(5, "PuOr"))
  colors1 <- c(getPalette1(1))#,"gray")
  colors2 <- c(getPalette2(2))
  colors <- c(colors2[1],colors2[2],colors1,"Black")
  
  
  
  df.Feedback <- subset(df,Type=="Net PSF")
  
  df.Feedback$Fit <- rep(c("High","Other1","Medium","Other2","Low"),length(index))
  df.Feedback$Col <- rep(c("High","Other","Medium","Other","Low"),length(index))


  df.agg <- aggregate(Value ~ Col + Fit, data = df.Feedback, FUN=mean)
  df.agg$Abundance <- aggregate(Abundance ~ Col + Fit, data = df.Feedback, FUN=mean)$Abundance
  df.agg$Min.val <- aggregate(Value ~ Col + Fit, data = df.Feedback, FUN=function(x) {
    mean(x) - std.error(x)
  })$Value
  df.agg$Max.val <- aggregate(Value ~ Col + Fit, data = df.Feedback, FUN=function(x) {
    mean(x) + std.error(x)
  })$Value
  df.agg$Min.abund <- aggregate(Abundance ~ Col + Fit, data = df.Feedback, FUN=function(x) {
    mean(x) - std.error(x)
  })$Abundance
  df.agg$Max.abund <- aggregate(Abundance ~ Col + Fit, data = df.Feedback, FUN=function(x) {
    mean(x) + std.error(x)
  })$Abundance
  
  feedPlot <- ggplot(df.agg,aes(x=Abundance,y=Value,color=Col))+
    labs(x="Tree species rel. abund.",y="Feedback value (I)", color="Fitness of\ntree species /\npreferred host") +
    geom_errorbar(aes(ymin = Min.val,ymax = Max.val),width=0.005) + 
    geom_errorbarh(aes(xmin = Min.abund,xmax = Max.abund),height=0.005) +
    guides(fill=FALSE) +
    geom_point(size=2) +
    scale_color_manual(values=colors) +
    theme_bw()
  return(feedPlot)  
  
}

#######################
#calculate.microbes.over.time.single.tree
#Returns data frame with average temporal dynamics of microbial abundances on a single host through time
#######################
calculate.microbes.over.time.single.tree <- function(modelOutput, #List of simulation results
                                                     mu.or.pa, #Mutualists or pathogens "mu" or "pa"
                                                     time.steps, #Number of time steps to illustrate in plot
                                                     step.range, #Number of time steps over which to calculate average dynamics
                                                     indices #Vector indicating which of the model results to average
                                                     ) {
  for (x in indices) {
    if (mu.or.pa=="mu") {
      mic.over.time <- modelOutput[[x]]$mutualists.over.time
      ylab = "Mutualist abundance"
    }
    else {
      mic.over.time <- modelOutput[[x]]$pathogens.over.time
      ylab = "Pathogen abundance"
    }
    mic.over.time <- mic.over.time[,,step.range]
    
    trees.over.time <- modelOutput[[x]]$trees.over.time[,,step.range]
    time.array <- modelOutput[[x]]$trees.over.time[,,step.range]
    
    for (i in c(2:length(trees.over.time[1,1,]))) {
      time.array[,,i] <- time.array[,,i-1] + trees.over.time[,,i]
      time.array[,,i][which(trees.over.time[,,i-1]-trees.over.time[,,i]==1)] <- 0
    }
    
    time <- time.steps
    tree.community <- modelOutput[[x]]$tree.community
    forest.matrix <- modelOutput[[x]]$forest.matrix
    dummy <-c(1,ceiling(length(unique(forest.matrix$TreeSpecies))/2),length(unique(forest.matrix$TreeSpecies)))
    dummy <- c(which(tree.community==sort(tree.community,decreasing=TRUE)[1]),
               which(tree.community==sort(tree.community,decreasing=TRUE)[dummy[2]]),
               which(tree.community==sort(tree.community,decreasing=TRUE)[dummy[3]]))
    dummy <-c(1,ceiling(length(unique(forest.matrix$TreeSpecies))/2),length(unique(forest.matrix$TreeSpecies)))
    dummy <- c(1,3,5)
    
    abund.mat <- array(numeric(),c(3,3,time))
    for (i in c(1:time)) {
      for (j in c(1:3)) {
        for (k in c(1:3)) {
          abund.mat[j,k,i] <- mean(mic.over.time[,dummy[k],][which(time.array[,dummy[j],]==i)])
        }
      }
    }
    df <- data.frame("Focal.tree"=rep(rep(c("High","Medium","Low"),3),time),
                     "Preferred.host"= rep(rep(c("High","Medium","Low"),each=3),time),
                     "Abundance"=as.vector(abund.mat),
                     "Time"=rep(c(1:time),each=9))
    df$Time <- df$Time*5
    
    if (x == min(indices)) {
      dff <- df
    }
    else {
      dff <- bind_rows(dff,df)
    }
  }
  df <- aggregate(Abundance ~ Focal.tree +
                    Preferred.host +
                    Time,data=dff,FUN=mean)
  df$std.error <- aggregate(Abundance ~ Focal.tree +
                              Preferred.host +
                              Time,data=dff,FUN=std.error)$Abundance
  
  if (mu.or.pa =="mu") {
    df$Guild <- "Mutualists"
  }
  else {
    df$Guild <- "Pathogens"
  }
  return(df)
}

#######################
#plot.microbes.over.time
#Returns plot of output from calculate.microbes.over.time.single.tree
#######################
plot.microbes.over.time <- function(mutualists.over.time, #Output from calculate.microbes.over.time.single.tree 
                                    pathogens.over.time #Output from calculate.microbes.over.time.single.tree 
                                    ) {
  
  df <- rbind(mutualists.over.time,pathogens.over.time)
  colourCount <- 3
  getPalette1 = colorRampPalette(brewer.pal(3, "Dark2"))
  getPalette2 = colorRampPalette(brewer.pal(5, "PuOr"))
  colors1 <- c(getPalette1(1))
  colors2 <- c(getPalette2(2))
  colors <- c(colors2[1],colors2[2],colors1)
  
  df$Focal.tree = factor(df$Focal.tree, levels=c('Low','Medium','High'))
  
  plot <- ggplot(df,aes(x=Time,y=Abundance,color=Preferred.host,linetype=Guild)) +
    scale_linetype_manual(values=c("dotted","dashed")) +
    facet_wrap(~Focal.tree) +
    geom_line(size=1) +
    geom_ribbon(aes(ymin=Abundance-std.error,ymax=Abundance+std.error,fill=Preferred.host),alpha=0.3,color=NA) +
    theme_bw() +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    guides(fill=FALSE) +
    labs(color = "Fitness of\ntree species/\npreferred host") +
    labs(x = "Time since recruitment (years)",y="Abundance")
  plot
}
#######################
#survival.over.time.single.tree
#Calculates and plots average survival probability dynamics of seedlings beneath three focal tree species
#######################
survival.over.time.single.tree <- function(modelOutput, #List of simulation results
                                           time.steps, #Number of time steps to illustrate in plot
                                           step.range, #Number of time steps over which to calculate average dynamics
                                           indices #Vector indicating which of the model results to average
                                           ) {
  for (x in indices) {
    forest.matrix <- modelOutput[[x]]$forest.matrix
    surv.over.time <- modelOutput[[x]]$survival.over.time[,,step.range]
    trees.over.time <- modelOutput[[x]]$trees.over.time[,,step.range]
    time.array <- modelOutput[[x]]$trees.over.time[,,step.range]
    for (i in c(2:length(trees.over.time[1,1,]))) {
      time.array[,,i] <- time.array[,,i-1] + trees.over.time[,,i]
      time.array[,,i][which(trees.over.time[,,i-1]-trees.over.time[,,i]==1)] <- 0
    }
    
    time <- time.steps
    dummy <-c(1,ceiling(length(unique(forest.matrix$TreeSpecies))/2),length(unique(forest.matrix$TreeSpecies)))
    tree.community <- modelOutput[[x]]$tree.community
    dummy <- c(which(tree.community==sort(tree.community,decreasing=TRUE)[1]),
               which(tree.community==sort(tree.community,decreasing=TRUE)[dummy[2]]),
               which(tree.community==sort(tree.community,decreasing=TRUE)[dummy[3]]))
    dummy <- c(1,3,5)
    surv.mat <- array(numeric(),c(3,3,time))

    for (i in c(1:time)) {
      for (j in c(1:3)) {
        for (k in c(1:3)) {
          surv.mat[j,k,i] <- mean(surv.over.time[,dummy[k],][which(time.array[,dummy[j],]==i)])
          #se.mat[j,k,i] <- std.error(surv.over.time[,dummy[k],][which(time.array[,dummy[j],]==i)])
        }
      }
    }
    
    df <- data.frame("Focal.tree"=rep(rep(c("High","Medium","Low"),3),time),
                     "Preferred.host"= rep(rep(c("High","Medium","Low"),each=3),time),
                     "Survival"=as.vector(surv.mat),
                     "Time"=rep(rep(c(1:time),each=9)))
    
    if (x == min(indices)) {
      dff <- df
    }
    else {
      dff <- bind_rows(dff,df)
    }
    
  }
  
  df <- aggregate(Survival ~ Focal.tree +
                    Preferred.host +
                    Time,data=dff,FUN=mean)
  df$std.error <- aggregate(Survival ~ Focal.tree +
                              Preferred.host +
                              Time,data=dff,FUN=std.error)$Survival
  
  colourCount <- 3
  getPalette1 = colorRampPalette(brewer.pal(3, "Dark2"))
  getPalette2 = colorRampPalette(brewer.pal(5, "PuOr"))
  colors1 <- c(getPalette1(1))#,"gray")
  colors2 <- c(getPalette2(2))
  colors <- c(colors2[1],colors2[2],colors1)
  
  df$Focal.tree = factor(df$Focal.tree, levels=c('Low','Medium','High'))
  
  df$Time <- df$Time*5
  
  plot <- ggplot(df,aes(x=Time,y=Survival,color=Preferred.host)) +
    facet_wrap(~Focal.tree) +
    geom_line() +
    geom_ribbon(aes(ymin=Survival-std.error,ymax=Survival+std.error,fill=Preferred.host),alpha=0.3,color=NA) +
    theme_bw() +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    labs(x="Time since recruitment (years)", ylab="Rel. surv. prob.",fill="Fitness of\ntree species\npreferred host")
  plot
}

#######################
#plot.tree.abundances.over.time
#Plots mean +/- se abundance of three focal species through time
#######################
plot.tree.abundances.over.time <- function(model.output, #List of simulation results
                                           step.range, #Time step range over which to plot results
                                           indices, #Vector indicating which simulation results to average
                                           xlim #Another way to constrain x axis
                                           ) {
  for (i in indices) {
    forest.matrix <- model.output[[i]]$forest.matrix
    tree.community <- model.output[[i]]$tree.community
    dummy <- c(1,3,5)
    tree.abundance.matrix <- (colSums(model.output[[i]]$trees.over.time)/colSums(colSums(model.output[[i]]$trees.over.time)))[dummy,step.range]
    tree.species <- rep(rownames(tree.abundance.matrix),each = ncol(tree.abundance.matrix))
    tree.abundances <- as.vector(t(tree.abundance.matrix))
    time = rep(c(1:ncol(tree.abundance.matrix)),nrow(tree.abundance.matrix))
    dat <- data.frame("tree.species"=tree.species,
                      "Abundances"=tree.abundances,
                      "time"=time)
    dummynames <- paste("TSp.",dummy,sep="")
    dat$tree.species <- as.character(dat$tree.species)
    dat$tree.species[which(dat$tree.species==dummynames[1])] <- "High"
    dat$tree.species[which(dat$tree.species==dummynames[2])] <- "Medium"
    dat$tree.species[which(dat$tree.species==dummynames[3])] <- "Low"
    dat$tree.species <- as.factor(dat$tree.species)
    dat$run.number <- i
    if (i==min(indices)) {
      datt <- dat
    }
    else{
      datt <- rbind(datt,dat)
    }
  }
  dat <- aggregate(Abundances~tree.species+time,data=datt,FUN=mean)
  dat$std.error <- aggregate(Abundances~tree.species+time,data=datt,FUN=sd)$Abundances

  dat <- subset(dat,dat$time%%10==0)
  colourCount = length(unique(tree.species))
  getPalette = colorRampPalette(brewer.pal(colourCount, "Dark2"))
  colors <- getPalette(colourCount)
  colors <- colors[c(2,3,1)]

  dat$time <- dat$time*5
  dat$time <- dat$time/1000
  dat$std.error.min <- dat$Abundances-dat$std.error
  dat$std.error.min[which(dat$std.error.min<0)] <- 0
  dat$std.error.max <- dat$Abundances+dat$std.error
  dat$std.error.max[which(dat$std.error.max>1)] <- 1
  
  
  ggplot(data=dat,aes(x=time, y=Abundances, colour=tree.species))+
    theme_bw()+
    xlim(xlim) +    
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    geom_ribbon(aes(ymin=std.error.min,ymax=std.error.max,fill=tree.species),alpha=0.5,color=NA) +
    geom_line(size=0.3,alpha=0.8) +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    scale_linetype_manual(values=rep('solid',8)) +
    labs(x = "Time (K years)",y="Prop. Abund.",col="Tree fitness",title="Tree abundance over time")
}








