pacman::p_load("INLA","Matrix","reshape2","rstan","ggdist","tidyverse")

##### Helper functions needed to run the scripts ###############################

Sum_Stan <- function(samples, digits = 2, ...){
  rstan::summary(object = samples,...)$summary %>% 
    round(digits = digits)
}

#### Functions for creation of BYM2 Model ######################################
CARData4Stan <- function(NeighborhoodMatrix){ #Region in Stan Matrix
N <- nrow(NeighborhoodMatrix) #Amount of Regions
N_edges <- sum(NeighborhoodMatrix) #Amount of Edges
node1 <- vector(mode="numeric", length=N_edges)
node2 <- vector(mode="numeric", length=N_edges)
iEdge <- 1 #Helpvector
for (j in 1:nrow(NeighborhoodMatrix)) {
  NumEdge <- sum(NeighborhoodMatrix[j,]) #Sum Neighbors
  TypeEdge <- which(NeighborhoodMatrix[j,]!=0) #Check which Region is Neighbor 
  node1[iEdge:(iEdge+NumEdge-1)] <- rep(j,NumEdge) 
  node2[iEdge:(iEdge+NumEdge-1)] <- TypeEdge 
  iEdge <- iEdge+NumEdge #Update Help Vector
}

node1Names <- rownames(NeighborhoodMatrix)[node1] 
node2Names <- rownames(NeighborhoodMatrix)[node2] 

return (list("N"=N,"N_edges"=N_edges,"node1"=node1,"node2"=node2,
             "node1Names"=node1Names, node2Names=node2Names))
} 

#Compute Scaling Factor for BYM2 Model
#Code taken from Morris (2019) 
#https://mc-stan.org/users/documentation/case-studies/icar_stan.html
ScalingFacBYM2 <- function(Nodes, AdjMat){
  #Add a small jitter to the diagonal for numerical stability 
  Q <- Matrix::Diagonal(Nodes$N, Matrix::rowSums(AdjMat)) - AdjMat
  # Compute the diagonal elements of the covariance matrix subject to the 
  # constraint that the entries of the ICAR sum to zero.
  #See the inla.qinv function help for further details.
  Q_pert <-  Q + Matrix::Diagonal(Nodes$N) * max(diag(Q)) * sqrt(.Machine$double.eps)
  
  # Compute the diagonal elements of the covariance matrix subject to the 
  # constraint that the entries of the ICAR sum to zero.
  #See the inla.qinv function help for further details.
  Q_inv <-  INLA::inla.qinv(Q_pert, constr=list(A = matrix(1,1,Nodes$N),e=0))
  
  #Compute the geometric mean of the variances, which are on the diagonal of Q.inv
  scaling_factor <-  exp(mean(log(diag(Q_inv))))
  
  return(scaling_factor)
}


#Calculation of Cohort Index for Mortality data
CohortIndex <- function(agegroup, time, maxAge, M){
  k <- M*(maxAge-agegroup)+time
  return(k=k)
}


##### Functions for Population Forecasts #######################################
#Leslie Matrix
Leslie <- function(L, f, SRB, l0, sex="Female") { 
  n <-  length(L) #Amount of Age Groups
  M <-  matrix(0, nrow = n, ncol =  n) #Leslie Matrix
  
  # lower diagonal has survivorship ratios
  for (i in 1:(n-1)) {
    M[i+1,i] <- L[i+1]/L[i] #Person Years Lived
  }
  M[n,n-1] <- M[n,n] <- L[n]/(L[n-1] + L[n]) #Compare. Preston p. 128 (Tx/Tx-1)
  
  # first row has net maternity contributions
  if(sex=="Female"){
    k <- (1/(1+SRB))*(L[1]/5*l0)
  }else { #for Males
    k <- (SRB/(1+SRB))*(L[1]/5*l0)
  }
  for(i in 1:(n-1)) {
    if(f[i] != 0 | f[i+1] != 0) {#Check if fertility is non zero
      M[1,i] <- k*(f[i]+f[i+1]*(L[i+1]/L[i]))*5/2 
    } else {
      #do nothing
    }
  }
  return(M)
}

#Creation of Life Table 
LifeTableFun <- function(Age, MortalityRate, radix=100000, sex="female"){

  LifeTable <- data.frame("Age"=Age,
                          "mx"=MortalityRate)
  
  last <- nrow(LifeTable) #Last element
  
  LifeTable <- LifeTable %>% mutate("n" = c(diff(Age), 0), #n 
                                    "a"=n/2)
  
  LifeTable[last,"a"] <- 1/LifeTable$mx[last] #last element
  
  #Probability of dying
  LifeTable <- LifeTable %>% mutate("qx" = n * mx/(1 + (n - a) * mx))
  
  LifeTable[last,"qx"] <- 1
  
  LifeTable <- LifeTable %>% mutate(
    "px"=1-qx, #Prob of surviving
    "lx"=radix * cumprod( c(1, px[-last])), #Number of people left alive
    "dx" = c(-diff(lx), lx[last]), #Number of deaths
    "Lx"=n*(lx - dx) +  (dx * a), #Total Person-Years Lived
    "Tx" =  rev(cumsum(rev(Lx))), #Person Years of life Remaining
    "ex" = Tx/lx) #Life Expectation
  
  LifeTable[last,"Lx"] <-  LifeTable$lx[last]/LifeTable$mx[last]
  
  return(LifeTable) 
}

#### Functions to Create Population Pyramid ####################################
## Helper function for nice, symmetric axis labels 
pretty_symmetric <- function(range, n = 5){
  range_1 <- c(-range[1], range[2])
  range_2 <- c(range[1], -range[2])
  
  pretty_vec_1 <- pretty(range_1)
  pretty_vec_2 <- pretty(range_2)
  
  pretty(
    c(pretty_vec_1, pretty_vec_2), 
    n = n
  )
  
}

PopulationPyramid <- function(dataF = PopTrajectories_Female, 
                              dataM = PopTrajectories_Male, 
                              Region = 1, 
                              YearUpper = 2049){
  
  
  
  #Transform into long format and summarise
  LongDataF <- reshape2::melt(dataF, value.name = "Pop", 
                              varnames = c("age_group","Year",
                                           "RegionNumber","Iteration"),
                              as.is = TRUE) %>% 
    reframe(ggdist::median_qi(Pop, .width = c(0.8, 0.95)),
            .by = c(age_group, Year, RegionNumber)) %>% 
    mutate(Sex = "female")  %>% 
    mutate(".width" = factor(.width)) %>% 
    mutate("age_group" = factor(age_group, levels = unique(age_group)))
  
  LongDataM <- reshape2::melt(dataM, value.name = "Pop", 
                              varnames = c("age_group","Year","RegionNumber","Iteration"),
                              as.is = TRUE) %>% 
    reframe(ggdist::median_qi(Pop, .width = c(0.8, 0.95)),
            .by = c(age_group, Year, RegionNumber)) %>% 
    mutate("Sex" = "male") %>% 
    mutate(".width" = factor(.width)) %>% 
    mutate("age_group" = factor(age_group, levels = unique(age_group))) %>% 
    mutate(across(where(is.numeric), ~ .* -1)) #multiply with minus one for males
  
  #combine both data sets
  JoinedData <- 
    full_join(x = LongDataF, 
              y = LongDataM) %>% 
    filter(Year %in% c(2024, YearUpper)) %>% 
    filter(RegionNumber == unique(RegionNumber)[{Region}]) %>% 
    mutate(AgeID = match(age_group, unique(age_group))) #%>%  #ID for second axis

  pop_range <- range(c(JoinedData$ymin, JoinedData$ymax))
  pop_range_seq <- pretty_symmetric(pop_range, n = 7)
  
  alpha_var <- 0.2
  alpha_bar <- 0.8
  width_obs <- 0.3
  width_fc <- 0.9
  
  legend_data <- data.frame(
    age_group = 0,
    y = 0,
    legend_label = factor(c(paste(YearUpper, "Median"), "95% PI", "80% PI", "2024"),
                          levels = c(paste(YearUpper, "Median"), "95% PI", "80% PI", "2024"))
  )
  
  # --- Plot ---

    Plot <- 
      ggplot(data = JoinedData, aes(x = AgeID)) +
      # --- Main plotting layers (suppress all default legends) ---
      # --- Add PI_s of 2048
      geom_errorbar(data = JoinedData %>% filter(Year == YearUpper) %>% arrange(rev(.width)),
                    stat = "identity", position = "identity", 
                    width = 0, linewidth = 13.5, alpha = alpha_bar,
                    aes(ymin = ymin, ymax = ymax,
                        color = .width),
                    show.legend = FALSE) +
      # --- add current year of 2023
      geom_bar(data = JoinedData %>% filter(Year == 2024),
               stat = "identity", 
               aes(fill = Year, y = y),
               position = "identity", width = width_obs,
               alpha = alpha_var,
               show.legend = FALSE) +
      # --- a small bar on top of geom bar
      geom_errorbar(data = JoinedData %>% filter(Year == 2024),
                    aes(ymin = ymin, ymax = ymax, 
                        col = Year),
                    alpha = alpha_var, linewidth = 1.2, 
                    width = width_obs + 0.1,
                    show.legend = FALSE) +
      # --- add median of future forecasts
      geom_bar(data = JoinedData %>% filter(Year == YearUpper),
               stat = "identity", position = "identity",
               aes(col = Year, y = y), fill = NA, width = width_fc,
               linewidth = 0.3,
               show.legend = FALSE) +
      # --- add small bar on top of geom bar
      geom_errorbar(data = JoinedData %>% filter(Year == YearUpper),
                    aes(ymin = y, ymax = y, col = Year),
                    linewidth = 1.2, width = width_fc + 0.1,
                    show.legend = FALSE)+
      # --- add empty dummy data for joined legend #
      geom_line(data = legend_data, 
                aes(x = age_group, y = y, linewidth = legend_label), 
                group = 1)+
      # --- adjust axis values and labels ---
      scale_y_continuous(breaks  = pop_range_seq, labels = abs) +
      expand_limits(y = range(pop_range_seq)) +
      scale_x_continuous(sec.axis = dup_axis(), #numeric axis for duplication
                         breaks = 1:18, 
                         labels = unique(JoinedData$age_group),
                         expand = c(0, -0.05))+ #remove space between axis and plot
      # --- Flip coordinates and scale ---
      coord_flip(clip = "off") +
      scale_color_manual(values = 
                           c("#4292C6","#97BCDD","#EE6100","black"))+
        scale_fill_manual(values = c("#EE6100"))+
      # --- Manual color legend with line-only keys ---
      guides(linewidth = 
               guide_legend(override.aes = 
                              list(linetype = "solid",
                                   color = c("black","#4292C6", #different order (same as factor levels)
                                             "#97BCDD","#EE6100"), 
                                   alpha = c(1, alpha_bar, alpha_bar, 0.7),
                                   linewidth = 2,
                                   fill = NA),
                            title = ""))+
      theme_bw()+
      theme(axis.title = element_blank(),
            axis.text = element_text(size = 20, face = "bold"),
            legend.text = element_text(size = 20),
            panel.border = element_blank(),
            axis.line = element_line(colour ="black"),
            legend.background = element_blank(),
            legend.position = "inside",
            legend.position.inside = c(0.9, 0.87))+
      annotate("text", x = 19.0, y = tail(pop_range_seq, 2)[1],
               label = "Females", size = 8, hjust = 1, vjust = 1) +  # Left-aligned title
      annotate("text", x = 19.0, y = head(pop_range_seq, 2)[2],
               label = "Males", hjust = 1, vjust = 1, size = 8)
  
    suppressWarnings(print(Plot))
  
  
}


