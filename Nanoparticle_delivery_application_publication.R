

# *** metrics ***

# The metric used for the optimization
SODI <- function(observed, predicted, comp.names =NULL){
  # Check if the user provided the correct input format
  if (!is.list(observed) || !is.list(predicted)){
    stop(" The observations and predictions must be lists")
  }
  # Check if the user provided equal length lists
  if (length(observed) != length(predicted)){
    stop(" The observations and predictions must have the same compartments")
  }
  Ncomp <- length(observed) # Number of compartments
  I <- rep(NA, Ncomp) # Compartment discrepancy index
  N_obs <- rep(NA, Ncomp) #Number of observations per compartment
  #loop over the compartments
  for (i in 1:Ncomp){
    Et <- 0 #relative error with observations
    St <- 0  #relative error with simulations
    N <- length(observed[[i]]) # number of observations for compartment i
    # Check if observations and predictions have equal length
    if(N != length(predicted[[i]])){
      stop(paste0("Compartment ",i," had different length in the observations and predictions"))
    }
    N_obs[i] <- N # populate the N_obs vector
    for (j in 1:N){
      # sum of relative squared errors (error = observed - predicted)
      Et <- Et + ( abs(observed[[i]][j] - predicted[[i]][j])  / observed[[i]][j] )  ^2
      St <- St + ( abs(observed[[i]][j] - predicted[[i]][j])  / predicted[[i]][j] )  ^2
    }
    
    # root mean of the square of observed values
    RMEt <- sqrt(Et/N)
    # root mean of the square of simulated values
    RMSt <- sqrt( St/N)
    
    I[i] <- (RMEt + RMSt)/2
  }
  # Total number of observations
  Ntot <- sum(N_obs)
  # Initialise the consolidated discrepancy index
  Ic <-0
  for (i in 1:Ncomp){
    # Give weight to compartments with more observations (more information)
    Ic <- Ic +  I[i]* N_obs[i]/Ntot
  }
  # Name the list of compartment discrepancy indices
  if ( !is.null(comp.names)){
    names(I) <- comp.names
  }else if (!is.null(names(observed))){
    names(I) <- names(observed)
  } else if (!is.null(names(predicted)) && is.null(comp.names) ){
    names(I) <- names(predicted)
  }
  return(Ic)
  #return(list(Total_index = Ic, Compartment_index= I))
}



#===============================================
#2. Function to create initial values for ODEs 
#===============================================

create.inits <- function( dose){
  NP_AS<-dose; NP_CV<-0;NP_CI <- 0; NP_OT <- 0; NP_EX = 0
  
  return(c( "NP_AS" = NP_AS, "NP_CV"=NP_CV,"NP_CI" = NP_CI,
            "NP_OT"=NP_OT,  "NP_EX" = NP_EX))
}
#==============
#3. ODEs System
#==============
Wu_extended_model <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    # Description:
    # k_AStCV: Administrtion Site to target Cell Vicinity
    # k_CVtAS: Cell Vicinity to  Administration Site
    # k_CVtCI: Cell Vicinity to  Cell Interior
    # k_AStOT: Administration Site to Off-Target Sites
    # k_OTtAS: Off-Target sites to Administration Site
    # k_CItCV: Cell Interior to Cell Vicinity
    # k_CVtOT: Cell Vicinity to Off-Target sites
    # k_OTtCV: Off-Target sites to Cell Vicinity

    # NP_AS: nanoparticles in Administration Site
    # NP_CV: nanoparticles in target Cell Vicinity
    # NP_CI: nanoparticles in target Cell Interior
    # NP_OT: nanoparticles in Off-Target sites
    # NP_EX: nanoparticles in excreta
    
    # Units:
    # k_TTC, k_C, k_ITC, k_ATC ---> 1/h
    # NP_AS, NP_CV, NP_CI, NP_OT ---> 1/

    k_CVtCI <- 0.400#Dubaj et al. 2022
    k_CItCV <- 0.0598  #Dubaj et al. 2022

    
    dNP_AS <-  -k_AStCV*NP_AS + k_CVtAS*NP_CV - k_AStOT*NP_AS + k_OTtAS*NP_OT-k_AStEX* NP_AS 
    dNP_CV <-  k_AStCV*NP_AS - k_CVtAS*NP_CV- k_CVtCI*NP_CV + k_CItCV*NP_CI -
               k_CVtOT*NP_CV + k_OTtCV*NP_OT
    dNP_CI <-  k_CVtCI*NP_CV - k_CItCV*NP_CI
    dNP_OT <-  k_AStOT*NP_AS -  k_OTtAS*NP_OT+ k_CVtOT*NP_CV -  k_OTtCV*NP_OT- k_OTtEX * NP_OT
    dNP_EX <- k_AStEX* NP_AS +  k_OTtEX * NP_OT
    
    NP_tot <- NP_AS + NP_CV + NP_OT + NP_CI+ NP_EX
    
    efficiency <- NP_CI/NP_tot
    
    return(list(c("dNP_AS" = dNP_AS,   "dNP_CV" = dNP_CV,
                  "dNP_CI" = dNP_CI, "dNP_OT" = dNP_OT, "dNP_EX" = dNP_EX), 
                "efficiency" = efficiency))
  })
}



obj_func <- function(x, dose, df, metric = "AAFE"){
  
  BodyBurden <- c( df$admin_site, df$lungs, df$off_target, df$excreta)

  
  parms <- c("k_AStCV" = exp(x[1]),  "k_CVtAS" = exp(x[2]),
             "k_AStOT" =  exp(x[3]) , 
             "k_OTtAS" = exp(x[4]),  
             "k_CVtOT" =  exp(x[5]),
             "k_OTtCV"=  exp(x[6]), "k_AStEX" =  exp(x[7]),  "k_OTtEX" =  exp(x[8]))
  
  sol_times <- c(seq(0,1, 0.001),seq(1.1,5, 0.1),  seq(6,28*24, 1))
  inits <- create.inits(unname(dose))
  solution <- data.frame(deSolve::ode(times = sol_times,  func = Wu_extended_model,
                                      y = inits, parms = parms, method="lsodes",
                                      rtol = 1e-7, atol = 1e-7))
  if(sum(solution$time %in% df$time) == length( df$time)){
    results <- c(solution[solution$time %in% df$time, "NP_AS"],
                 (solution[solution$time %in% df$time, "NP_CV"]+ solution[solution$time %in% df$time, "NP_CI"]),
                 solution[solution$time %in% df$time, "NP_OT"],
                 solution[solution$time %in% df$time, "NP_EX"])
  }else{
    results <-   BodyBurden *100
    
  }
  
  
  # Find the position of the current PFAS in the PFAS_names vector
  if(metric == "AAFE"){
    score <- AAFE(BodyBurden, results) 
  }else if (metric =="rmse"){
    score <- rmse(BodyBurden, results)
  }else if(metric == "SODI"){
    score <- SODI(list(BodyBurden), list(results))
  }       
  return(score)
}

##############################################################################
################################################################################
# Load the data for PFAS concentration
data_concentration <- t(openxlsx::read.xlsx ('PEG-AU-NPs.xlsx', 
                                             sheet = "Kozics_concentration", 
                                             colNames = TRUE, rowNames = TRUE))
data_mass <- t(openxlsx::read.xlsx ('PEG-AU-NPs.xlsx', 
                                    sheet = "Kozics_mass", 
                                    colNames = TRUE, rowNames = TRUE)[1:5,])


dose_per_weight <- 0.7 #mg/kg
dose <- dose_per_weight*229
#The percentage of gold recovery in the analyzed samples
#detected 1 h after i.v. injection accounted for approximately 68.5% of the total injected
#dose.
#dose <- rowSums(data_mass)[1]/0.685 
df <- data.frame(time = c(1, 4, 24, 7*24, 28*24), lungs = unname(data_mass[,"lungs"]), 
                 off_target = unname(data_mass[,"liver"]+ data_mass[,"spleen"] + 
                                             data_mass[,"kidneys"]),  
                 admin_site = unname(data_mass[,"blood"]),
                 excreta = rep(NA, 5))

df$excreta <- dose - (df$lungs+df$off_target+df$admin_site)

opts <- list( "algorithm" = "NLOPT_LN_SBPLX", #"NLOPT_LN_NEWUOA",NLOPT_LN_SBPLX , #"NLOPT_LN_BOBYQA" #"NLOPT_LN_COBYLA"
              "xtol_rel" = 1e-7, 
              "ftol_rel" = 1e-7,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0 ,
              "maxeval" = 5000,
              "print_level" = 1)

N_pars <- 8
# Define initial values of fitted parameters to provide to the optimization routine
x0 <-  rep(1,N_pars)
set.seed(1221312)
optimization<- nloptr::nloptr(x0 = x0,
                              eval_f = obj_func,
                              opts = opts,
                              dose = dose,
                              df = df,
                              metric = "SODI")

parms <- c("k_AStCV" =exp(optimization$solution[1]), 
           "k_CVtAS" = exp(optimization$solution[2]),
           "k_AStOT" = exp(optimization$solution[3]) , 
           "k_OTtAS" = exp(optimization$solution[4]), "k_CVtOT" =  exp(optimization$solution[5]),
           "k_OTtCV"=  exp(optimization$solution[6]), "k_AStEX" =  exp(optimization$solution[7]), 
           "k_OTtEX" =  exp(optimization$solution[8]))


sol_times <- seq(0,28*24, 1 )
inits <- create.inits(unname(dose))
solution <- data.frame(deSolve::ode(times = sol_times,  func = Wu_extended_model,
                                    y = inits, parms = parms, method="lsodes",
                                    rtol = 1e-7, atol = 1e-7))


library(ggplot2)
# Lungs plot
p1 <- ggplot()+
  geom_line(data = solution, aes(x=time, y= (NP_CV+NP_CI), color = "Model Predictions" ), size=1.7)+
  geom_point(data = df, aes(x=time  , y=lungs, color = "Experimental Observations" ), size=5)+
  labs(title = "Cell Vicinity & Cell Interior", y = "PEG-AU NPs mass (ug)", x = "Time (hours)")+
  scale_color_manual("",values = c( "Model Predictions"  = "#000000",
                                    "Experimental Observations" = "#000000"), 
                     guide = guide_legend(override.aes =
                                            list(shape = c(16, NA),
                                                 linetype = c(0,1))))+
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5,size=20), 
        axis.title.y =element_text(hjust = 0.5,size=18),
        axis.text.y=element_text(size=16),
        axis.title.x =element_text(hjust = 0.5,size=18),
        axis.text.x=element_text(size=16),
        legend.title=element_text(hjust = 0.5,size=18), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1.5, 'cm'))

########
# Off-target sites plot
p2 <- ggplot()+
  geom_line(data = solution, aes(x=time, y=NP_OT, color = "Model Predictions" ), size=1.7)+
  geom_point(data = df, aes(x=time  , y=off_target, color = "Experimental Observations" ), size=5)+
  
  labs(title = "Off target sites", y = "PEG-AU NPs mass (ug)", x = "Time (hours)")+
  scale_color_manual("",values = c( "Model Predictions"  = "#000000",
                                    "Experimental Observations" = "#000000"), 
                     guide = guide_legend(override.aes =
                                            list(shape = c(16, NA),
                                                 linetype = c(0,1))))+
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5,size=20), 
        axis.title.y =element_text(hjust = 0.5,size=18),
        axis.text.y=element_text(size=16),
        axis.title.x =element_text(hjust = 0.5,size=18),
        axis.text.x=element_text(size=16),
        legend.title=element_text(hjust = 0.5,size=18), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1.5, 'cm'))



########
# Excreta plot
p3 <- ggplot()+
  geom_line(data = solution, aes(x=time, y=NP_EX, color = "Model Predictions" ), size=1.7)+
  geom_point(data = df, aes(x=time  , y=excreta, color = "Experimental Observations" ), size=5)+
  labs(title = "Excreta", y = "PEG-AU NPs mass (ug)", x = "Time (hours)")+
  scale_color_manual("",values = c( "Model Predictions"  = "#000000",
                                    "Experimental Observations" = "#000000"), 
                     guide = guide_legend(override.aes =
                                            list(shape = c(16, NA),
                                                 linetype = c(0,1))))+
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5,size=20), 
        axis.title.y =element_text(hjust = 0.5,size=18),
        axis.text.y=element_text(size=16),
        axis.title.x =element_text(hjust = 0.5,size=18),
        axis.text.x=element_text(size=16),
        legend.title=element_text(hjust = 0.5,size=18), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1.5, 'cm'))

########
# Administration Site plot
p4 <- ggplot()+
  geom_line(data = solution, aes(x=time, y=NP_AS, color = "Model Predictions" ), size=1.7)+
  geom_point(data = df, aes(x=time  , y=admin_site, color = "Experimental Observations" ), size=5)+
  labs(title = "Administration Site", y = "PEG-AU NPs mass (ug)", x = "Time (hours)")+
  scale_color_manual("",values = c( "Model Predictions"  = "#000000",
                                    "Experimental Observations" = "#000000"), 
                     guide = guide_legend(override.aes =
                                            list(shape = c(16, NA),
                                                 linetype = c(0,1))))+
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5,size=20), 
        axis.title.y =element_text(hjust = 0.5,size=18),
        axis.text.y=element_text(size=16),
        axis.title.x =element_text(hjust = 0.5,size=18),
        axis.text.x=element_text(size=16),
        legend.title=element_text(hjust = 0.5,size=18), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1.5, 'cm'))


final_plot<-ggpubr::ggarrange(p4,p3, p2, p1, ncol=2, nrow=2, 
                              common.legend = TRUE, legend="bottom")
plot.margin=unit(c(0,0,0,0), "pt")


dev.off()


