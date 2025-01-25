# Set your Working directory
setwd("C:/Insert_your_working_directory_here")

# Function to install packages if not already present
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

# List of packages to install
packages_to_install <- c("ggplot2", "deSolve")

# Call the installation function
install_if_missing(packages_to_install)

# Load the data 
df <- read.csv('example.csv')
dose <- 160 #ug

# Parameter optimisation is realised in log space
parameters <- list(
  k_AStCV = log(1.0),
  k_CVtAS = log(1.0),
  k_AStOT = log(1.0),
  k_OTtAS = log(1.0),
  k_CVtOT = log(1.0),
  k_OTtCV = log(1.0),
  k_AStEX = log(1.0),
  k_OTtEX = log(1.0),
  k_CVtCI = log(1.0),
  k_CItCV = log(1.0),
  k_CItOT = log(1.0),
  k_DEG = log(1.0)
)

# Alter the list, enlarge or reduce, based on the NP and available data
fixed_params <- list(
  "k_DEG" = log(0), # Comment out to apply NP biodegradation
  "k_CItOT"  = log(0), #Comment out if target cell is liver
  "k_CVtCI" = log(0.400), 
  "k_CItCV" = log(0.0598)
  
)

# Update parameters with fixed values
parameters[names(fixed_params)] <- fixed_params

# Identify parameters to estimate
parameters_to_estimate <- setdiff(names(parameters), names(fixed_params))
N_pars <- length(parameters_to_estimate)

# Extract only the parameters to estimate
init_params <- parameters[parameters_to_estimate]

  
opts <- list( "algorithm" = "NLOPT_LN_SBPLX", #"NLOPT_LN_NEWUOA",NLOPT_LN_SBPLX , #"NLOPT_LN_BOBYQA" #"NLOPT_LN_COBYLA"
              "xtol_rel" = 1e-7, 
              "ftol_rel" = 1e-7,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0 ,
              "maxeval" = 5000,
              "print_level" = 1)


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
  NP_AS<-dose; NP_CV<-0;NP_CI <- 0; NP_OT <- 0; NP_EX = 0;NP_DEG = 0 
  
  return(c( "NP_AS" = NP_AS, "NP_CV"=NP_CV,"NP_CI" = NP_CI,
            "NP_OT"=NP_OT,  "NP_EX" = NP_EX, "NP_DEG" = NP_DEG))
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
    # k_DEG: degradation rate

    # NP_AS: nanoparticles in Administration Site
    # NP_CV: nanoparticles in target Cell Vicinity
    # NP_CI: nanoparticles in target Cell Interior
    # NP_OT: nanoparticles in Off-Target sites
    # NP_EX: nanoparticles in excreta
    
    # Units:
    # k_TTC, k_C, k_ITC, k_ATC ---> 1/h
    # NP_AS, NP_CV, NP_CI, NP_OT ---> ug


    
    dNP_AS <-  -k_AStCV*NP_AS + k_CVtAS*NP_CV - k_AStOT*NP_AS + k_OTtAS*NP_OT - 
                k_AStEX* NP_AS - k_DEG*NP_AS 
    dNP_CV <-  k_AStCV*NP_AS - k_CVtAS*NP_CV- k_CVtCI*NP_CV + k_CItCV*NP_CI -
               k_CVtOT*NP_CV + k_OTtCV*NP_OT - k_DEG*NP_CV
    dNP_CI <-  k_CVtCI*NP_CV - k_CItCV*NP_CI - k_CItOT*NP_CI -  k_DEG*NP_CI
    dNP_OT <-  k_AStOT*NP_AS -  k_OTtAS*NP_OT+ k_CVtOT*NP_CV -  k_OTtCV*NP_OT- 
               k_OTtEX * NP_OT + k_CItOT*NP_CI - k_DEG*NP_OT
    dNP_EX <- k_AStEX* NP_AS +  k_OTtEX * NP_OT
    dNP_DEG <- k_DEG*NP_OT + k_DEG*NP_CI + k_DEG*NP_CV+ k_DEG*NP_AS 
    
    NP_tot <- NP_AS + NP_CV + NP_OT + NP_CI+ NP_EX
    
    efficiency <- NP_CI/NP_tot
    
    return(list(c("dNP_AS" = dNP_AS,   "dNP_CV" = dNP_CV,
                  "dNP_CI" = dNP_CI, "dNP_OT" = dNP_OT, "dNP_EX" = dNP_EX, "dNP_DEG" = dNP_DEG), 
                "efficiency" = efficiency))
  })
}



obj_func <- function(x, dose, df, fixed_params,parameters_to_estimate, metric = "AAFE"){
  
  BodyBurden <- c( df$admin_site, df$target, df$off_target, df$excreta)

  
  estimated_params <- exp(unlist(x))
  names(estimated_params) <- parameters_to_estimate
  
  parms <- c(exp(unlist(fixed_params)), estimated_params)
  sol_times <- c(seq(0,1, 0.001),seq(1.1,5, 0.1),  seq(6,df$time[length(df$time)], 1))
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

# Define initial values of fitted parameters to provide to the optimization routine
x0 <-  unlist(parameters[parameters_to_estimate])
set.seed(1221312)
optimization<- nloptr::nloptr(x0 = x0,
                              eval_f = obj_func,
                              opts = opts,
                              dose = dose,
                              df = df,
                              fixed_params = fixed_params,
                              parameters_to_estimate = parameters_to_estimate,
                              metric = "SODI")

parameters[parameters_to_estimate] <- optimization$solution
parms<- exp(unlist(parameters))

write.csv(data.frame("parameters" = unlist(parms)), "parameters.csv")

sol_times <- seq(0,28*24, 1 )
inits <- create.inits(unname(dose))
solution <- data.frame(deSolve::ode(times = sol_times,  func = Wu_extended_model,
                                    y = inits, parms = parms, method="lsodes",
                                    rtol = 1e-7, atol = 1e-7))


library(ggplot2)
# Lungs plot
p1 <- ggplot()+
  geom_line(data = solution, aes(x=time, y= (NP_CV+NP_CI), color = "Model Predictions" ), size=1.7)+
  geom_point(data = df, aes(x=time  , y=target, color = "Experimental Observations" ), size=5)+
  labs(title = "Cell Vicinity & Cell Interior", y = "NP mass (ug)", x = "Time (hours)")+
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
  
  labs(title = "Off target sites", y = "NP mass (ug)", x = "Time (hours)")+
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
  labs(title = "Excreta", y = "NP mass (ug)", x = "Time (hours)")+
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
  labs(title = "Administration Site", y = "NP mass (ug)", x = "Time (hours)")+
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
print(final_plot)



