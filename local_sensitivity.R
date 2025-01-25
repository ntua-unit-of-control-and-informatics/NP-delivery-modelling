# Working directory

setwd("C:/Users/ptsir/Desktop/projects/Papers/Nanoparticle_delivery_application")


# *** metrics ***

# The metric used for the optimization
mse_custom <- function(observed, predicted){
  mean((observed - predicted)^2)
}

mape <- function(observed, predicted){
  mean(abs(observed-predicted)*100/observed)
}

rmse <- function(observed, predicted){
  sqrt(mean((observed-predicted)^2)) 
}

AAFE <- function(observations, predictions, times=NULL){
  y_obs <- unlist(observations)
  y_pred <- unlist(predictions)
  # Total number of observations
  N<- length(y_obs)
  log_ratio <- rep(NA, N) 
  for ( i in 1:N){
    log_ratio[i] <- abs(log((y_pred[i]/y_obs[i]), base = 10))
  }
  aafe <- 10^(sum(log_ratio)/N) 
  return(aafe)
}

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
Rat_model <- function(time, inits, params){
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
    # k_AStEX: Administration site to Excreta
    # k_OTtEX: Off-Target sites to Excreta
    
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
    
    efficiency <- 100*NP_CI/NP_tot
    
    return(list(c("dNP_AS" = dNP_AS,   "dNP_CV" = dNP_CV,
                  "dNP_CI" = dNP_CI, "dNP_OT" = dNP_OT, "dNP_EX" = dNP_EX), 
                "efficiency" = efficiency))
  })
}


# Function for estimating perturbations of response based on changes in size
perturbations <- function(param_index, dose =  160.3){
  parms <- c("k_AStCV" = 0.01783763, "k_CVtAS" = 0.5687722, "k_AStOT" =  0.01587592 , 
             "k_OTtAS" =0, "k_CVtOT" =  1.532472, "k_OTtCV"=   0.008731313,
             "k_AStEX" =   0.07578454, 
             "k_OTtEX" =  0)
  #Store perturbed params in an appropriate vector 
  perturbed_params <- parms
  inits <<- create.inits(unname(dose))
  changes <- c(1, 1.5, 1.99 , 0.5, 0.01) #fold change\
  sol_times <- seq(0,28*24, 1 )
  #initialise result to store ODE solutions
  result <- data.frame(time = sol_times, baseline = rep(NA, length(sol_times)),
                       plus50 = rep(NA, length(sol_times)),  plus99= rep(NA, length(sol_times)),
                       minus50  = rep(NA, length(sol_times)), minus99 = rep(NA, length(sol_times)))
  for (i in 1:5){
  #change the selected parameter
  perturbed_params[param_index] <- parms[param_index]*changes[i]
  
  result[i+1] <- data.frame(deSolve::ode(times = sol_times,  func = Rat_model,
                                      y = inits, parms = perturbed_params, method="lsodes",
                                      rtol = 1e-7, atol = 1e-7))$efficiency
  }
  return(result)
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

df_k_AStCV <- perturbations(param_index = 1)
df_k_CVtAS <- perturbations(param_index = 2)
df_k_AStOT <- perturbations(param_index = 3)
df_k_OTtAS <- perturbations(param_index = 4)
df_k_CVtOT <- perturbations(param_index = 5)
df_k_OTtCV <- perturbations(param_index = 6)
df_k_AStEX <- perturbations(param_index = 7)
df_k_OTtEX <- perturbations(param_index = 8)

cls <- c("+50%" = "#E69F00", "-50%" = "#56B4E9",
         "No perturbation" ="#000000", "-99%" = "#009E73", "+99%" ="#CC79A7")
library(ggplot2)
# Lungs plot
p1 <- ggplot()+
  geom_line(data = df_k_AStCV, aes(x=time, y= baseline, color = "No perturbation" ), size=1.7)+
  geom_line(data = df_k_AStCV, aes(x=time, y= plus50, color = "+50%" ), size=1.7)+
  geom_line(data = df_k_AStCV, aes(x=time, y= minus50, color = "-50%" ), size=1.7)+
  geom_line(data = df_k_AStCV, aes(x=time, y= plus99, color = "+99%" ), size=1.7)+
  geom_line(data = df_k_AStCV, aes(x=time, y= minus99, color = "-99%" ), size=1.7)+
  labs(title =  expression("k"[AStCV]), y = "Efficiency (%)", x = "Time (hours)")+
  ylim(0, 4.5)+
  scale_color_manual("",values = cls)+
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5,size=20), 
        axis.title.y =element_text(hjust = 0.5,size=18),
        axis.text.y=element_text(size=16),
        axis.title.x =element_text(hjust = 0.5,size=18),
        axis.text.x=element_text(size=16),
        legend.title=element_text(hjust = 0.5,size=18), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1.5, 'cm'))

p2 <- ggplot()+
  geom_line(data = df_k_CVtAS, aes(x=time, y= baseline, color = "No perturbation" ), size=1.7)+
  geom_line(data = df_k_CVtAS, aes(x=time, y= plus50, color = "+50%" ), size=1.7)+
  geom_line(data = df_k_CVtAS, aes(x=time, y= minus50, color = "-50%" ), size=1.7)+
  geom_line(data = df_k_CVtAS, aes(x=time, y= plus99, color = "+99%" ), size=1.7)+
  geom_line(data = df_k_CVtAS, aes(x=time, y= minus99, color = "-99%" ), size=1.7)+
  ylim(0,4)+
  labs(title = expression("k"[CVtAS]), y = "Efficiency (%)", x = "Time (hours)")+
  scale_color_manual("",values = cls)+
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5,size=20), 
        axis.title.y =element_text(hjust = 0.5,size=18),
        axis.text.y=element_text(size=16),
        axis.title.x =element_text(hjust = 0.5,size=18),
        axis.text.x=element_text(size=16),
        legend.title=element_text(hjust = 0.5,size=18), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1.5, 'cm'))

p3 <- ggplot()+
  geom_line(data = df_k_AStOT, aes(x=time, y= baseline, color = "No perturbation" ), size=1.7)+
  geom_line(data = df_k_AStOT, aes(x=time, y= plus50, color = "+50%" ), size=1.7)+
  geom_line(data = df_k_AStOT, aes(x=time, y= minus50, color = "-50%" ), size=1.7)+
  geom_line(data = df_k_AStOT, aes(x=time, y= plus99, color = "+99%" ), size=1.7)+
  geom_line(data = df_k_AStOT, aes(x=time, y= minus99, color = "-99%" ), size=1.7)+
  ylim(0,4)+
  labs(title = expression("k"[AStOT]), y = "Efficiency (%)", x = "Time (hours)")+
  scale_color_manual("",values = cls)+
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5,size=20), 
        axis.title.y =element_text(hjust = 0.5,size=18),
        axis.text.y=element_text(size=16),
        axis.title.x =element_text(hjust = 0.5,size=18),
        axis.text.x=element_text(size=16),
        legend.title=element_text(hjust = 0.5,size=18), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1.5, 'cm'))


p4 <- ggplot()+
  geom_line(data = df_k_OTtAS, aes(x=time, y= baseline, color = "No perturbation" ), size=1.7)+
  geom_line(data = df_k_OTtAS, aes(x=time, y= plus50, color = "+50%" ), size=1.7)+
  geom_line(data = df_k_OTtAS, aes(x=time, y= minus50, color = "-50%" ), size=1.7)+
  geom_line(data = df_k_OTtAS, aes(x=time, y= plus99, color = "+99%" ), size=1.7)+
  geom_line(data = df_k_OTtAS, aes(x=time, y= minus99, color = "-99%" ), size=1.7)+
  ylim(0, 4.5)+
  labs(title = expression("k"[OTtAS]), y = "Efficiency (%)", x = "Time (hours)")+
  scale_color_manual("",values = cls)+
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5,size=20), 
        axis.title.y =element_text(hjust = 0.5,size=18),
        axis.text.y=element_text(size=16),
        axis.title.x =element_text(hjust = 0.5,size=18),
        axis.text.x=element_text(size=16),
        legend.title=element_text(hjust = 0.5,size=18), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1.5, 'cm'))

p5 <- ggplot()+
  geom_line(data = df_k_CVtOT, aes(x=time, y= baseline, color = "No perturbation" ), size=1.7)+
  geom_line(data = df_k_CVtOT, aes(x=time, y= plus50, color = "+50%" ), size=1.7)+
  geom_line(data = df_k_CVtOT, aes(x=time, y= minus50, color = "-50%" ), size=1.7)+
  geom_line(data = df_k_CVtOT, aes(x=time, y= plus99, color = "+99%" ), size=1.7)+
  geom_line(data = df_k_CVtOT, aes(x=time, y= minus99, color = "-99%" ), size=1.7)+
  ylim(0, 4.5)+
  labs(title = expression("k"[CVtOT]), y = "Efficiency (%)", x = "Time (hours)")+
  scale_color_manual("",values = cls)+
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5,size=20), 
        axis.title.y =element_text(hjust = 0.5,size=18),
        axis.text.y=element_text(size=16),
        axis.title.x =element_text(hjust = 0.5,size=18),
        axis.text.x=element_text(size=16),
        legend.title=element_text(hjust = 0.5,size=18), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1.5, 'cm'))

p6 <- ggplot()+
  geom_line(data = df_k_OTtCV, aes(x=time, y= baseline, color = "No perturbation" ), size=1.7)+
  geom_line(data = df_k_OTtCV, aes(x=time, y= plus50, color = "+50%" ), size=1.7)+
  geom_line(data = df_k_OTtCV, aes(x=time, y= minus50, color = "-50%" ), size=1.7)+
  geom_line(data = df_k_OTtCV, aes(x=time, y= plus99, color = "+99%" ), size=1.7)+
  geom_line(data = df_k_OTtCV, aes(x=time, y= minus99, color = "-99%" ), size=1.7)+
  ylim(0, 4.5)+
  labs(title = expression("k"[OTtCV]), y = "Efficiency (%)", x = "Time (hours)")+
  scale_color_manual("",values = cls)+
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5,size=20), 
        axis.title.y =element_text(hjust = 0.5,size=18),
        axis.text.y=element_text(size=16),
        axis.title.x =element_text(hjust = 0.5,size=18),
        axis.text.x=element_text(size=16),
        legend.title=element_text(hjust = 0.5,size=18), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1.5, 'cm'))

p7 <- ggplot()+
  geom_line(data = df_k_AStEX, aes(x=time, y= baseline, color = "No perturbation" ), size=1.7)+
  geom_line(data = df_k_AStEX, aes(x=time, y= plus50, color = "+50%" ), size=1.7)+
  geom_line(data = df_k_AStEX, aes(x=time, y= minus50, color = "-50%" ), size=1.7)+
  geom_line(data = df_k_AStEX, aes(x=time, y= plus99, color = "+99%" ), size=1.7)+
  geom_line(data = df_k_AStEX, aes(x=time, y= minus99, color = "-99%" ), size=1.7)+
  ylim(0, 4.5)+
  labs(title = expression("k"[AStEX]), y = "Efficiency (%)", x = "Time (hours)")+
  scale_color_manual("",values = cls)+
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5,size=20), 
        axis.title.y =element_text(hjust = 0.5,size=18),
        axis.text.y=element_text(size=16),
        axis.title.x =element_text(hjust = 0.5,size=18),
        axis.text.x=element_text(size=16),
        legend.title=element_text(hjust = 0.5,size=18), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1.5, 'cm'))

p8 <- ggplot()+
  geom_line(data = df_k_OTtEX, aes(x=time, y= baseline, color = "No perturbation"  ), size=1.7)+
  geom_line(data = df_k_OTtEX, aes(x=time, y= plus50, color = "+50%" ), size=1.7)+
  geom_line(data = df_k_OTtEX, aes(x=time, y= minus50, color = "-50%" ), size=1.7)+
  geom_line(data = df_k_OTtEX, aes(x=time, y= plus99, color = "+99%" ), size=1.7)+
  geom_line(data = df_k_OTtEX, aes(x=time, y= minus99, color = "-99%" ), size=1.7)+
  ylim(0, 4.5)+
  labs(title = expression("k"[OTtEX]), y = "Efficiency (%)", x = "Time (hours)")+
  scale_color_manual("",values = cls)+
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5,size=20), 
        axis.title.y =element_text(hjust = 0.5,size=18),
        axis.text.y=element_text(size=16),
        axis.title.x =element_text(hjust = 0.5,size=18),
        axis.text.x=element_text(size=16),
        legend.title=element_text(hjust = 0.5,size=18), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1.5, 'cm'))

final_plot<-ggpubr::ggarrange(p1,p2, p3, p5, p6, p7, ncol=2, nrow=3, 
                              common.legend = TRUE, legend="bottom")
plot.margin=unit(c(0,0,0,0), "pt")

# Save the plot with dynamically adjusted dimensions
ggsave("efficiency.png", plot = final_plot,
       device = 'png', dpi = 300,
       width = 13,
       height = 10,
       units = "in")
dev.off()





