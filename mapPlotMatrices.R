library(deSolve)
library(ggplot2)
library(reshape2)
library(Matrix)
library(rgdal)
library(plyr)
library(maptools)
library(RColorBrewer)
library(sp)

solveSIRCohortModel <- function(pars,         # parameter vector
                                vtimes,       # simulation time vector
                                CE,           # Effective contact matrix
                                InitPop,      # Initial population matrix (by area)
                                InitStocks,   # Initial stocks (area x stock)
                                num_areas,    # Total number of areas 
                                num_stocks)   # Total number of stocks per area
{
  # create the beta matrix for use in the simulation
  beta <- CE/InitPop  
  
  # This is the callback function for ode
  deriv <- function(t, state, p){
    with(as.list(c(state, p)),{ 
      plotmap <- function()
      {
        tryCatch({
          cols <- brewer.pal(4, "Reds")
          
          #where to save the image, and name of each image: 
          #t is used to keep track of the time step for each map.
          #name Infected Kerry/Images Kerry chosen, as that is the
          #county where the infection is starting in this example
          fname <- paste("C://Users//Administrator//Desktop//SEIR//ImagesKerry15Sep//Infected_Kerry",t, ".jpg", sep="")
        
          plot.new()
    
          #range of colours - green for no infection, red for high
          cols <- colorRampPalette(c("green", "red"))
          #creating the image - data to use, column required, how many colours, image title
          temp <- spplot(d, zcol="Infected", colorkey=TRUE, col.regions = cols(100), main=paste("Infected ","Kerry ",t,sep=""))
          
          png(filename=fname)
          #save the image
          print(temp)
          
          #as soon as plot hits a condition it has 0 constituencies to plot for
          #it crashes, instead, we want it to finish the map
        },
        error=function(cond) {
     
          return(NA)
        },
        warning=function(cond) {
       
          # Choose a return value in case of warning
          return(NULL)
        },
        finally={
      
        dev.off()
    
        }
        )    
      }
      # convert all state variables to an areas(rows) by stocks (cols) matrix
      states<-matrix(state,
                     nrow=num_areas,
                     ncol=num_stocks)
      
      # extract required state vectors for flow calculations 
      Susceptible <- states[,1]
      Exposed  <- states[,2]
      Infected    <- states[,3]
      Recovered 	<- states[,4]
      ProportionInfected <- Infected/InitPop
      # matrix operation to calculate lambda values
      Lambda      <- beta %*% Infected
      
      # calculate all flows
      IR         <- Lambda * Susceptible
      ER	     <- Exposed * e
      RR         <- Infected / RDelay
      
      # specify net flows for each stock
      dS_dt      <- -IR
      dE_dt	     <- IR -ER
      dI_dt      <- ER - RR
      dR_dt      <- RR
      d <- spCbind(d,Susceptible)
      d <- spCbind(d,Exposed)
      d <- spCbind(d,Infected)
      d <- spCbind(d,Recovered)
      d <- spCbind(d, ProportionInfected)
      
      
      plotmap()
    #details to return for the model, and extra details wanted 
      return (list(c(dS_dt,dE_dt,dI_dt,dR_dt),ProportionInfected))  
    })
  }
  
  return (ode(y=InitStocks, 
              times=vtimes, 
              func = deriv, 
              parms=pars, 
              method="euler"))
}
setwd("C:/Users/Administrator/Desktop/SEIR/")
d <- readOGR(dsn="census", "Census2011_Constituencies_2013")
# list of areas in the model (or cohorts)
areas                 <-d$NAME

# list of model stocks (SIR Model) S01_ etc needed for sorting operation
stocks                <-c("S01_Susceptible","S02_Exposed","S03_Infected","S04_Recovered")

# Sorted cartesian product of stocks x areas to get combinations of stock x area
allStocks             <- sort(apply(expand.grid(stocks,areas),1,function(x) paste(x,collapse="_")))

# Sample effective contact rates (area x area matrix)
C <- read.csv("Ce_non.csv")
CE                    <- as.matrix(C)

# Initial sub-population values and names provided
InitPopulations        <- d$Total2011
names(InitPopulations) <- areas


#where we're starting the infection. 
#our data has 40 consituencies, we want the
#infection to start in the 5th in the list (Kerry)
G1 <- rep(0, 4)
G2 <- c(1)
G3 <- rep(0, 35)
#Infection vector
I <- c(G1, G2, G3)
#Susceptible vector
S <- InitPopulations-I
#Exposed vector
E <- rep(0, length(areas))
#Recovered vector
R <- rep(0, length(areas))


# Initial state values: orders (Sx4) (Ix4) (Rx4)
InitStocks            <- c(S,  # sequence of init susceptible
                           E,  # sequence of init exposed
                           I,  # sequence of init infected 
                           R)  # sequence of init recovered
#name the compartments
names(InitStocks)     <- allStocks
#time sequence for the model
times=seq(0,50,by=0.01)
#parameters (Recovery delay, and exposure delay)
params=c(RDelay=5, e=.2)
#call the model
out <- data.frame(solveSIRCohortModel(params,             # vector for the model parameters
                                      times,              # vector for solution times
                                      CE,                 # square matrix of effective contacts
                                      InitPopulations,     # vector of initial area populations
                                      InitStocks,         # vector of initial stock values 
                                      length(areas),      # number of areas
                                      length(stocks)))    # number of stocks


#reshape into tidy data for plotting graphs
out2<-melt(out,id.vars="time")
#plot the graphs (same as graphPlotMatrices...)
ggplot(out2,aes(time,value,color=variable))+geom_line()