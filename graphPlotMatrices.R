library(deSolve)
library(ggplot2)
library(reshape2)
library(rgdal)
library(sp)
setwd("C:/Users/Administrator/Desktop/SEIR")
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
  write.csv(beta, "betacheckgraph2.csv")
  # This is the callback function for ode
  deriv <- function(t, state, p){
    with(as.list(c(state, p)),{
      
      # convert all state variables to an areas(rows)
      # by stocks (cols) matrix
      states<-matrix(state,
                     nrow=num_areas,
                     ncol=num_stocks)
      
      # extract required state vectors for flow calculations 
      Susceptible <- states[,1]
      Exposed    	<- states[,2]
      Infected    <- states[,3]
      Recovered 	<- states[,4]
     
      # matrix operation to calculate lambda values
      Lambda      <- beta %*% Infected

      # calculate all flows
      IR         <- Lambda * Susceptible
      ER  	     <- Exposed * e
      RR         <- Infected / RDelay
      
      # specify net flows for each stock
      dS_dt      <- -IR
      dE_dt	     <- IR -ER
      dI_dt      <- ER - RR
      dR_dt      <- RR
      return (list(c(dS_dt,dE_dt,dI_dt,dR_dt)))  
    })
  }
  
  return (ode(y=InitStocks, 
              times=vtimes, 
              func = deriv, 
              parms=pars, 
              method="euler"))
}
d <- readOGR(dsn="census", "Census2011_Constituencies_2013")
# list of areas in the model (or cohorts)
areas                 <-d$NAME

# list of model stocks (SIR Model) S01_ etc needed for sorting operation
stocks                <-c("S01_Susceptible","S02_Exposed","S03_Infected","S04_Recovered")

# Sorted cartesian product of stocks x areas to get combinations of stock x area
allStocks             <- sort(apply(expand.grid(stocks,areas),1,function(x) paste(x,collapse="_")))
print("allStocks")
print(allStocks)
print(class(allStocks))
# Sample effective contact rates (area x area matrix)

C <- read.csv("Ce_non.csv")
CE  <- as.matrix(C)

# Initial sub-population values and names provided
InitPopulations        <- d$Total2011
names(InitPopulations) <- areas

I1 <- c(1)
I2 <- rep(0, length(areas) -1)
I <- c(I1, I2)
S <- InitPopulations-I
E <- rep(0, length(areas))
R <- rep(0, length(areas))


# Initial state values: orders (Sx4) (Ix4) (Rx4)
InitStocks            <- c(S,  # sequence of init susceptible
                           E,  # sequence of init exposed
                           I,  # sequence of init infected 
                           R)  # sequence of init recovered

names(InitStocks)     <- allStocks
#times=seq(0,50,by=0.01)
times=seq(0,20,by=0.01)
params=c(RDelay=2, e=.2)
out <- data.frame(solveSIRCohortModel(params,             # vector for the model parameters
                                      times,              # vector for solution times
                                      CE,                 # square matrix of effective contacts
                                      InitPopulations,     # vector of initial area populations
                                      InitStocks,         # vector of initial stock values 
                                      length(areas),      # number of areas
                                      length(stocks)))    # number of stocks


#print(names(out))
out2<-melt(out,id.vars="time")
#head(out2)
#View(out)
#View(out2)
print(names(d))
names(out2)[names(out2)=="value"] <- "population"
ggplot(out2,aes(time,population,color=variable))+geom_line()
#write.csv(out, "trainingData.csv")