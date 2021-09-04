filepath="~/AmNatNote/Outputs"
setwd(filepath)


# Functions 

'%!in%' <- function(x,y)!('%in%'(x,y))
mat.torus <- function(Matrix,Rad,xcord,ycord){      # Torus of Space
  
  dm    <- nrow(Matrix)                               # dimension of matrix (n x n marix, so rows = col)
  
  Crown.Pot       <-  c(seq(1,Rad,by=2)^2)-1    # arbitrarily set to 10, because that many neighbors gets silly biologically, just used for computation 
  Crowns          <-  seq(0,length(Crown.Pot))
  Crown           <-  Crowns[min(which(Crown.Pot >= (Rad-1)^2))]  #returns crown exension total (even if not whole crown)
  
  rows  <- c(  (xcord-Crown):(xcord+Crown)       )    # figure out which rows of matrix crown falls in
  cols  <- c(  (ycord-Crown):(ycord+Crown)       )    # figure out which columns of matrix crown falls in
  
  rows[which(rows<1)]    <- rows[which(rows<1)] + dm  # if crown extends to a row value less than 1, go to opposite side of torus 
  rows[which(rows>dm)]   <- rows[which(rows>dm)] - dm # if crown extends to a row value greater than dm, go to opposite side of torus
  
  cols[which(cols<1)]    <- cols[which(cols<1)] + dm  # if crown extends to a column value less than 1, go to opposite side of torus 
  cols[which(cols>dm)]   <- cols[which(cols>dm)] - dm # if crown extends to a column value greater than dm, go to opposite side of torus
  
  JC_Matrix              <- Matrix[rows,cols ]        # returns subset of matrix / trees that are in JCE zone + extras 
  
  return(JC_Matrix)
  
}


############## create Spatial Grid of Environment Types ####################


#Basic Species info
dm <- 200
S <- 150
S.list <- seq(1:S)

set.seed(181)
S.list <- seq(1:S)
TimeSteps <- 25000
dd <- sample(1:S, dm*dm, replace = TRUE)
Mat.S <- matrix(dd,nrow=dm,ncol=dm)
df.Props                <- data.frame( matrix(NA,ncol=S+1,nrow=(1+TimeSteps) ))
df.Props[,1]            <- seq(1:(TimeSteps+1))
df.Props[1,2:(S+1)]     <- c(table(Mat.S))/(dm*dm)

A <-  seq(.5,.5,length=S)
#A <-  seq(0,0,length=S)


set.seed(150)
Y <- rlnorm(S,mean=0,sd=.5)
names(Y) <- seq(1:S)

d <- 1
Dist.Rate <- .025
R <- .2

# decay JCE
v <- 10
# Dispersal 
#Az <- 10

#DispKern <- function(x){
#  2*pi*x / (pi*exp(Az)*(1 + (x^2)*exp(-Az))^2 )
#}

Ae <- 1
DispKern <- function(x){
  (2*pi*x /(2*pi*Ae^2))*exp(-x/Ae)
}



d_A <- 1/sqrt(.15)
Rad_D <- 10
Disp.Local <- integrate(DispKern, lower = 0, upper = d_A/2)$value
Vec.Disp <- sapply(1:Rad_D, function(x) (1/(4*x))* integrate(DispKern, lower = d_A*(x-.5), upper = d_A*(x+.5))$value )
Disp.Tail <- 1-integrate(DispKern, lower = 0, upper =  d_A*(Rad_D+.5))$value
Vec.Disp.Final <-c(Disp.Local,Vec.Disp)
Dists <- seq(1:(Rad_D+1))*d_A
#plot(Vec.Disp.Final~Dists)

Rad <- Rad_D*2 + 1
Dummy.Frame <- matrix(1,ncol=Rad,nrow=Rad)
Dummy.Vec   <- data.frame(which(Dummy.Frame==1, arr.ind=TRUE))
Dummy.Dist  <-  abs(Dummy.Vec$row  -  (Rad-1)/2 ) + abs(Dummy.Vec$col  -   (Rad-1)/2)    


Distance.Func <- function(FC,Alpha,M,Yb,Xb,Y,P){
  
  # Conspecific Dispersal
  CONS                     <-  data.frame(which(M==FC, arr.ind=TRUE))
  DISTCONS                 <-  abs(CONS$row  - Xb) + abs(CONS$col  - Yb) 
  DISP.Seed                <-  sum(as.vector(table(factor(DISTCONS,levels=0:Rad_D))) * Vec.Disp.Final) + P*Disp.Tail
  
  # CNDD
  PredCons                 <- Alpha* (  sum(exp(-DISTCONS*d_A/v)) )^b
  
  # Seeds 
  
  Seeds    <- Y*DISP.Seed*exp(-PredCons)
  
  return(Seeds)
  
}

b <- 2.5


for(mm in 1:TimeSteps){
  
  Mat.S2 <- Mat.S  
  P      <- c(table(factor(Mat.S2, levels = 1:S)))/(dm*dm)   # P goes to proportion of each species in environment
  
  Prob.Dist                             <- matrix(runif(dm*dm),ncol=dm) # Matrix that defines hte probability that each species is disurbed
  Prob.Dist[Prob.Dist >= Dist.Rate]     <- NA # Spcies with draw less than Dist.Rate are distriubed 
  
  df.Rep                      <- which(!is.na(Prob.Dist), arr.ind=TRUE) # saves the indexes of each location that is distrubed
  
  x.val  <- df.Rep[,1]   # x coordinates of disturbance
  y.val  <- df.Rep[,2]   # y coordinates of distriubance
  
  Replaceb <- length(x.val)  # total number of distrubances
  
  Replacements <- sapply(1:Replaceb, function(x){  # function that determines which speices replaces disturbed patches (apply function loops over all distrubed)
    
    
    #Rad <- 21
    Rad <- Rad_D*2 + 1
    M   <- mat.torus(Mat.S2,Rad,x.val[x],y.val[x])
    Xb   <- (Rad+1)/2
    Yb   <- (Rad+1)/2
    
    # Predation <- sapply(1:S,function(z) Distance.Func(x.val[x],y.val[x],dm,Mat.S2,z,A[z],R,D_A,Rad))
    #Seeds <- rep(0,length=S)
    #Seeds            <-   sapply(1:S,function(z) Distance.Func(z,A[z],M,M2,Yb,Xb,Y[z],P[z]) )
    #Seeds            <-   sapply(1:S,function(z) Distance.Func.Rec(z,A[z],M,M2,Yb,Xb,Y[z],P[z]) )
    Seeds             <-   sapply(1:S,function(z) Distance.Func(z,A[z],M,Yb,Xb,Y[z],P[z]) )
    
    
    Total.Seeds         <- sum(Seeds) # total scaled number of seeds in local patch
    Vec.Probs           <- Seeds/Total.Seeds  # probability that each species wins lottery         
    
    Vec.Probs.Ordred    <- Vec.Probs[order(as.numeric(names(Vec.Probs)))]  # order previous probability values (1 to S)
    
    Vec.Sum             <- cumsum(Vec.Probs.Ordred) # creates probability intervals to determine which species wins 
    
    prob.rep <- runif(1) # draw from uniform distribution to determine the winner of the lottery
    
    Replacement <- as.numeric(names(Vec.Sum[min(which(Vec.Sum > prob.rep))])) # store winner of lottery
    return(Replacement) # return winner
  }
  )
  
  Mat.S[df.Rep] <- Replacements # put winner of lottery into correct location
  df.Props[mm+1,2:(S+1)]     <- c(table(factor(Mat.S, levels = 1:S)))/(dm*dm) # Store proportion of each species at each time step
  
} # Code that runs the simulation (notes inside)





df.PropsM <- as.matrix(df.Props)


write.csv(df.PropsM,"TS_Super_v10_D_1_b_2p5b.csv",quote=F,row.names=F)
write.csv(Mat.S,"DIST_Super_v10_D_1_b_2p5b.csv",quote=F,row.names=F)






