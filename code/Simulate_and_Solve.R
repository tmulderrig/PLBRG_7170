simulate_and_solve<-function(BV,A,VarA,VarE){
  set.seed(NULL)
  #setting block effects
  blockEffects=c(100,140)
  #setting up full incidence matrix to simulate y
  Xf=matrix(0,60,32)
  #loop to form full incidence matrix (Xf)
  Xf[1:30,1]=1
  Xf[31:60,2]=1

  for(i in c(1:30)){
    Xf[i,(i+2)]=1
    Xf[i+30,(i+2)]=1

  }
  # generate y values - the variance for varieties = 1 and the residual variance = 4
  y=Xf[,1:2]%*%blockEffects+Xf[,3:32]%*%BV+rnorm(60,0,VarE**.5)

  # setting up incidence matrices
  X=matrix(0,60,31)
  #Add an overall mean and remove one level of the block effects
  X[,1]=1
  X[31:60,2]=1


  # set up X and Z matrix columns associate with the variety effect
  Z=Xf[,3:32]

  #for fixed effects we need to remove 1 level of the variety effect
  X[,3:31]=Xf[,4:32]

  # setting up LHS and RHS for OLS
  LHS_OLS=t(X)%*%X
  RHS_OLS=t(X)%*%y

  # setting up LHS and RHS for mixed models
  LHS_MM=matrix(0,32,32)
  #calculating alpha (residual variance/random effect variance)
  alpha=VarE/VarA
  LHS_MM[1:2,1:2] = t(X[,1:2])%*%X[,1:2]
  LHS_MM[1:2,3:32] = t(X[,1:2])%*%Z
  LHS_MM[3:32,1:2] = t(Z)%*%X[,1:2]
  #modifying to use the inverse of the A matrix
  LHS_MM[3:32,3:32] = t(Z)%*%Z + alpha*solve(A)
  # note that alpha is generally unknown (unless you simulated the data) so mixed model software has to estimate alpha
  RHS_MM=matrix(0,32,1)
  RHS_MM[1:2,1]=t(X[,1:2])%*%y
  RHS_MM[3:32,1]=t(Z)%*%y

  #solving for OLS and MM
  sol_OLS=solve(LHS_OLS)%*%RHS_OLS
  sol_MM=solve(LHS_MM)%*%RHS_MM

  # Plot OLS vs MM solutions
  OLS=rep(0,30)
  #Adding mean to OLS solutions for variety
  OLS[2:30]=rep(sol_OLS[1],29)+sol_OLS[3:31]
  OLS[1]=sol_OLS[1]

  #Adding mean to OLS solutions for variety
  MM=rep(sol_MM[1],30)+sol_MM[3:32]
  LSM=OLS
  sol=cbind(LSM,MM)
  return(sol)


}
