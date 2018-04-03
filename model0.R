require(deSolve)

#Two equation I-W model
lever=function(t, n, parameters){ 
  with(as.list(parameters),{
    
    I = n[1]
    W = n[2]
      
    dIdt = beta*W*(1-I) - (mu0+mu)*I
    dWdt = lambda*I - (delt0 + delt)*W
    
    return(list(c(dIdt, dWdt)))
  })
}

#parameter values here (currently pulled out of thin air)
parameters = c(beta = 0.1,    #[insert parameter source here] 
               dt = 0.05,     #[insert parameter source here]
               lambda = 0.3,  #[insert parameter source here]
               delt0 = 0.1,   #[insert parameter source here]
               delt = 0,      #[insert parameter source here]
               mu0 = 0.1,     #[insert parameter source here]
               mu = 0)        #[insert parameter source here]

#Estimate r0 from parameter values
r0 = (parameters['beta'] * parameters['lambda']) / 
      ((parameters['delt0'] + parameters['delt']) * (parameters['mu0'] + parameters['mu']))

#Run the model for a year with random starting values
nstart = c(W = 0.5, I = 0.5)

time = c(0:365)

run1 = as.data.frame(ode(nstart, time, lever, parameters))

eq.vals = run1[dim(run1)[1],c(2:3)]

op.fx = function(max.time,   #time to run the model
                 start.I,    #starting proportion of infected individuals
                 start.W,    #starting value of environmental reservoir
                 #delt,       #value for added treatment of environmental reservoir
                 #mu,         #value for added treatment of infected individuals  
                 par,        #vector with delt and mu values to use in optim
                 c1,         #cost/effort associated with human treatment
                 c2,         #cost/effort associated with human treatment
                 c3,         #cost/effort associated with environmental treatment
                 c4){        #cost/effort associated with environemntal treatment
  t = c(0:max.time)    #total time to run the model

  s = c(I = start.I, W = start.W)

  parameters['delt'] = par[1]  #reset delt parameter
  parameters['mu'] = par[2]      #reset mu parameter

  m1 = ode(s, t, lever, parameters)   #run the model

#objective function
  to.min = sum(m1[,2], 
               c1*par[1]*max.time, c2*max.time*par[1]^2, 
               c3*max.time*par[2], c4*max.time*par[2]^2)

  return(to.min)
}         
  
op.fx(max.time = 365, start.I = (1-1/r0), start.W = 2,
      par = c(0.0001, 0.0001), c1 = 1, c2 = 1, c3 = 1, c4 = 1) 

op1 = optim(par = c(0.0001, 0.0001), #optimize over two intervention parameters
            fn = op.fx,         #use the function which returns the objective function
            #other values to pass to op.fx
            max.time = 365, start.I = (1-1/r0), start.W = 2,
            c1 = 1, c2 = 1, c3 = 1, c4 = 1) #need to add controls on intervention parameters

time.tests = c(1,5,10,25,50)*365  #look at different time scales for optimization

par.fill = matrix(ncol = 2, nrow = length(time.tests))  #

for(t in 1:length(time.tests)){
  par.fill[t,] = optim(par = c(0.0001, 0.0001), #optimize over two intervention parameters
        fn = op.fx, #use the function which returns the objective function
        #other values to pass to op.fx
        max.time = time.tests[t], start.I = (1-1/r0), start.W = 2,
        c1 = 1, c2 = 1, c3 = 1, c4 = 1)$par
    print(t)
}