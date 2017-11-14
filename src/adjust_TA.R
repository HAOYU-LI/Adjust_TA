source("./discrepancyCriteria.R")
source("./fcts_optimalLHSDesigns.R")
source("./lhsDesign.R")

#####Function to extract elements from lists composed with lists#####
extract_list<-function(l){return(l[[1]])}
#####discrepADTA_LHS#####
#####L2 DISCREPANCY LHS VIA ESE OPTIMIZATION#####

#---------------------------------------------------------------------------|
#args :  design     : the design                                            |
#        T0    : the initial temperature                                    |
#        inner_it  : number of iterations for inner loop                    |
#        J     : number of new proposed LHS in inner loop                   |
#        it    : number of iterations for outer loop                        |
#        criterion: "C2", "W2" or "L2star"                               |
#output        : a list containing all the input arguments plus:            |
#       low L2_discrepancy design                                           |
#       vector of criterion values along the iterations                     |
#       vector of temperature values along the iterations                   |
#       vector of acceptation probability values along the iterations       |
#depends :  discrepancyL2_EP_ESE, discrepancyW2_EP_ESE, discrepancyL2_STAR  |
#           discrepancyC2_EP_ESE, discrepancyCriteria                       |
#---------------------------------------------------------------------------|

para_c = 0.03
discrepADTA_LHS<-function(design,T0=0.005*discrepancyCriteria(design,type='C2')[[1]],inner_it=100,J=50,it=2,criterion="C2")
{
  #
  #Adjusted TA algorithm is implemented based on code of DiceDesign
  #Only latin hypercube design is implemented.
  #
  #
  m<-design
  crit <- NULL ; temp <- NULL ; proba <- NULL
  
  ############################adjust_TA_Treshold_sequence#######################
  threshold_seq = c(T0)
  for(t in c(2:(it-1)))
  {
    T_cur = (it-t)/it*threshold_seq[t-1]
    if(T_cur <= para_c * T0)
    {
      T_cur = (it-t-1)/(it-t)*threshold_seq[t-1]
    }
    threshold_seq = c(threshold_seq,T_cur)
  }
  threshold_seq = c(threshold_seq,0)
  ############################adjust_TA_Treshold_sequence#######################
  
  if(criterion=="C2")
  {
    d<-ncol(m)
    # Temperature<-T0
    Best<-m
    dC2<-discrepancyCriteria(m,type='C2')[[1]]                   
    best<-dC2
    crit <- dC2
    
    for (q in 1:it)
    {
      
      BOLD<-Best
      bold<-best                                   # Best=new LHS built at every step        
      # BOLD= new LHS built at each iteration q
      ni<-0
      count<-0
      na<-0
      while(count<=inner_it)                      # inner loop
      {
        count<-count+1
        
        
        modulo<-count%%d                         # d : number of columns of m
        l<-list(m)
        l<-rep(l,J) #liste de liste
        
        g<-lapply(l,discrepancyC2_EP_ESE,k=modulo+1,p=discrepancyCriteria(m,type='C2')[[1]]^2)
        values<-lapply(g,extract_list) 
        k<-which.min(values)
        a<-values[[k]]  
        
        Delta<-a-dC2
        
        if((Delta)<=(threshold_seq[q]))# higher is the temperature, higher is the probability of accepting a bad design.
          # if Delta is low, the probability is high of accepting a bad design.   
          # if Delta>Temperature, m is never accept.
          
          
        {m<-g[[k]][[2]]
        
        dC2<-a
        
        na<-na+1
        if(a<=best)  
        {Best<-m
        best<-a
        ni<-ni+1}                       #if optimization is ok, ni=ni+1
        }
        crit <- c(crit,best)
      }
      
      v1<-na/inner_it    # v1<-acceptance ratio
      v2<-ni/inner_it    # v2<-optimization ratio
      
      temp <- c(temp,rep(threshold_seq[q],inner_it)) ; proba <- c(proba,rep(v1,inner_it))   
      
      # if (best-bold<0){f<-1
      # if(v1>=0.1 & v2<=v1)
      # {Temperature<-0.8*Temperature}
      # else {if (v1>=0.1 & v2==v1){}
      #   else {Temperature<-Temperature/0.8}
      # }
      # }
      # 
      # # if the criteria is optimized after the inner loop, then f equals 1
      # else {f<-0
      # if (v1<=0.1){Temperature<-Temperature/0.7}
      # else {Temperature<-Temperature*0.9}
      # }
      
      # else, it is the exploratory step
      
    }
    List.res <- list(design,T0,inner_it,J,it,criterion,Best,crit,temp,proba)
    names(List.res) <-c("InitialDesign","TO","inner_it","J","it","criterion","design","critValues","tempValues","probaValues") 
    
  }

  return(List.res)
}


lhsDesign <- function(n, dimension, randomized=TRUE, seed=NULL){
  
  # arguments : n         = number of points
  #             dimension = number of variables
  #  	    randomized = logical for randomized or centered points
  #             seed       = value of the random seed
  # output : a list containing the arguments and the LHS design
  
  
  # if no seed is provided in argument, choice of the seed for 'runif' and 'sample'
  if (is.null(seed)){
    seed <- as.numeric(Sys.time())
  }
  set.seed(seed)
  
  # Randomized LHS: U[0,1]-sampling of n x dim values
  if (randomized) ran = matrix(runif(n*dimension),nrow=n,ncol=dimension) 
  # Centered LHS
  else ran = matrix(0.5,nrow=n,ncol=dimension) 
  
  x = matrix(0,nrow=n,ncol=dimension)  # initializing matrix x
  
  for (i in 1:dimension) {
    idx = sample(1:n)        # vector of permutations of [1 to n]
    P = (idx-ran[,i]) / n    # vecteur of probabilities
    x[,i] <- P  }
  
  # Outputs:
  return(list(n=n,dimension=dimension,design=x,randomized=randomized,seed=seed))
}



