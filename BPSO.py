def BPSO(noP, Max_Iteration, BPSO_num, CostFunction, noV):
    """function that calculate the optimum using Bynary Particle Swarm"""
    import numpy as np
    import random
    import math
    
    #Initial Parameters for PSO

    wMax=0.9           # Max inirtia weight
    wMin=0.4           # Min inirtia weight
    c1=5               # cognitive weight
    c2=5               # social weight
    Vmax=6

    #Velocity vector
    Velocity = np.zeros([noP,noV])       
    
    #Position vector, Uniformly random O or 1 
    Position = np.random.randint(2, size=[noP,noV])
      
    ########## Cognitive Component #############
    pBestScore = np.ones(noP)*math.inf
    pBest = np.zeros([noP,noV])
    ######### Social Component ################
    gBestScore = math.inf         #initial condition is infinity
    gBest = np.zeros(noV)
    
    ################################
    NumbConv=200             # Number of equal fitness
    ConvergenceCurve = np.ones(Max_Iteration)*math.nan #Convergence vector
    IterNumb=0               # counter
    stepEqual = 0
    
    while IterNumb<Max_Iteration:
        IterNumb+=1
        for i in range(noP): 
            fitness = CostFunction(Position[i,:])
        
            if pBestScore[i]>fitness:
                pBestScore[i]=fitness;
                pBest[i,:] = Position[i,:]
            if gBestScore>fitness:
                gBestScore=fitness
                gBest=Position[i,:]
        # update the W of PSO (inertia weight)
        w = wMax-IterNumb*((wMax-wMin)/Max_Iteration)
            
        # Update the Velocity and Position of particles
        for i in range(noP):
            for j in range(noV):
                Velocity[i,j]= w*Velocity[i,j] + c1*random.uniform(0,1)*(pBest[i,j]-Position[i,j])+c2*random.uniform(0,1)*(gBest[j]-Position[i,j]);           
                if Velocity[i,j]>Vmax:
                    Velocity[i,j] = Vmax
                if Velocity[i,j]<-Vmax:
                    Velocity[i,j]=-Vmax
            
                if BPSO_num==1:
                    s=1/(1+math.exp(-2*Velocity[i,j])) #S1 transfer function
                            
                if BPSO_num==2:
                #s=1/(1+math.exp(-Velocity[i,j]))   %S2 transfer function   
                #s=1/(1+exp(-(noV/45)*Velocity[i,j]))  #S2 transfer function
                    s=1/(1+math.exp(-(0.013*noV+0.13)*Velocity[i,j]))  #S2 transfer function
                if BPSO_num==3:
                    s=1/(1+math.exp(-Velocity[i,j]/2)) #S3 transfer function              
                if BPSO_num==4:
                    s=1/(1+math.exp(-Velocity[i,j]/3))  #S4 transfer function   
            
                if BPSO_num<=4:  # S-shaped transfer functions
                    if random.uniform(0,1)<s:  # Equation (4) and (8)
                        Position[i,j]=1
                    else:
                        Position[i,j]=0
            
                if BPSO_num==5:
                    s=abs(math.erf(((math.sqrt(math.pi)/2)*Velocity[i,j]))) #V1 transfer function
           
                if BPSO_num==6:
                    s=abs(math.tanh(Velocity[i,j])) #V2 transfer function
           
                if BPSO_num==7:
                    s=abs(Velocity[i,j]/math.sqrt((1+Velocity[i,j]**2))) #V3 transfer function
           
                if BPSO_num==8:
                    s=abs((2/math.pi)*math.atan((math.pi/2)*Velocity[i,j])) #V4 transfer function (VPSO)       
                if BPSO_num>4 & BPSO_num<=8: #V-shaped transfer functions
                    if random.uniform(0,1)<s: # Equation (10)
                        Position[i,j]= 1 - Position[i,j] 

        ConvergenceCurve[IterNumb]= gBestScore;
        if ConvergenceCurve[IterNumb-1] == gBestScore:
            stepEqual += 1
        else:
            stepEqual = 0
        if stepEqual == NumbConv:
            break
        
        
    return(gBest,gBestScore,ConvergenceCurve,IterNumb)