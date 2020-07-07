def PSEQEA(noP,Max_iteration,CostFunction,noV):
    """function that calculate the optimum using Particle Swarm Embeded Quantum Evolutionary algoritym"""
    import numpy as np
    import random
    import math
    
    ############# weight variable #################
    c1=2       # cognitive weight
    c2=2       # social weight

    ########## Angle for the Quantum operator ################
    Theta = np.zeros([noP,noV])  #angle will replace Velocity vector when using QEA
  
    ############# initial q-bits all with the same propability ###########
    qbit = np.ones([2,noP,noV])/np.sqrt(2)

    ############# initialize Position vector ####################
    Position = np.random.randint(2, size=[noP,noV])

    ############# Cognitive component #############
    pBestScore = np.ones(noP)*math.inf
    pBest = np.zeros([noP,noV])

    ######### Social Component ################
    gBestScore = math.inf         #initial condition is infinity
    gBest = np.zeros(noV)

    ########### Q-bit opertor ##################
    def OpQ(theta):
        return (np.cos(theta), -np.sin(theta), np.sin(theta), np.cos(theta))
  
    #Calculate cost for each particle

    NumbConv=200             # Number of equal fitness
    ConvergenceCurve = np.ones(Max_iteration)*math.nan #Convergence vector
    IterNumb=0               # counter
    stepEqual = 0

    while IterNumb<Max_iteration:
        IterNumb+=1
        for i in range(noP): 
            fitness = CostFunction(Position[i,:])
        
            if pBestScore[i]>fitness:
                pBestScore[i]=fitness;
                pBest[i,:] = Position[i,:].copy()
            if gBestScore>fitness:
                gBestScore=fitness
                gBest=Position[i].copy()
        
        #Calculating the rotation angle for each q-bit
        for i in range(noP):
            for j in range(noV):
                Theta[i,j]= c1*(pBest[i,j]-Position[i,j])+c2*(gBest[j]-Position[i,j])

         # updating the q-bits
        for i in range(noP):
            for j in range(noV):
                Qmat = OpQ(Theta[i,j])
                qbitTemp1 = qbit[0,i,j]
                qbitTemp2 = qbit[1,i,j]
                qbit[0,i,j]= Qmat[0]*qbitTemp1+Qmat[1]*qbitTemp2
                qbit[1,i,j]= Qmat[2]*qbitTemp1+Qmat[3]*qbitTemp2

        #updating the states
        for i in range(noP):
            for j in range(noV):
                if random.uniform(0,1) < qbit[1,i,j]**2:
                    Position[i,j]=1;
                else:
                    Position[i,j]=0;

        ConvergenceCurve[IterNumb]= gBestScore;
        if ConvergenceCurve[IterNumb-1] == gBestScore:
            stepEqual += 1
        else:
            stepEqual = 0
        if stepEqual == NumbConv:
            break

    return(gBest,gBestScore,ConvergenceCurve,IterNumb)