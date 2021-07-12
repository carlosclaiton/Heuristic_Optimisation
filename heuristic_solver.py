"""
Author: Carlos C. N. Kuhn
email: carlos.kuhn@dst.defence.gov.au
Last Modified: 10 Jul 2021
"""

import numpy as np
import random
import math

class HeuristicSolver:
    """
    This class contains the implementations for differents ways to compute the best portfolio:
    It will use a fitness fucntion to be solved by stochastic and reinforcement learning algorithms.
    """
    
    def __init__(self):
        """
        Initialize

        """
        return

    def bpso_solver(noP, Max_Iteration, BPSO_num, CostFunction, noV, Initialguess):
        """
        function that calculate the optimum using Bynary Particle Swarm
        """
        
        #Initial Parameters for PSO
        wMax=0.9           # Max inirtia weight
        wMin=0.4           # Min inirtia weight
        c1=2               # cognitive weight
        c2=3               # social weight
        Vmax=6

        #Velocity vector
        Velocity = np.zeros([noP,noV])       

        #Position vector, Uniformly random O or 1 
        Position = np.random.randint(2, size=[noP,noV])
        Position[0]=Initialguess

        ########## Cognitive Component #############
        pBestScore = - np.ones(noP)*math.inf
        pBest = np.zeros([noP,noV])
        ######### Social Component ################
        gBestScore = - math.inf         #initial condition is infinity
        gBest = np.zeros(noV)

        ################################
        NumbConv=200             # Number of equal fitness
        ConvergenceCurve = np.ones(Max_Iteration)*math.nan #Convergence vector
        IterNumb=0               # counter
        stepEqual = 0

        while IterNumb<(Max_Iteration-1):
            IterNumb+=1
            for i in range(noP): 
                fitness = CostFunction(Position[i,:])[0]

                if pBestScore[i] < fitness:
                    pBestScore[i]=fitness;
                    pBest[i,:] = Position[i,:].copy()
                if gBestScore < fitness:
                    gBestScore = fitness
                    gBest = Position[i].copy()

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
                        s=1/(1+math.exp(-Velocity[i,j]))   #S2 transfer function   
                        #s=1/(1+exp(-(noV/45)*Velocity[i,j]))  #S2 transfer function
                        #s=1/(1+math.exp(-(0.013*noV+0.13)*Velocity[i,j]))  #S2 transfer function, work for the Knapsack problem
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


    def pseqea_solver(noP,Max_iteration,CostFunction,noV,Initialguess):
        """
        function that calculate the optimum using Particle Swarm Embeded Quantum Evolutionary algoritym
        """

        ############# weight variable #################
        c1=2       # cognitive weight
        c2=2       # social weight

        ########## Angle for the Quantum operator ################
        Theta = np.zeros([noP,noV])  #angle will replace Velocity vector when using QEA

        ############# initial q-bits all with the same propability ###########
        qbit = np.ones([2,noP,noV])/np.sqrt(2)

        ############# initialize Position vector ####################
        Position = np.random.randint(2, size=[noP,noV])
        Position[0]=Initialguess

        ############# Cognitive component #############
        pBestScore = -np.ones(noP)*math.inf
        pBest = np.zeros([noP,noV])

        ######### Social Component ################
        gBestScore = -math.inf         #initial condition is infinity
        gBest = np.zeros(noV)

        ########### Q-bit opertor ##################
        def OpQ(theta):
            return (np.cos(theta), -np.sin(theta), np.sin(theta), np.cos(theta))

        #Calculate cost for each particle

        NumbConv=200             # Number of equal fitness
        ConvergenceCurve = np.ones(Max_iteration)*math.nan #Convergence vector
        IterNumb=0               # counter
        stepEqual = 0

        while IterNumb<(Max_iteration-1):
            IterNumb+=1
            for i in range(noP): 
                fitness = CostFunction(Position[i,:])[0]

                if pBestScore[i]<fitness:
                    pBestScore[i]=fitness;
                    pBest[i,:] = Position[i,:].copy()
                if gBestScore<fitness:
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
    
    def q_learning_solver(NumAgents,NUM_ITERATIONS,CostFunction,noV,Initialguess):
        """
        function that calculate the optimum using Q-learning reinforcement learning algoritym
        """
        # start condition for the decision variables
        xVar = np.random.randint(2, size=[NumAgents,noV])
        xVar[0]=Initialguess

        # learning rate settings
        alpha = 0.8   #0.8; 
        gamma = 0.2   #0.5;

        # build a state action matrix by finding all valid states from maze
        # we have two actions (flip state or stay, 0 or 1);
        Q = np.zeros([NumAgents,NUM_ITERATIONS])
        Best = np.ones(NUM_ITERATIONS)*math.nan

        NumbConv = 200
        IterNumb=0
        stepEqual = 0

        while IterNumb < (NUM_ITERATIONS-1): 
            IterNumb += 1
            # calculate the fitness for each xVar
            fitness = np.array([(CostFunction(x)[0]) for x in xVar])

            #find the agent with the better fitness
            Best[IterNumb] = max(fitness) 
            indexB = fitness.argsort()[::-1]

            # Getting the rewarding
            for i in range(NumAgents):
                if i == indexB[0]:
                    rewardVal = 2
                    Q[i,IterNumb] = Q[i,IterNumb] + alpha*(rewardVal+gamma*max(Q[:,IterNumb]) - Q[i,IterNumb])
                elif i == indexB[1] or i==indexB[2]:
                    rewardVal = 1
                    Q[i,IterNumb] = Q[i,IterNumb] + alpha*(rewardVal+gamma*max(Q[:,IterNumb]) -Q[i,IterNumb])
                else:
                    rewardVal = 0;
                    Q[i,IterNumb] = Q[i,IterNumb] + alpha*(rewardVal+gamma*max(Q[:,IterNumb]) -Q[i,IterNumb])

            #indexQ = (-Q[:,IterNumb]).argsort()
            indexQ = Q[:,IterNumb].argsort()[::-1]
            for i in range(NumAgents):
                if i == indexQ[0]:
                    xVar[i] = xVar[i]
                elif i == indexQ[1] or i == indexQ[2]:
                    indFlip = np.random.randint(noV)
                    xVar[i,indFlip] = 1-xVar[i,indFlip]
                else:
                    indFlip = np.random.randint(noV)
                    xVar[i] = xVar[indexQ[0]]
                    xVar[i,indFlip] = 1-xVar[i,indFlip]

            if Best[IterNumb-1] == Best[IterNumb]:
                stepEqual += 1
            else:
                stepEqual = 0
            if stepEqual == NumbConv:
                break

        BestxVar = xVar[indexQ[0]]
        BestCost = CostFunction(BestxVar)[0]

        return(BestxVar,BestCost,Best,IterNumb) 
   
