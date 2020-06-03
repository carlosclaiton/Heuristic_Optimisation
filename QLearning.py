def QLearning(NumAgents,NUM_ITERATIONS,CostFunction,noV):
    """function that calculate the optimum using Q-learning algoritym"""
    import numpy as np
    import random
    import math
    
    # start condition for the decision variables
    xVar = np.random.randint(2, size=[noV,NumAgents])

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
        fitness = np.ones(NumAgents)*math.nan
        for i in range(NumAgents):
            fitness[i] = 1/CostFunction(xVar[:,i])

        Best[IterNumb] = max(fitness) #find the agent with the better fitness
        indexB = (-fitness).argsort()
    
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

        indexQ = (-Q[:,IterNumb]).argsort()
        for i in range(NumAgents):
            if i == indexQ[0]:
                xVar[:,i] = xVar[:,i]
            elif i == indexQ[1] or i == indexQ[2]:
                indFlip = np.random.randint(noV)
                xVar[indFlip,i] = 1-xVar[indFlip,i]
            else:
                indFlip = np.random.randint(noV)
                xVar[:,i] = xVar[:,indexQ[0]]
                xVar[indFlip,i] = 1-xVar[indFlip,i]
        
        if Best[IterNumb-1] == Best[IterNumb]:
            stepEqual += 1
        else:
            stepEqual = 0
        if stepEqual == NumbConv:
            break

    BestxVar = xVar[:,indexQ[0]]
    BestCost = CostFunction(BestxVar)
     
    return(BestxVar,BestCost,1/Best,IterNumb)