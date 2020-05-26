class ModelBS:
  
    def __init__(self, DepLevel, NumberInitiave, Budget):
        
        import numpy as np
        import math
        import random
        
        #Level of the dependence between Initiatives
        DepLevelArray = np.array([0,0.25,0.5,0.75,1]); 
        
         #Initiatives 
        Initiative = np.array(['O','M','N','D','J','I','B','L','E','A','H','G','C','K','F']);
        
        nameload = 'data' +str(NumberInitiave)+'in.npy'
        data = np.load(nameload)
        Cost = data[0]
        S1 = data[1]
        S2 = data[2]
        S3 = data[3]
        S4 = data[4]
        
        #Creat the dependence matrices
        
        # Scenario 1
        alphaS1 = np.ones(NumberInitiave)
        alphaS1[[0,6]] = DepLevelArray[DepLevel];      # O and B
        
        DepMatrixS1 = np.zeros([NumberInitiave,NumberInitiave])
        DepMatrixS1[0,13] = 1;      # O needs K
        DepMatrixS1[6,13] = 1;      # B needs K
        
        betaS1 = np.zeros([NumberInitiave,NumberInitiave]);
        totalRow = np.sum(DepMatrixS1,axis=1)
        for i in range(NumberInitiave):
            for j in range(NumberInitiave):
                if totalRow[i] == 0:
                    betaS1[i,j] = (1-alphaS1[i])*DepMatrixS1[i,j]
                else:
                    betaS1[i,j] = (1-alphaS1[i])*DepMatrixS1[i,j]/totalRow[i]
       
        # Scenario 2
        alphaS2 = np.ones(NumberInitiave)

        DepMatrixS2 = np.zeros([NumberInitiave,NumberInitiave]);    # No dependences

        betaS2 = np.zeros([NumberInitiave,NumberInitiave]);
        totalRow = np.sum(DepMatrixS2,axis=1)
        for i in range(NumberInitiave):
            for j in range(NumberInitiave):
                if totalRow[i] == 0:
                    betaS2[i,j] = (1-alphaS2[i])*DepMatrixS2[i,j]
                else:
                    betaS2[i,j] = (1-alphaS2[i])*DepMatrixS2[i,j]/totalRow[i]
                
        # Scenario 3
        alphaS3 = np.ones(NumberInitiave)
        alphaS3[[4,7,10]] = DepLevelArray[DepLevel];      # J, L, H are dependent initiatives
        
        DepMatrixS3 = np.zeros([NumberInitiave,NumberInitiave]);
        DepMatrixS3[4,[6,8,11]] = 1;      # J needs B,E,G
        DepMatrixS3[7,[6,8,11]] = 1;      # L needs B,E,G
        DepMatrixS3[10,[6,8,11]] = 1;     # H needs B,E,G
        
        betaS3 = np.zeros([NumberInitiave,NumberInitiave]);
        totalRow = np.sum(DepMatrixS3,axis=1)
        for i in range(NumberInitiave):
            for j in range(NumberInitiave):
                if totalRow[i] == 0:
                    betaS3[i,j] = (1-alphaS3[i])*DepMatrixS3[i,j]
                else:
                     betaS3[i,j] = (1-alphaS3[i])*DepMatrixS3[i,j]/totalRow[i]

        
        # Scenario 4
        alphaS4 = np.ones(NumberInitiave)

        DepMatrixS4 = np.zeros([NumberInitiave,NumberInitiave]);    # No dependences

        betaS4 = np.zeros([NumberInitiave,NumberInitiave]);
        totalRow = np.sum(DepMatrixS4,axis=1)
        for i in range(NumberInitiave):
            for j in range(NumberInitiave):
                if totalRow[i] == 0:
                    betaS4[i,j] = (1-alphaS4[i])*DepMatrixS4[i,j]
                else:
                    betaS4[i,j] = (1-alphaS4[i])*DepMatrixS4[i,j]/totalRow[i]
        
        
        BetaMatrices = np.array([betaS1,betaS2,betaS3,betaS4])
        AlphaMatrix = np.array([alphaS1,alphaS2,alphaS3,alphaS4])
        
        self.BetaMatrices = BetaMatrices
        self.Initiative = Initiative
        self.Cost = Cost
        self.AlphaMatrix = AlphaMatrix
        self.DepLevel = DepLevel
        self.NumberInitiave = NumberInitiave   
        # Probability of a scneario
        self.ScenarioWeights = np.array([0.3,0.1,0.4,0.2]);
        self.ScenarioMatrix = np.transpose(np.array([S1,S2,S3,S4]))
        self.Budget=Budget
 
    # Define the method (function)  that calculates the cost
    def MyCost(self,x):
        
        import numpy as np
        
        # Calculate the value expended to purchase
        OptBudget = np.abs(self.Budget - np.sum(self.Cost*x)); 
        
        Beta_alpha = np.zeros([self.NumberInitiave,len(self.ScenarioWeights)])
        #Define the value using the Beta Matrix
        for i in range(len(self.ScenarioWeights)):  
            Beta_alpha[:,i] = self.BetaMatrices[i,:,:].dot(x) + self.AlphaMatrix[i,:]
         
        
        SceBeta = self.ScenarioMatrix*Beta_alpha;
        
        OptTemp = np.ones(len(self.ScenarioWeights))
        for i in range(len(self.ScenarioWeights)):
            OptTemp[i] = max(SceBeta[:,i]*x)
        
        # Calculate the Score 
        OptScore = OptTemp.dot(self.ScenarioWeights);
        if OptScore == 0.:
            Opt = 1/(OptScore+0.00000001) + 1000*OptBudget
        else:
            Opt = 1/OptScore + 1000*OptBudget
        
        return(Opt)
        