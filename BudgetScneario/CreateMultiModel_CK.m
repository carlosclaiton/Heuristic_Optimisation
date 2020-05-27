%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates the model to solve the Budget Scenario problem   %
% using a Particle Swarm algorithm, BPSO coded by Seyedali Mirjalili.    %
%                                                                        %
% Data used here is from Priority List - Revisited  by N. Order (2009)   %
%                                                                        %
% Developer: Carlos Kuhn                                                 %
% Contact Info: carlos.kuhn@dst.defence.gov.au                           %
% Created 20/03/2020                                                     %
%
% Adapted to deal with modifications on the model, the idea we feed the
% model with the number of the initiative we want to remove from the model,
% that way I can call this function inside of  loop on the main code.
%
% Carlos Kuhn - 23/0/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model=CreateMultiModel_CK(Budget,Dependence,indRemov)
    
   
    ScenarioWeights = [0.3,0.1,0.4,0.2]';                   % probability of scenario
    Initiative = [{'O'},{'M'},{'N'},{'D'},{'J'},{'I'},{'B'},{'L'},{'E'},{'A'},...
    {'H'},{'G'},{'C'},{'K'},{'F'}]';                        % initiatives 
    Cost = [1,1,3,4,7,6,8,11,13,12,14,17,14,20,16];         % cost of each initiative
    S1 = [10,7,6,8,4,8,8,10,8,7,1,9,5,6,1]';                % score for scenario 1
    S2 = [1,1,3,7,3,7,8,5,4,2,5,5,7,4,7]';                  % score for scenario 2
    S3 = [8,3,5,1,7,4,3,6,9,9,9,8,5,9,9]';                  % score for scenario 3
    S4 = [2,4,6,2,9,1,10,9,7,1,9,4,7,4,1]';                 % score for scenario 4

    [Cost,IndexCost] = sort(Cost);
    Initiative =Initiative(IndexCost);
    S1=S1(IndexCost);
    S2=S2(IndexCost);
    S3=S3(IndexCost); 
    S4=S4(IndexCost);
    
    ScenarioMatrix = [S1,S2,S3,S4];                         %arrange the scores in a Matrix v_ij
% 
%% Getting the Cost Matrix 
%    CostMat = ones(length(Cost)-1,length(Cost))*NaN;
%    InitiativeMatrix = cell(length(Cost),1);
%    for i = 1:length(Cost)
%        tempCost=Cost;
%        tempCost(i) = [];
%        CostMat(:,i) = tempCost';
%        
%        tempInitiative  = Initiative;
%        tempInitiative(i) = [];
%        InitiativeMatrix{i} = tempInitiative;
%    end
%  ???? Wrong approach, that way it will generate way to much data, 
%  ???? it would be good to paralel processing    
    
    
    
%% Dependence Matrices
    DepLevel = [0,0.25,0.5,0.75,1];  %Dependence Levels use DepLevel(5) -> Independent
                                     %                      DepLevel(4) -> Tangential
                                     %                      DepLevel(3) -> Associated
                                     %                      DepLevel(2) -> Dependente
                                     %                      DepLevel(1) -> Mandatory

%%DEPENDENCE MATRICES
    DepMatrix = zeros(length(Initiative)); %That is the D matrix Taylor Model

%% Scenario 1
    alphaS1 = ones(1,length(Initiative));
    alphaS1([1,7]) = DepLevel(Dependence);      % O and B

    DepMatrixS1 = DepMatrix;
    DepMatrixS1(1,15) = 1;      % O needs K
    DepMatrixS1(7,15) = 1;      % B needs K
    
    alphaS1(indRemov) = [];     % removing the initiative
    DepMatrixS1(indRemov,:) = []; % removing the coluum associate to the initiative
    DepMatrixS1(:,indRemov) = []; % removing the rows associate to the initiative
    
    betaS1 = zeros(length(alphaS1));
    for i = 1:length(DepMatrixS1)
        betaS1(:,i) = (1-alphaS1)'.*DepMatrixS1(:,i)./sum(DepMatrixS1,2);
    end
%Replace NaN for Zeros
    betaS1(isnan(betaS1)) = 0;
    
%%Scenario 2
    alphaS2 = ones(1,length(Initiative));
    DepMatrixS2 = DepMatrix;    % No dependence
    
    alphaS2(indRemov) = [];     % removing the initiative
    DepMatrixS2(indRemov,:) = []; % removing the coluum associate to the initiative
    DepMatrixS2(:,indRemov) = []; % removing the rows associate to the initiative
    
    betaS2 = zeros(length(alphaS2));
    for i = 1:length(DepMatrixS2)
        betaS2(:,i) = (1-alphaS2)'.*DepMatrixS2(:,i)./sum(DepMatrixS2,2);
    end
%Replace NaN for Zeros
    betaS2(isnan(betaS2)) = 0;
       
%%Scenario 3
    alphaS3 = ones(1,length(Initiative));
%     alphaS3([6,8,11]) = 0;
    alphaS3([6,8,11]) = DepLevel(Dependence);      % J, L, H are dependent initiatives

    DepMatrixS3 = DepMatrix;
    DepMatrixS3(6,[7,10,14]) = 1;      % J needs B,E,G
    DepMatrixS3(8,[7,10,14]) = 1;      % L needs B,E,G
    DepMatrixS3(11,[7,10,14]) = 1;     % H needs B,E,G 
    
    alphaS3(indRemov) = [];     % removing the initiative
    DepMatrixS3(indRemov,:) = []; % removing the coluum associate to the initiative
    DepMatrixS3(:,indRemov) = []; % removing the rows associate to the initiative
    
    betaS3 = zeros(length(alphaS3));
    for i = 1:length(DepMatrixS3)
        betaS3(:,i) = (1-alphaS3)'.*DepMatrixS3(:,i)./sum(DepMatrixS3,2);
    end
%Replace NaN for Zeros
    betaS3(isnan(betaS3)) = 0;
    
%%Scenario 4
    alphaS4 = ones(1,length(Initiative));
    DepMatrixS4 = DepMatrix;    % No dependence
    
    alphaS4(indRemov) = [];     % removing the initiative
    DepMatrixS4(indRemov,:) = []; % removing the coluum associate to the initiative
    DepMatrixS4(:,indRemov) = []; % removing the rows associate to the initiative
    
    betaS4 = zeros(length(alphaS4));
    for i = 1:length(DepMatrixS4)
        betaS4(:,i) = (1-alphaS4)'.*DepMatrixS4(:,i)./sum(DepMatrixS4,2);
    end
%Replace NaN for Zeros
    betaS4(isnan(betaS4)) = 0;                                       

%%Create the 4 pages nxn Beta Matrix
    BetaMatrices =  cat(3,betaS1,betaS2,betaS3,betaS4);
    AlphaMatrix = [alphaS1',alphaS2',alphaS3',alphaS4'];
    
    %% Removing the initiatives from the model
    Cost(indRemov)=[];
    Initiative(indRemov)=[];
    ScenarioMatrix(indRemov,:) = [];
 
    
%%THIS IS WHAT THE MODEL WILL PASS THROUGH
    model.Budget = Budget;
    model.Initiative = Initiative;    
    model.Cost=Cost;
    model.ScenarioMatrix =ScenarioMatrix;
    model.BetaMatrices=BetaMatrices;
    model.ScenarioWeights =ScenarioWeights;
    model.AlphaMatrix = AlphaMatrix;
    model.Dependence = Dependence;


end