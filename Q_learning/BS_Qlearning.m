% This code uses the Bellmann Equation to determine the action of the
% agent. We have 10 agents that interact to eacho other to decide what is
% the best solucion subset to maximize the fiteness. 
%
% Developer: Carlos C. N. Kuhn
% email: carlosclaitonkuhn@gmail.com
% Created at 08/04/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sample usage

% addpath(genpath('/Users/carlosclaitonnoschangkuhn/Google Drive/MatLabCodes/OptiizationAlgorithm/CK_Implementations/BudgetScneario/'))
% savepath
% 
% clc
% clear all

% DepLevel = 1;          % Dependence Levels use DepLevel = 5 -> Independent
%                         %                      DepLevel = 4 -> Tangential
%                         %                      DepLevel = 3 -> Associated
%                         %                      DepLevel = 2 -> Dependente
%                         %                      DepLevel = 1 -> Mandatory
% 
% Budget = 20;            % Budget
% model=CreateModel_CK(Budget,DepLevel); % here is where the model is stablished 
% noV = length(model.Cost); % number of decision variables, match with the size of cost vector
% CostFunction=@(x) MyCost_CK(x,model); % Modify or replace Mycost.m according to your cost funciton
% NUM_ITERATIONS = 50000;
% NumAgents = 10; %number of agents to learn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BestxVar,BestCost,Best,IterNumb]=BS_Qlearning(NumAgents,NUM_ITERATIONS,CostFunction,noV)

% start condition for the decision variables
xVar = randi([0,1],[noV,NumAgents]);

% learning rate settings
alpha = 0.8;%0.8; 
gamma = 0.2;%0.5;

% build a state action matrix by finding all valid states from maze
% we have two actions (flip state or stay, 0 or 1);
Q = zeros(NumAgents,NUM_ITERATIONS);
Best = ones(NUM_ITERATIONS,1)*NaN;

NumbConv = 200;
IterNumb=1;
stepEqual = 0;


while IterNumb < NUM_ITERATIONS % && stepEqual < NumbConv
    IterNumb=IterNumb+1;
    %%calculate the fitness for each xVar
    fitness = ones(size(xVar,2),1)*NaN;
    for i=1:NumAgents
        fitness(i) = 1/CostFunction(xVar(:,i)');
    end
    [Best(IterNumb),~] = max(fitness);%find the agent with the better fitnees
    [~,indexB] = sort(fitness,'descend');
    
    %%Determine the number of steps at the same optimum
%     stepEqual = sum(Best == Best(IterNumb));
    %Getting the rewarding
    for i=1:NumAgents
        if i == indexB(1)
            rewardVal = 2;
            Q(i,IterNumb) = Q(i,IterNumb) + alpha*(rewardVal+gamma*max(Q(:,IterNumb)) - Q(i,IterNumb));
        elseif i == indexB(2) || i==indexB(3)
            rewardVal = 1;
            Q(i,IterNumb) = Q(i,IterNumb) + alpha*(rewardVal+gamma*max(Q(:,IterNumb)) - Q(i,IterNumb));
        else
            rewardVal = 0;
            Q(i,IterNumb) = Q(i,IterNumb) + alpha*(rewardVal+gamma*max(Q(:,IterNumb)) - Q(i,IterNumb));
        end
    end
    
%     [~,indexQ] = max(Q(:,j)) ;          %The action will be select a random position to flipp
    [~,indexQ] = sort(Q(:,IterNumb),'descend');
    for i=1:NumAgents
        if i == indexQ(1)
            xVar(:,i) = xVar(:,i);
        elseif i == indexQ(2) || i == indexQ(3)
            indFlip = randi([1,noV]);
            xVar(indFlip,i) = 1-xVar(indFlip,i);
        else
            indFlip = randi([1,noV]);
%             indFlip = randi([1,noV],[1,floor(noV/6)]);
            xVar(:,i) = xVar(:,indexQ(1));
            xVar(indFlip,i) = 1-xVar(indFlip,i);
        end   
    end
    
    if Best(IterNumb) - Best(IterNumb-1) < 10^-9
        stepEqual = stepEqual+ 1;
    else
        stepEqual = 0;
    end
    
    if stepEqual > NumbConv % || IterNumb > Max_iteration-1
        break
    end
    
end

BestxVar = [xVar(:,indexQ(1))]';
BestCost = 1/CostFunction(BestxVar);


% figure()
% p1=plot(Best);
% p1.LineStyle  = '-'; p1.LineWidth = 0.5; p1.Marker = 'o'; 
% p1.MarkerSize = 8; p1.MarkerFaceColor = 'b';
% grid on
% ylabel('best score ','FontSize', 16)
% xlabel('iteration ','FontSize', 16)
% title('Convergence Curve', 'FontSize',18)

end


