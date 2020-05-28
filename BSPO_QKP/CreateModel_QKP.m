%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates the model to solve the Quadratic Knapsack problem %
% using a Particle Swarm algorithm, BPSO coded by Seyedali Mirjalili.     %
%                                                                         %
% Data used here is from Priority List - Revisited  by N. Order (2009)    %
%                                                                         %
% Developer: Carlos Kuhn                                                  %
% Contact Info: carlos.kuhn@dst.defence.gov.au                            %
% Created 27/03/2020                                                      %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model=CreateModel_QKP(Budget,Dependence)

%     Initiative = [{'O'},{'M'},{'N'},{'D'},{'J'},{'I'},{'B'},{'L'},{'E'},{'A'},...
%     {'H'},{'G'},{'C'},{'K'},{'F'}]';                        % Capabilities 
%     Cost = [1,1,3,4,7,6,8,11,13,12,14,17,14,20,16];         % cost of each capability (N Order, 2009)
%     ValueCap = [10,7,6,8,4,8,8,10,8,7,1,9,5,6,1]';          % value of each capability (N Order, 2009)
    Initiative = [{'P1'},{'P2'},{'P3'},{'P4'},{'P5'},{'P6'},{'P7'},{'P8'},{'P9'},{'P10'},...
    {'P11'},{'P12'},{'P13'},{'P14'},{'P15'}]';  
    Cost = [34.5,3,16.3,32.8,44.3,13.86,4.71,4.65,10.5,3.3,10.2,30.6,15.1,16.46,4.8];         % cost of each capability (Greiner et al (2003))
    ValueCap = [0.073,0.069,0.07,0.092,0.07,0.066,0.056,0.055,0.052,0.068,0.062,0.072,0.069,0.065,0.061]';          % value of each capability (Greiner et al (2003))

%     [Cost,IndexCost] = sort(Cost);
%     Initiative =Initiative(IndexCost);
%     ValueCap=ValueCap(IndexCost);

%% Dependence Matrices
    DepLevel = [0,0.25,0.5,0.75,1];  %Dependence Levels use DepLevel(5) -> Independent
                                     %                      DepLevel(4) -> Tangential
                                     %                      DepLevel(3) -> Associated
                                     %                      DepLevel(2) -> Dependente
                                     %                      DepLevel(1) -> Mandatory

%%DEPENDENCE MATRICES

   %Populates the dependence matrices with random integers
%     for i = 1:size(DepMatrix,1)
%         for j  = 1:size(DepMatrix,2)
%             DepMatrix(i,j)=randi([0,5]);
%         end
%     end

    DepMatrix = zeros(length(Initiative));
    alpha = ones(1,length(Initiative));
    
    beta = zeros(length(Initiative));
    for i = 1:length(DepMatrix)
        beta(:,i) = (1-alpha)'.*DepMatrix(:,i)./sum(DepMatrix,2);
    end
%Replace NaN for Zeros
    beta(isnan(beta)) = 0;
    
%%Create the 4 pages nxn Beta Matrix
    BetaMatrices = beta;
    AlphaMatrix = alpha';
    
%%THIS IS WHAT THE MODEL WILL PASS THROUGH
    model.Budget = Budget;
    model.Initiative = Initiative;    
    model.Cost=Cost;
    model.ValueCap =ValueCap;
    model.BetaMatrices=BetaMatrices;
    model.AlphaMatrix = AlphaMatrix;
    model.Dependence = Dependence;


end