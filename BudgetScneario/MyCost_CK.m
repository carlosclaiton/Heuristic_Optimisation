%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates the cost function for the Budget Scenario Problem %
% using the model presented by Tailor - Algorithmic complexity of two     %
% defence budget problems (2015)                                          %
%                                                                         %
%                                                                         %
% Developer: Carlos Kuhn                                                  %
% Contact Info: carlos.kuhn@dst.defence.gov.au                            %
% Created 20/03/2020                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Opt, sol]=MyCost_CK(x,model)

%     model=CreateModelBidData_CK(1,1,15);
    Cost = model.Cost;  % get the cost for each initiative
    Budget = model.Budget;
    ScenarioMatrix = model.ScenarioMatrix;
    BetaMatrices = model.BetaMatrices;
    ScenarioWeights = model.ScenarioWeights;
    AlphaMatrix = model.AlphaMatrix;


%%X is a vector that contains 1 or 0 as decision variables
%     x = zeros(1,15);
%     x([1,7,15])=1;


%% CALCULATE THE BEST VALUE BY MONEY 
    OptBudget = abs(Budget - sum(Cost.*x)); 
    
%% GET THE MAXIMUM BASED ON SCENARIOS WEIGHT
    
%%Calculate the Beta Matrix (note each collum became the beta for scenarios
    BetaXalpha = ones(length(Cost),length(ScenarioWeights))*NaN;
    
    for i = 1:length(ScenarioWeights)
        BetaXalpha(:,i) = BetaMatrices(:,:,i)*x'+AlphaMatrix(:,i);
    end

    
    SceBeta = ScenarioMatrix.*BetaXalpha;
    OptTemp = ones(1,length(ScenarioWeights))*NaN;
    for i = 1:length(ScenarioWeights)
        OptTemp(i) = max(SceBeta(:,i).*x');
    end

% Find the optimum with the available initiatives     
    OptScore = OptTemp*ScenarioWeights;
    Opt = 1/(OptScore) + 1000*OptBudget;
    
%% RETURN SOLUCTION
    sol.OptBudget=OptBudget;
    sol.OptScore=OptScore;
    sol.Opt=Opt;
    
end