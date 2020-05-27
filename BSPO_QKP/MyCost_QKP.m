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

function [Opt, sol]=MyCost_QKP(x,model)

    Cost = model.Cost;  % get the cost for each initiative
    Budget = model.Budget;
    ValueCap = model.ValueCap;
    BetaMatrices = model.BetaMatrices;
    AlphaMatrix = model.AlphaMatrix;

% X is a vector that contains 1 or 0 as decision variables
%  x = ones(1,15);
% %  x([6,10])=1; Debuging features
%  x([1,7,15])=1;


%% CALCULATE THE BEST VALUE BY MONEY 
    OptBudget = abs(Budget - sum(Cost.*x)); 
    
%% GET THE MAXIMUM BASED ON SCENARIOS WEIGHT
    
%%Calculate the Beta Matrix (note each collum became the beta for scenarios
    BetaXalpha = BetaMatrices*x'+AlphaMatrix;

    SceBeta = ValueCap.*BetaXalpha;
    OptScore = sum(SceBeta.*x');
    Opt = 1/(OptScore) + 5*OptBudget;
%     Opt = (1-OptScore) + 100*OptBudget;

    
%% RETURN SOLUCTION
    sol.OptBudget=OptBudget;
    sol.OptScore=OptScore;
    sol.Opt=Opt;
    
end