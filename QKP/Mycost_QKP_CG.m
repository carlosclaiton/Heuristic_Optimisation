%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates the cost function for Cubic Graph problem         %
%                                                                         %
%                                                                         %
% Developer: Carlos Kuhn                                                  %
% Contact Info: carlos.kuhn@dst.defence.gov.au                            %
% Created 30/04/2020                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Opt, sol]=Mycost_QKP_CG(x,model)

    Cost = model.Cost;  % get the cost for each initiative
    Budget = model.Budget;
    Value = model.Value;
    Edge = model.Edge;

% X is a vector that contains 1 or 0 as decision variables
%  x = ones(1,15);
% %  x([6,10])=1; Debuging features
%  x([1,7,15])=1;


%%CALCULATE THE BEST VALUE BY MONEY 
    OptBudget = abs(Budget - sum(Cost.*x')); 
 
%%CALCULATE THE PROFIT
    Profit = sum(Value.*x');

%%CALCULATE THE PROFIT IF BOTH ITEM i AND j ARE ADDED.
    Profit2=0;    
    for i=1:length(Cost)
        for j=1:length(Cost)
            Profit2 = Profit2+Edge(i,j)*x(i)*x(j);
        end
    end
            
            
%%Calculate the Beta Matrix (note each collum became the beta for scenarios
    OptScore = Profit+Profit2/2;
    Opt = 1/(OptScore) + 100*OptBudget;
%     Opt = (1-OptScore) + 100*OptBudget;

    
%% RETURN SOLUCTION
    sol.OptBudget=OptBudget;
    sol.OptScore=OptScore;
    sol.Opt=Opt;
    
end
