%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS CODE SOLVE THE HYPOTHETICAL PROBLEM PROPOST BY N. ORDER (2007 AND
% 2009 PRIORITY LIST) USING A NUMERICAL APPROACH WITH DEPENDENCE 
% BETWEEN THE INITIATIES.
% MAIN IDEA IS TO UNDERSTAND THE PROBLEM.
% 
% History
% Created by Carlos Kuhn
% Created at 12/03/2020
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GENERATING THE TABLE PRESENTED ON THE PAPER

% CHOOSE THE BUDGET  FOR TEST MODE (COMMENT THE FUNCTION AND END OF
% FUCNTION)
% clear all;
% budget = 10;
function [OptSubFound,optimun] = PLMatrix_Dependence(budget)
ScenarioWeights = [0.3,0.1,0.4,0.2]';
Initiative = [{'O'},{'M'},{'N'},{'D'},{'J'},{'I'},{'B'},{'L'},{'E'},{'A'},...
    {'H'},{'G'},{'C'},{'K'},{'F'}]';
Cost = [1,1,3,4,7,6,8,11,13,12,14,17,14,20,16]';
S1 = [10,7,6,8,4,8,8,10,8,7,1,9,5,6,1]';
S2 = [1,1,3,7,3,7,8,5,4,2,5,5,7,4,7]';
S3 = [8,3,5,1,7,4,3,6,9,9,9,8,5,9,9]';
S4 = [2,4,6,2,9,1,10,9,7,1,9,4,7,4,1]';
% ScenarioMatrix = [S1,S2,S3,S4];
% Score = ScenarioMatrix*ScenarioWeights;
% ValueRatio = sum(Score,2)./Cost;
Rank = 1:15;
Rank = Rank';
% DataPriority = table(Initiative,Cost,S1,S2,S3,S4,Score,ValueRatio,Rank);

%% SORT IN ASCENT ORDER BY COST
[Cost,IndexCost] = sort(Cost);
Initiative =Initiative(IndexCost);
S1=S1(IndexCost);
S2=S2(IndexCost);
S3=S3(IndexCost); 
S4=S4(IndexCost);
ScenarioMatrix = [S1,S2,S3,S4];
Score = ScenarioMatrix*ScenarioWeights;
ValueRatio = sum(Score,2)./Cost;

Rank=Rank(IndexCost);
DataPriority = table(Initiative,Cost,S1,S2,S3,S4,Score,ValueRatio,Rank);
% disp(DataPriority)


%% Now to create the 3 dimensional matrix (three possible combinations)

% NumberComb = nchoosek(15,3);
sizetest = length(Cost(Cost<=budget));
CostMatrix3 = zeros(sizetest,sizetest,sizetest);
for i=1:sizetest
    for j = 1:sizetest
        for k = 1:sizetest
            if j==i && k==j
                CostMatrix3(i,j,k) = Cost(i);
            elseif i==k 
                CostMatrix3(i,j,k) = Cost(i)+Cost(j);
            elseif i==j || k==j
                CostMatrix3(i,j,k) = Cost(i)+Cost(k);
            else
                CostMatrix3(i,j,k) = Cost(i)+Cost(j)+Cost(k);
            end
        end
    end    
end

%% Approach using linear indeces 
% ind = find(CostMatrix3>(budget-0.5) & CostMatrix3<=budget);
ind = find(CostMatrix3<=budget);
sz = size(CostMatrix3);
[row,col,dim] = ind2sub(sz,ind);
%put back in the original matrix form to identify the combinations
combsSubsetIndeces = [row,col,dim];
Combs = cell(length(row),1);

for i=1:length(row)
    Combs{i} = {unique(combsSubsetIndeces(i,:))};
end
%% we need to create the depency table here
% O and B depend on K for Scenario 1 S1
% if Budget < Cost of K we need to remove it O  or B from the optimum score



TotalSubSetScore = ones(length(row),1)*NaN;
for  i= 1:length(row)
    if ismember(1,cell2mat(Combs{i})) && ismember(7,cell2mat(Combs{i}))&& ismember(15,cell2mat(Combs{i}))
        S1D = S1;
    elseif ismember(1,cell2mat(Combs{i})) && ismember(7,cell2mat(Combs{i})) %&& budget<Cost(15)
        S1D = [0,7,6,8,8,4,0,10,7,8,1,5,1,9,6]';
    elseif ismember(1,cell2mat(Combs{i})) %&& budget<Cost(15)
        S1D = [0,7,6,8,8,4,8,10,7,8,1,5,1,9,6]';
    elseif ismember(7,cell2mat(Combs{i})) %&& budget<Cost(15)
        S1D = [10,7,6,8,8,4,0,10,7,8,1,5,1,9,6]';
    else
        S1D = S1;
    end

    if ismember(7,cell2mat(Combs{i})) && ismember(10,cell2mat(Combs{i})) && ...
           ismember(14,cell2mat(Combs{i}))% && budget<sum(Cost([7,10,14]))
        S3D =  S3;
    elseif ismember(6,cell2mat(Combs{i})) && ismember(8,cell2mat(Combs{i})) && ...
            ismember(11,cell2mat(Combs{i}))% && budget<sum(Cost([7,10,14]))
        S3D =  [8,3,5,1,4,0,3,0,9,9,0,5,9,8,9]';
    elseif ismember(8,cell2mat(Combs{i})) && ...
            ismember(11,cell2mat(Combs{1}))% && budget<sum(Cost([7,10,14]))
        S3D =  [8,3,5,1,4,7,3,0,9,9,0,5,9,8,9]';
    elseif ismember(6,cell2mat(Combs{i})) && ...
            ismember(11,cell2mat(Combs{i})) %&& budget<sum(Cost([7,10,14]))
        S3D =  [8,3,5,1,4,0,3,6,9,9,0,5,9,8,9]';
    elseif ismember(6,cell2mat(Combs{i})) && ...
            ismember(8,cell2mat(Combs{i})) %&& budget<sum(Cost([7,10,14]))
        S3D =  [8,3,5,1,4,0,3,0,9,9,9,5,9,8,9]';    
    elseif ismember(6,cell2mat(Combs{i})) %&& budget<sum(Cost([7,10,14]))
        S3D =  [8,3,5,1,4,0,3,6,9,9,9,5,9,8,9]';
    elseif ismember(8,cell2mat(Combs{i})) %&& budget<sum(Cost([7,10,14]))
        S3D =  [8,3,5,1,4,7,3,0,9,9,9,5,9,8,9]';
    elseif ismember(11,cell2mat(Combs{i})) %&& budget<sum(Cost([7,10,14]))
        S3D =  [8,3,5,1,4,7,3,6,9,9,0,5,9,8,9]';
    else
        S3D = S3;
    end
    
    ScenarioMatrix = [S1D,S2,S3D,S4];
    SubScenario = ScenarioMatrix(Combs{i}{1},:);
    TotalSubSetScore(i) = max(SubScenario,[],1)*ScenarioWeights;
end

%% Selection the optimum
optimun = max(TotalSubSetScore);
OptSubIndex = find(TotalSubSetScore==optimun);
OptSub = Combs(OptSubIndex);
% OptSubFound = unique(cell2mat(OptSub));


OptSub2 = cellfun(@(x) [x{:}], OptSub, 'un',0); % convert to cell array of numeric array
[m, tf] = padcat(OptSub2{:}); % concatenate, pad rows with NaNs
m(~tf) = 0 ;% replace NaNs by zeros

OptSubFound = unique(m,'rows');

% disp([Initiative(OptSubFound)' ' Score : ' num2str(optimun) ]);

end


%% find the comb
% Combs2 = cellfun(@(x) [x{:}], Combs, 'un',0); % convert to cell array of numeric array
% [m343, tf] = padcat(Combs2{:}); % concatenate, pad rows with NaNs
% m343(~tf) = 0 ;% replace NaNs by zeros

