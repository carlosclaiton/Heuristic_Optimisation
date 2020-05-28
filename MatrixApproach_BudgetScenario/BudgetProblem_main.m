%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS CODE SOLVE THE HYPOTHETICAL PROBLEM PROPOST BY N. ORDER (2007 AND
% 2009 PRIORITY LIST) 
%
% Using teh function PLMatrix that solve the problem when there is no
% dependence between the Initiatives
% 
% History
% Created by Carlos Kuhn
% Created at 12/03/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%FEED THE INITIATIVE AND COST
clear all
tic
Initiative = [{'O'},{'M'},{'N'},{'D'},{'J'},{'I'},{'B'},{'L'},{'E'},{'A'},...
    {'H'},{'G'},{'C'},{'K'},{'F'}]';
Cost = [1,1,3,4,7,6,8,11,13,12,14,17,14,20,16];
[Cost,IndexCost] = sort(Cost);
Initiative =Initiative(IndexCost);
budget = 1:29;

%%CALL THE FUNCTION THAT FOUND THE OPTMUM WITH ZERO DEPENDENCE
OptSubFound = cell(length(budget),1);
labels = cell(length(budget),1);
optimun = ones(length(budget),1)*NaN;
for i = budget
    [OptSubFound{i},optimun(i)] = PLMatrix(i);
    if size(OptSubFound{i},1)>2
        labels1 = { Initiative(OptSubFound{i}(1,OptSubFound{i}(1,:)~=0)) , ...
            {' or '}, Initiative(OptSubFound{i}(2,OptSubFound{i}(2,:)~=0)), ...
             {' or '}, Initiative(OptSubFound{i}(3,OptSubFound{i}(3,:)~=0)) };
        labels{i} = cat(1, labels1{:});
    elseif size(OptSubFound{i},1)>1
        labels1 = { Initiative(OptSubFound{i}(1,OptSubFound{i}(1,:)~=0)) , ...
            {' or '}, Initiative(OptSubFound{i}(2,OptSubFound{i}(2,:)~=0)) };
        labels{i} = cat(1, labels1{:});
    else
        labels{i} = Initiative(OptSubFound{i});
    end    
end

% figure(1)
% plot(budget,optimun,'-o')
% text(budget,optimun,labels,'VerticalAlignment','top','HorizontalAlignment','right')
% grid on
% grid(gca,'minor')
% xlabel('Cost Bound', 'FontSize', 16)
% ylabel('Optimal Score', 'FontSize', 16)
% legend({'Dependence Level = Zero'}, 'FontSize', 16)


%%SOLVING THE SAME BUDGET PROBLEM WITH DEPENDENCE

OptSubFoundD = cell(length(budget),1);
optimunD = ones(length(budget),1)*NaN;
labelsD = cell(length(budget),1);
for i = budget
    [OptSubFoundD{i},optimunD(i)] = PLMatrix_Dependence(i);
    if size(OptSubFoundD{i},1)>1
        labels1 = { Initiative(OptSubFoundD{i}(1,OptSubFoundD{i}(1,:)~=0)) , ...
            {' or '}, Initiative(OptSubFoundD{i}(2,OptSubFoundD{i}(2,:)~=0)) };
        labelsD{i} = cat(1, labels1{:});
    else
        labelsD{i} = Initiative(OptSubFoundD{i});
    end    
end
%%PRINTING THE OPTIMUM FOR GIVEN BUDGET

b=29; %BUDGET TO BE TESTED
for i=1:size(OptSubFoundD{b},1)
    disp(Initiative(OptSubFoundD{b}(i,OptSubFoundD{b}(i,:)~=0))')
end
%%

% figure(3)
% plot(budget,optimunD,'-o')
% text(budget,optimunD,labelsD,'VerticalAlignment','top','HorizontalAlignment','right')
% grid on
% grid(gca,'minor')
% xlabel('Cost Bound', 'FontSize', 16)
% ylabel('Optimal Score', 'FontSize', 16)
% legend({'Dependence Level = 4'}, 'FontSize', 16)
toc

%% PLOTTING BOTH RESULTS ON THE SAME GRAFIC

load('BetaModel.mat')
figure(4)
plot(budget,optimun,'-or')
hold on
%text(budget,optimun,labels,'VerticalAlignment','bottom','HorizontalAlignment','right')
plot(budget,optimunD,'-db')
%text(budget,optimunD,labelsD,'VerticalAlignment','top','HorizontalAlignment','right')
plot(BudgetArray,optimunBeta(:,5),'-sg')
plot(BudgetArray,optimunBeta(:,1),'-*m')
hold off
grid on
grid(gca,'minor')
xlabel('Cost Bound', 'FontSize', 16)
ylabel('Optimal Score', 'FontSize', 16)
legend({'Dep Level = 0','Dep Level = 4','Dep Level = 0, Beta Model','Dep Level = 4 Beta Model'}, 'FontSize', 16)
ylim([4,10])
