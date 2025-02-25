%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  BPSO, VPSO, PSEQEA and QL source codes version 2.0               %
%  Developed in MATLAB R2015b                                       %
%                                                                   %
%  Author and programmer: Carlos Kuhn                               %
%   email: Carlos.Kuhn@dst.defence.gov.au                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
% close all
clc

%ADDIND THE PATH WHERE IT HAS THE FUCTIONS CALLED
addpath(genpath('/Users/carlosclaitonnoschangkuhn/Google Drive/MatLabCodes/OptiizationAlgorithm/CK_Implementations/BSPO_BudgetScenario/'))
addpath(genpath('/Users/carlosclaitonnoschangkuhn/Google Drive/MatLabCodes/OptiizationAlgorithm/CK_Implementations/Q_learning/'))
savepath

Max_iteration=5000;  % Maximum number of iterations
NumbRuns = 10;     % Nubmer of runs

%%%%///////////////DEPENDENCE LEVEL//////////////////////////////////////
DepLevel = 1;          % Dependence Levels use DepLevel = 5 -> Independent
                        %                      DepLevel = 4 -> Tangential
                        %                      DepLevel = 3 -> Associated
                        %                      DepLevel = 2 -> Dependente
                        %                      DepLevel = 1 -> Mandatory

BudgetArray =6;%5;                     % Budget array
OptAnalytical = 7.1;%6.8;                     % analytical Optimum score

NumbInitiative = 15;%,50,100];%,150];%300,400,500]';

gBest = cell(4,length(BudgetArray));   % store the best initiatives
gBestScore = ones(4,length(BudgetArray),NumbRuns)*NaN;  % keep the best score

optimunBeta = ones(length(BudgetArray),length(DepLevel))*NaN;

ConvergenceCurves=ones(4,Max_iteration,NumbRuns)*NaN;
IteNumber = ones(4,NumbRuns,length(NumbInitiative))*NaN;
UsedBudget =  ones(4,NumbRuns,length(NumbInitiative))*NaN;

timeEvol = ones(4,length(BudgetArray))*NaN;

for DPI = DepLevel     
    noP = 30; % Number of particles on the Swarm or number of q-bit or number of agents
    ind = 0;                    
    Budget = BudgetArray;
    for NumIni = NumbInitiative %for Budget = BudgetArray
        ind = ind+1;
        %Uses the data from Order paper (Priority List)
%         model=CreateModel_CK(Budget,DPI); % here is where the model is stablished 

%       %Uses the data Generated randomly 
        model=CreateModelBidData_CK(Budget,DPI,NumIni); % here is where the model is stablished 
        
        noV = length(model.Cost); % number of decision variables, match with the size of cost vector
        CostFunction=@(x) MyCost_CK(x,model); % Modify or replace Mycost.m according to your cost funciton

        % BPSO with s-shaped family of transfer functions
        tic
        for k=1:NumbRuns
            [gBest{1,ind}, gBestScore(1,ind,k) ,ConvergenceCurves(1,:,k),IteNumber(1,k,ind)]=BPSO(noP,Max_iteration,2,CostFunction,noV);
            UsedBudget(1,k,ind) = sum(gBest{1,ind}.*model.Cost);
        % disp(['Budget: ' num2str(Budget) ' Initiatives: ' cell2mat(model.Initiative(gBest{1,ind}~=0)') ' Score= ' num2str(1/gBestScore(1,ind))])
        end
        timeEvol(1,ind)=toc;

        % BPSO with v-shaped family of transfer functions
        for k=1:NumbRuns
            [gBest{2,ind}, gBestScore(2,ind,k) ,ConvergenceCurves(2,:,k),IteNumber(2,k,ind)]=BPSO(noP,Max_iteration,8,CostFunction,noV);
            UsedBudget(2,k,ind) = sum(gBest{2,ind}.*model.Cost);
        % disp(['Budget: ' num2str(Budget) ' Initiatives: ' cell2mat(model.Initiative(gBest{1,ind}~=0)') ' Score= ' num2str(1/gBestScore(1,ind))])
        end
        timeEvol(2,ind)=toc;
%     
        % PSEQEA
        tic
        for k=1:NumbRuns
            [gBest{3,ind}, gBestScore(3,ind,k) ,ConvergenceCurves(3,:,k),IteNumber(3,k,ind)]=PSEQEA(noP,Max_iteration,CostFunction,noV);
            UsedBudget(3,k,ind) = sum(gBest{3,ind}.*model.Cost);
        % disp(['Budget: ' num2str(Budget) ' Initiatives: ' cell2mat(model.Initiative(gBest{2,ind}~=0)') ' Score= ' num2str(1/gBestScore(2,ind))])
        end
        timeEvol(3,ind)=toc;
%         
        % Q-learning
        tic
        for k=1:NumbRuns
            [gBest{4,ind}, gBestScore(4,ind,k) ,ConvergenceCurves(4,:,k),IteNumber(4,k,ind)]=BS_Qlearning(noP,Max_iteration,CostFunction,noV);
            UsedBudget(4,k,ind) = sum(gBest{4,ind}.*model.Cost);
        % disp(['Budget: ' num2str(Budget) ' Initiatives: ' cell2mat(model.Initiative(gBest{2,ind}~=0)') ' Score= ' num2str(1/gBestScore(2,ind))])
        end
        timeEvol(4,ind)=toc;
%         
        
    end
    
end

%% comparing the both BSPO methods by convergency curves
ConvergenceCurves(1,:,:) = 1./ConvergenceCurves(1,:,:); % changing the output format
ConvergenceCurves(2,:,:) = 1./ConvergenceCurves(2,:,:);
ConvergenceCurves(3,:,:) = 1./ConvergenceCurves(3,:,:);
%%Taking the mean value on the number of runs
MeanConv = mean(ConvergenceCurves,3);
StdConv = std(ConvergenceCurves,[],3);

figure(1)
col=hsv(size(ConvergenceCurves,1));
count=0;
marker = ['o','s','d','*','v'];
for i= 1:size(ConvergenceCurves,1)
%     ydata = ConvergenceCurves(i,:);
    ydata = MeanConv(i,:);
    xdata = 1:size(ConvergenceCurves,2);
    ydataErr = StdConv(i,:);
%     p1=plot(xdata,ydata,'Color',col(i,:));
    p1 = errorbar(xdata,ydata,ydataErr,ydataErr,'Color',col(i,:));
    p1.LineStyle  = '-'; p1.LineWidth = 0.5; p1.Marker = marker(i); 
    p1.MarkerSize = 8; p1.MarkerFaceColor = col(i,:);
    hold on
end
    
set(gca,'FontSize',14); %set the axis font size
%set(gca,'xtick',0:0.25:1.5)
box on;
grid on;
xlabel('Iteration', 'FontSize', 18)
ylabel('Optimal Score', 'FontSize', 18)
legend({'BPSO','VPSO','PSEQEA','QL'},'FontSize',18,'Location','southeast');
hold off
title(['Convergence Curve, Budget: ' num2str(BudgetArray)],'FontSize', 18)


%% Getting the stats 
Methode = {'BPSO';'VPSO';'PSEQEA';'QL'};
MeanScore = ones(4,ind)*NaN;
StdScore = ones(4,ind)*NaN;
ModeScore = ones(4,ind)*NaN;
MeanBudget = ones(4,ind)*NaN;
MeanIterationNumber = ones(4,ind)*NaN;
Hits = ones(4,ind)*NaN;
for i=1:ind
    testeBSPO = 1./reshape(gBestScore(1,i,:),[1,NumbRuns]);
    testeVSPO = 1./reshape(gBestScore(2,i,:),[1,NumbRuns]);
    testePSEQEA = 1./reshape(gBestScore(3,i,:),[1,NumbRuns]);
    testeQL  = reshape(gBestScore(4,i,:),[1,NumbRuns]);

    %Combine all data
    CompMatrix = [testeBSPO',testeVSPO',testePSEQEA',testeQL'];
    MeanScore(:,i) = mean(CompMatrix);
    StdScore(:,i) = std(CompMatrix);
    MeanBudget(:,i) = mean(UsedBudget(:,:,i),2);
    MeanIterationNumber(:,i) = mean(IteNumber(:,:,i),2);

    ModeScore(1,i) = round(mode(testeBSPO),2);
    ModeScore(2,i) = round(mode(testeVSPO),2);
    ModeScore(3,i) = round(mode(testePSEQEA),2);
    ModeScore(4,i) = round(mode(testeQL),2);

%     Hits(1,i) = sum(round(testeBSPO,2)==ModeData(3,i));
%     Hits(2,i) = sum(round(testeVSPO,2)==ModeData(3,i));
%     Hits(3,i) = sum(round(testePSEQEA,2)==ModeData(3,i));
%     Hits(4,i) = sum(round(testeQL,2)==ModeData(3,i));
    
    Hits(1,i) = sum(round(testeBSPO,2)==OptAnalytical);
    Hits(2,i) = sum(round(testeVSPO,2)==OptAnalytical);
    Hits(3,i) = sum(round(testePSEQEA,2)==OptAnalytical);
    Hits(4,i) = sum(round(testeQL,2)==OptAnalytical);

%     disp([ 'Budge: ' num2str(BudgetArray(i)) ',  Initiatives: ' cell2mat(model.Initiative(gBest{3,i}~=0)') ])
    disp([ 'Budget: ' num2str(Budget) ', Number of Initiatives: ' num2str(NumbInitiative(i)), ' Initiative BSPO: ' cell2mat(model.Initiative(gBest{3,i}~=0)') ])
    DataStat = table(Methode, MeanBudget(:,i), MeanScore(:,i), StdScore(:,i), Hits(:,i), timeEvol(:,i), MeanIterationNumber(:,i));
    disp(DataStat)
%     disp(['BPSO: Mean Budget = ' num2str(MeanBudget(1,i)) ' Mean Score = ' num2str(MeanScore(1,i)) ' std Score = ' num2str(StdScore(1,i)) ' Hits: ' num2str(Hits(1,i))  '  time (sec): ' num2str(timeEvol(1,i)) ' iterations: ' num2str(MeanIterationNumber(1,i)) ])
%     disp(['VPSO: Mean Budget = ' num2str(MeanBudget(2,i)) ' Mean Score = ' num2str(MeanScore(2,i)) ' std Score= ' num2str(StdScore(2,i)) ' Hits: ' num2str(Hits(2,i)) ' time (sec): ' num2str(timeEvol(2,i)) ' iterations: ' num2str(MeanIterationNumber(2,i)) ])
%     disp(['PSEQEA: Mean Budget = ' num2str(MeanBudget(3,i)) ' Mean Score = ' num2str(MeanScore(3,i)) ' std Score = ' num2str(StdScore(3,i)) ' Hits: ' num2str(Hits(3,i)) ' time (sec): ' num2str(timeEvol(3,i)) ' iterations: ' num2str(MeanIterationNumber(3,i)) ])
%     disp(['QL: Mean Budget = ' num2str(MeanBudget(4,i)) ' Mean = ' num2str(MeanScore(4,i)) ' std = ' num2str(StdScore(4,i)) ' Hits: ' num2str(Hits(4,i)) ' time (sec): ' num2str(timeEvol(4,i)) ' iterations: ' num2str(MeanIterationNumber(4,i))])
end

%%
% save ('BetaModelVaryInitiatives500iterations50_Att2.mat','BudgetArray','ModeScore','gBest','gBestScore','CompMatrix',...
%       'Hits','MeanScore','StdScore','timeEvol', 'MeanBudget', 'IteNumber', 'UsedBudget');

%% Plots

figure(2)
bar(NumbInitiative,Hits')
box on;
grid on;
xlabel('Number of Initiatives', 'FontSize', 18)
ylabel('Number of Hits', 'FontSize', 18)
ylim([0,55])
legend({'BPSO','VPSO','PSEQEA','QL'},'FontSize',18,'Location','northeast');

figure(3)
bar(NumbInitiative,MeanScore')
% bar(model_series, 'grouped');
hold on
% Find the number of groups and the number of bars in each group
ngroups = size(MeanScore', 1);
nbars = size(MeanScore', 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 2.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
% for i = 1:nbars
%     % Calculate center of each bar
%     x = [15,50,100,150,300,400,500] - 3.6*groupwidth + (12*i-1) * groupwidth / (2*nbars);
%     errorbar(x, MeanData(i,:), StdData(i,:), 'k', 'linestyle', 'none');
% end
hold off
box on;
grid on;
xlabel('Number of Initiatives', 'FontSize', 18)
ylabel('Mean Best Score', 'FontSize', 18)
ylim([0,10])
legend({'BPSO','VPSO','PSEQEA','QL'},'FontSize',18,'Location','northeast');


