%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  BPSO, VPSO  and PSEQEA source codes version 1.0                  %
%  Developed in MATLAB R2015b                                       %
%                                                                   %
%  Author and programmer: Carlos Kuhn                               %
%   email: Carlos.Kuhn@dst.defence.gov.au                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all 
% close all
clc

%ADDIND THE PATH WHERE IT HAS THE FUCTIONS CALLED
addpath(genpath('/Users/carlosclaitonnoschangkuhn/Google Drive/MatLabCodes/OptiizationAlgorithm/CK_Implementations/BSPO_BudgetScenario/'))
addpath(genpath('/Users/carlosclaitonnoschangkuhn/Google Drive/MatLabCodes/OptiizationAlgorithm/CK_Implementations/Q_learning/'))
savepath

Max_iteration=300;  % Maximum number of iterations
NumbRuns = 10;     % Nubmer of runs

%%%%///////////////DEPENDENCE LEVEL//////////////////////////////////////
DepLevel = 1;          % Dependence Levels use DepLevel = 5 -> Independent
                        %                      DepLevel = 4 -> Tangential
                        %                      DepLevel = 3 -> Associated
                        %                      DepLevel = 2 -> Dependente
                        %                      DepLevel = 1 -> Mandatory

BudgetArray = 5;%28;%10;%5;                     % Budget array
optAnalitics = 6.8;%9.1;%8.1;%6.8;                     % analytical Optimum score


gBest = cell(4,length(BudgetArray));   % store the best initiatives
gBestScore = ones(4,length(BudgetArray),NumbRuns)*NaN;  % keep the best score

optimunBeta = ones(length(BudgetArray),length(DepLevel))*NaN;

ConvergenceCurves=zeros(4,Max_iteration);

% dimension 2 we change the neumber of initiatives, dimension one is the 4
% model 
timeEvol = ones(4,10)*NaN;
dataIndex = 0;
%%
indexIni = [4 5 6 7 8 9 10 11 12 13 14]; %initiatives to be removed from the model
dataIndex = 1+dataIndex;

for DPI = DepLevel     
    noP = 30; % Number of particles on the Swarm or number of q-bit or number of agents
    ind = 0;                    
    for Budget = BudgetArray
        ind = ind+1;
%         model=CreateModel_CK(Budget,DPI); % here is where the model is stablished 
     %   IndexIni contem the index of Initiatives want ro remove from
     %   the origional model
        model=CreateMultiModel_CK(Budget,DPI,indexIni); % here is where the model is stablished 

        noV = length(model.Cost); % number of decision variables, match with the size of cost vector
        CostFunction=@(x) MyCost_CK(x,model); % Modify or replace Mycost.m according to your cost funciton

        % BPSO with s-shaped family of transfer functions
        tic
        for k=1:NumbRuns
            [gBest{1,ind}, gBestScore(1,ind,k) ,ConvergenceCurves(1,:)]=BPSO(noP,Max_iteration,4,CostFunction,noV);
        % disp(['Budget: ' num2str(Budget) ' Initiatives: ' cell2mat(model.Initiative(gBest{1,ind}~=0)') ' Score= ' num2str(1/gBestScore(1,ind))])
        end
        timeEvol(1,dataIndex)=toc;

        % BPSO with v-shaped family of transfer functions
        for k=1:NumbRuns
            [gBest{2,ind}, gBestScore(2,ind,k) ,ConvergenceCurves(2,:)]=BPSO(noP,Max_iteration,8,CostFunction,noV);
        % disp(['Budget: ' num2str(Budget) ' Initiatives: ' cell2mat(model.Initiative(gBest{1,ind}~=0)') ' Score= ' num2str(1/gBestScore(1,ind))])
        end
        timeEvol(2,dataIndex)=toc;
    
        % PSEQEA
        tic
        for k=1:NumbRuns
            [gBest{3,ind}, gBestScore(3,ind,k) ,ConvergenceCurves(3,:)]=PSEQEA(noP,Max_iteration,CostFunction,noV);
        % disp(['Budget: ' num2str(Budget) ' Initiatives: ' cell2mat(model.Initiative(gBest{2,ind}~=0)') ' Score= ' num2str(1/gBestScore(2,ind))])
        end
        timeEvol(3,dataIndex)=toc;
        
        % Q-learning
        tic
        for k=1:NumbRuns
            [gBest{4,ind}, gBestScore(4,ind,k) ,ConvergenceCurves(4,:)]=BS_Qlearning(noP,Max_iteration,CostFunction,noV);
        % disp(['Budget: ' num2str(Budget) ' Initiatives: ' cell2mat(model.Initiative(gBest{2,ind}~=0)') ' Score= ' num2str(1/gBestScore(2,ind))])
        end
        timeEvol(4,dataIndex)=toc;
        
        
    end
    
end

%%comparing the both BSPO methods by convergency curves
ConvergenceCurves(4,:) = 1./ConvergenceCurves(4,:); % correting the output
figure(1)
col=hsv(size(ConvergenceCurves,1));
count=0;
marker = ['o','s','d','*','v'];
for i= 1:size(ConvergenceCurves,1)
    ydata = 1./ConvergenceCurves(i,:);
    xdata = 1:size(ConvergenceCurves,2);
    p1=plot(xdata,ydata,'Color',col(i,:));
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

% save ('BetaModel.mat','BudgetArray','optimunBeta')

%%Getting the stats 

testeBSPO = 1./reshape(gBestScore(1,1,:),[1,NumbRuns]);
testeVSPO = 1./reshape(gBestScore(2,1,:),[1,NumbRuns]);
testePSEQEA = 1./reshape(gBestScore(3,1,:),[1,NumbRuns]);
testeQL  = reshape(gBestScore(4,1,:),[1,NumbRuns]);

%Combine all data
CompMatrix = [testeBSPO',testeVSPO',testePSEQEA',testeQL'];
MeanData = mean(CompMatrix);
StdData = std(CompMatrix);

HitsBSPO = sum(round(testeBSPO,2)==optAnalitics);
HitsVSPO = sum(round(testeVSPO,2)==optAnalitics);
HitsPSEQEA = sum(round(testePSEQEA,2)==optAnalitics);
HitsQL = sum(round(testeQL,2)==optAnalitics);

disp(['BPSO: Mean = ' num2str(MeanData(1)) ' std = ' num2str(StdData(1)) ' Hits: ' num2str(HitsBSPO)  ' time (sec): ' num2str(timeEvol(1,dataIndex))])
disp(['VPSO: Mean = ' num2str(MeanData(2)) ' std = ' num2str(StdData(2)) ' Hits: ' num2str(HitsVSPO) ' time (sec): ' num2str(timeEvol(2,dataIndex))])
disp(['PSEQEA: Mean = ' num2str(MeanData(3)) ' std = ' num2str(StdData(3)) ' Hits: ' num2str(HitsPSEQEA) ' time (sec): ' num2str(timeEvol(3,dataIndex))])
disp(['QL: Mean = ' num2str(MeanData(4)) ' std = ' num2str(StdData(4)) ' Hits: ' num2str(HitsQL) ' time (sec): ' num2str(timeEvol(4,dataIndex))])

% disp('Method  Mean   Std  hits time (sec) ')
% disp(['BPSO ' num2str(MeanData(1)) '  ' num2str(StdData(1)) '  ' num2str(HitsBSPO)  '  ' num2str(timeEvol(1))])
% disp(['VPSO ' num2str(MeanData(2)) '  ' num2str(StdData(2)) ' ' num2str(HitsVSPO) '  ' num2str(timeEvol(2))])
% disp(['PSEQEA ' num2str(MeanData(3)) '  ' num2str(StdData(3)) '  ' num2str(HitsPSEQEA) '  ' num2str(timeEvol(3))])
% disp(['QL ' num2str(MeanData(4)) '  ' num2str(StdData(4)) '  ' num2str(HitsQL) '  ' num2str(timeEvol(4))])
% 

%% Ploting the Evolution time


figure(2)
col=hsv(size(timeEvol,1));
count=0;
marker = ['o','s','d','*','v'];
for i= 1:size(ConvergenceCurves,1)
    ydata = timeEvol(i,:);
    xdata = [15,14,12,10,8,6,4,NaN,NaN,NaN];
    p1=plot(xdata,ydata,'Color',col(i,:));
    p1.LineStyle  = '-'; p1.LineWidth = 0.5; p1.Marker = marker(i); 
    p1.MarkerSize = 8; p1.MarkerFaceColor = col(i,:);
    hold on
end
    
set(gca,'FontSize',14); %set the axis font size
%set(gca,'xtick',0:0.25:1.5)
box on;
grid on;
xlabel('Number of Initiatives', 'FontSize', 18)
ylabel('Time of 10 runs [sec]', 'FontSize', 18)
legend({'BPSO','VPSO','PSEQEA','QL'},'FontSize',18,'Location','southeast');
hold off
title(['300 Iteration, Budget: ' num2str(BudgetArray)],'FontSize', 18)



