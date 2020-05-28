%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  BPSO, VPSO  and PSEQEA source codes version 1.0                  %
%  Developed in MATLAB R2015b                                       %
%                                                                   %
%  Author and programmer: Carlos Kuhn                               %
%   email: Carlos.Kuhn@dst.defence.gov.au                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
% close all
clc
% 
% addpath(genpath('/Users/carlosclaitonnoschangkuhn/Google Drive/MatLabCodes/OptiizationAlgorithm/CK_Implementations/BSPO_BudgetScenario/'))
% savepath
    

DepLevel = 1;          % Dependence Levels use DepLevel = 5 -> Independent
                        %                      DepLevel = 4 -> Tangential
                        %                      DepLevel = 3 -> Associated
                        %                      DepLevel = 2 -> Dependente
                        %                      DepLevel = 1 -> Mandatory

% BudgetArray = 10;                     % Budget array
% BudgetArray = 196.46;                     % Budget (Greiner et al (2003)
BudgetArray = 190.28;                     % Budget (Greiner et al (2003)
numberRun =50;
gBest = cell(3,length(BudgetArray),numberRun);   % store the best initiatives
gBestScore = ones(3,length(BudgetArray),numberRun)*NaN;  % keep the best score

% optimunBeta = ones(length(BudgetArray),length(DepLevel))*NaN;

Max_iteration=1000; % Maximum number of iterations

ConvergenceCurves=zeros(3,Max_iteration);

for DPI = DepLevel     
    noP = 50; % Number of particles on the Swarm 
    ind = 0;                    
    for Budget = BudgetArray
        ind = ind+1;
        model=CreateModel_QKP(Budget,DPI); % here is where the model is stablished 
        noV = length(model.Cost); % number of decision variables, match with the size of cost vector
        CostFunction=@(x) MyCost_QKP(x,model); % Modify or replace Mycost.m according to your cost funciton

        %BPSO with s-shaped family of transfer functions
        tic
        for k=1:numberRun
            [gBest{1,ind,k}, gBestScore(1,ind,k) ,ConvergenceCurves(1,:)]=BPSO(noP,Max_iteration,4,CostFunction,noV);
%             disp(['Budget: ' num2str(Budget) ' Initiatives: ' cell2mat(model.Initiative(gBest{1,ind}~=0)') ' Score= ' num2str(1/gBestScore(1,ind))])
        end
         toc

    %     %BPSO with v-shaped family of transfer functions
%         for k=1:50
%             [gBest{2,ind}, gBestScore(2,ind,k) ,ConvergenceCurves(2,:)]=BPSO(noP,Max_iteration,8,CostFunction,noV);
% %         disp(['Budget: ' num2str(Budget) ' Initiatives: ' cell2mat(model.Initiative(gBest{1,ind}~=0)') ' Score= ' num2str(1/gBestScore(1,ind))])
%         end
%          toc
    
%       PSEQEA
%         tic
%         for k=1:numberRun
%             [gBest{3,ind,k}, gBestScore(3,ind,k) ,ConvergenceCurves(3,:)]=PSEQEA(noP,Max_iteration,CostFunction,noV);
% %        disp(['Budget: ' num2str(Budget) ' Initiatives: ' cell2mat(model.Initiative(gBest{2,ind}~=0)') ' Score= ' num2str(1/gBestScore(2,ind))])
%         end
%         toc
    end
    
end


%comparing the both BSPO methods by convergency curves
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
legend({'BSPO','VPSO','PSEQEA'},'FontSize',18,'Location','southeast');
hold off
title(['Convergence Curve, Budget: ' num2str(BudgetArray)],'FontSize', 18)

% save ('BetaModel.mat','BudgetArray','optimunBeta')

%% Getting the stats 
CostReq = ones(size(gBest,3),1)*NaN;
for i=1:size(gBest,3)
    CostReq(i) = sum(model.Cost(gBest{1,1,1}~=0));
end
    

testeBSPO=1./reshape(gBestScore(1,1,:),[1,numberRun]);
% testeBSPO=1-reshape(gBestScore(1,1,:),[1,50]);
% testeVSPO=1./reshape(gBestScore(2,1,:),[1,50]);
% testePSEQEA=1./reshape(gBestScore(3,1,:),[1,50]);
Std(1) = std(testeBSPO);
% Std(2) = std(testeVSPO);
% Std(3) = std(testePSEQEA);
Mean = mean(1./gBestScore,3);
% Mean = mean(testeBSPO);
% 
optAnalitics = round(mode(testeBSPO),4);
% 
HitsBSPO = sum(round(testeBSPO,4)==optAnalitics);
% HitsVSPO = sum(round(testeVSPO,2)==optAnalitics);
% HitsPSEQEA = sum(round(testePSEQEA,2)==optAnalitics);
% 
disp(['BSPO ; Budget Mean = ' num2str(mean(CostReq))  ', Portfolio Mean = ' num2str(Mean(1)) ', sdt = ' num2str(Std(1)) ', Hits: ' num2str(HitsBSPO) ])
disp(['Mode = ' num2str(optAnalitics)  ])

% disp(['VSPO ; Mean = ' num2str(Mean(2)) ' sdt = ' num2str(Std(2)) ' Hits: ' num2str(HitsVSPO) ])
% disp(['PSEQEA ; Mean = ' num2str(Mean(3)) ' sdt = ' num2str(Std(3)) ' Hits: ' num2str(HitsPSEQEA) ])
% 
% wang no dependence clain to be
gbestWang = ones(15,1);
gbestWang([5,9])=0;
sum(model.Cost(gbestWang~=0));
