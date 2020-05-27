%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code uses the Bellmann Equation to determine the action of the
% agent. We have 10 agents that interact to eacho other to decide what is
% the best solucion subset to maximize the fiteness. 
%
% Developer: Carlos C. N. Kuhn
% email: carlosclaitonkuhn@gmail.com
% Created at 23/04/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
% close all
% clc


DepLevel = 1;          % Dependence Levels use DepLevel = 5 -> Independent
                        %                      DepLevel = 4 -> Tangential
                        %                      DepLevel = 3 -> Associated
                        %                      DepLevel = 2 -> Dependente
                        %                      DepLevel = 1 -> Mandatory

BudgetArray = 1:30;                     % Budget array
InitiativeNumber = 15;
gBest = cell(length(BudgetArray));   % store the best initiatives
gBestScore = ones(length(BudgetArray),InitiativeNumber)*NaN;  % keep the best score
IterNumb = ones(length(BudgetArray),InitiativeNumber)*NaN;
SubInit = cell(length(BudgetArray),InitiativeNumber);
optimunBeta = ones(length(BudgetArray),length(DepLevel))*NaN;

Max_iteration=2000; % Maximum number of iterations

ConvergenceCurves=zeros(1,Max_iteration);
tic
for DPI = DepLevel     
    noP = 30; % Number of particles on the Swarm / or number of q-bit
    ind = 0;                    
    for Budget = BudgetArray
        ind = ind+1;
        for  indexIni=1:InitiativeNumber
            model=CreateMultiModel_CK(Budget,DPI,indexIni); % here is where the model is stablished 
            noV = length(model.Cost); % number of decision variables, match with the size of cost vector
            CostFunction=@(x) MyCost_CK(x,model); % Modify or replace Mycost.m according to your cost funciton

            %BPSO with s-shaped family of transfer functions
            [gBest{ind,indexIni}, gBestScore(ind,indexIni) ,ConvergenceCurves(1,:),IterNumb(ind,indexIni)]=BS_Qlearning(noP,Max_iteration,CostFunction,noV);
%             disp(['Budget: ' num2str(Budget) ', Initiatives: ' cell2mat(model.Initiative(gBest{ind,indexIni}~=0)') ', Score= ' num2str(1/gBestScore(ind,indexIni))])
            SubInit{ind,indexIni} = cell2mat(model.Initiative(gBest{ind,indexIni}~=0)');
        end
    end
end
toc
%% Filter to get the unique soluction per budget
optimum = ones(1,length(BudgetArray))*NaN;
IndOpt = ones(length(BudgetArray),InitiativeNumber)*NaN;
SubIntOpt = cell(length(BudgetArray),1);
scoreMat = gBestScore;
for i=1:length(BudgetArray)
    optimum(i) = max(scoreMat(i,:));
    IndOpt(i,:) = scoreMat(i,:)==optimum(i);
    teste={SubInit{i,IndOpt(i,:)~=0}};
    SubIntOpt{i} = unique(teste);
    
    disp(['Budget: ' num2str(BudgetArray(i))...
        ', Score= ' num2str(optimum(i)) ', Initiatives:' SubIntOpt{i} ])

end
%%
figure(1)
ydata =optimum;
xdata =BudgetArray;
p1=plot(xdata,ydata,'Color','b');
p1.LineStyle  = '-'; p1.LineWidth = 0.5; p1.Marker = 'o'; 
p1.MarkerSize = 8; p1.MarkerFaceColor = 'b';    
set(gca,'FontSize',16); %set the axis font size
box on;
grid on;
xlabel('Cost Bound', 'FontSize', 18)
ylabel('Optimal Score', 'FontSize', 18)

%% 
% figure(1)
% col=hsv(size(optimunBeta,2));
% count=0;
% marker = ['o','s','d','*','v'];
% for i= DepLevel
%     count = count+1;
%     ydata =optimunBeta(:,i);
%     xdata =BudgetArray;
%     p1=plot(xdata,ydata,'Color',col(count,:));
%     p1.LineStyle  = '-'; p1.LineWidth = 0.5; p1.Marker = marker(count); 
%     p1.MarkerSize = 8; p1.MarkerFaceColor = col(count,:);
%     leg1{count} = ['Dep Level  = ' num2str(DepLevel(count) - 1) ]; 
%     hold on
% end
%     
% set(gca,'FontSize',14); %set the axis font size
% %set(gca,'xtick',0:0.25:1.5)
% box on;
% grid on;
% xlabel('Cost Bound', 'FontSize', 18)
% ylabel('Optimal Score', 'FontSize', 18)
% legend(leg1,'FontSize',18,'Location','northeast');
% hold off

% COMPARE BOTH BSPO IMPLEMENTATION

% figure(2)
% col=hsv(size(ConvergenceCurves,1));
% count=0;
% marker = ['o','s','d','*','v'];
% for i= 1:size(ConvergenceCurves,1)
%     ydata = 1./ConvergenceCurves(i,:);
%     xdata = 1:size(ConvergenceCurves,2);
%     p1=plot(xdata,ydata,'Color',col(i,:));
%     p1.LineStyle  = '-'; p1.LineWidth = 0.5; p1.Marker = marker(i); 
%     p1.MarkerSize = 8; p1.MarkerFaceColor = col(i,:);
%     hold on
% end
%     
% set(gca,'FontSize',14); %set the axis font size
% %set(gca,'xtick',0:0.25:1.5)
% box on;
% grid on;
% xlabel('Iteration', 'FontSize', 18)
% ylabel('Optimal Score', 'FontSize', 18)
% legend({'BSPO','PSEQEA'},'FontSize',18,'Location','southeast');
% hold off
% title(['Convergence Curve, Budget: ' num2str(BudgetArray)],'FontSize', 18)

% save ('BetaModel.mat','BudgetArray','optimunBeta')

%% Getting the stats

% testeBSPO=1./reshape(gBestScore(1,1,:),[1,50]);
% testePSEQEA=1./reshape(gBestScore(2,1,:),[1,50]);
% Std(1) = std(testeBSPO);
% Std(2) = std(testePSEQEA);
% Mean = mean(1./gBestScore,3);
% 
% HitsBSPO = sum(round(testeBSPO,2)==8.10);
% HitsPSEQEA = sum(round(testePSEQEA,2)==8.10);
% 
% disp(['BSPO ; Mean = ' num2str(Mean(1)) ' sdt = ' num2str(Std(1)) ' Hits: ' num2str(HitsBSPO) ])
% disp(['PSEQEA ; Mean = ' num2str(Mean(2)) ' sdt = ' num2str(Std(2)) ' Hits: ' num2str(HitsPSEQEA) ])
% 


