%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  BPSO and VPSO source codes version 1.0                           %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper: S. Mirjalili and A. Lewis, "S-shaped versus         %
%               V-shaped transfer functions for binary Particle     %
%               Swarm Optimization," Swarm and Evolutionary         %
%               Computation, vol. 9, pp. 1-14, 2013.                %
%                                                                   %
%   20/03/20 - Adapeted to solve the Budget Scenario problem        %
%   Carlos C. N. Kuhn                                               % 
%   email: Carlos.Kuhn@dst.defence.gov.au                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
% close all
% clc

tic
DepLevel = 1;           %Dependence Levels use DepLevel = 5 -> Independent
                        %                      DepLevel = 4 -> Tangential
                        %                      DepLevel = 3 -> Associated
                        %                      DepLevel = 2 -> Dependente
                        %                      DepLevel = 1 -> Mandatory

BudgetArray = 10;        % entre the list of the budget 


Max_iteration=300;      % Maximum number of iterations
Max_SubIteration = 30;  % Maximum of sub-iterations to find the multiples minimum

gBest1 = cell(Max_SubIteration,length(BudgetArray));
gBestScore1 = ones(Max_SubIteration,length(BudgetArray))*NaN;
optimunBeta = ones(length(BudgetArray),length(DepLevel))*NaN;
InitOpt = cell(1,length(BudgetArray));
count = 0;
for DPI = DepLevel;     
    count = count+1;
    count1 = 0;
    for Budget = BudgetArray
        count1 = count1+1;
        model=CreateModel_CK(Budget,DPI); % here is where the model is stablished 
        CostFunction=@(x) MyCost_CK(x,model); % Modify or replace Mycost.m according to your cost funciton

        noP=30;                 % Number of particles on the Swarm
        noV = length(model.Cost); % number of decision variables, match with the size of cost vector

        ConvergenceCurves=zeros(Max_SubIteration,Max_iteration);

        %BPSO with s-shaped family of transfer functions 
        TranFunOption = 4; % 1 to 4 it is a s-shaped
                           % 4 to 8 it is a v-shaped
        for it = 1:Max_SubIteration
            [gBest1{it,count1}, gBestScore1(it,count1) ,ConvergenceCurves(it,:)] = BPSO(noP,Max_iteration,TranFunOption,CostFunction,noV);
        end
        
        InitOpt{count1} = unique(cell2mat(gBest1(:,count1)),'rows');
        %Display the initiatives
        for i = 1:size(InitOpt{count1},1)
            disp(['Budget: ' num2str(Budget) ', Initiatives: ' cell2mat(model.Initiative(InitOpt{count1}(i,:)~=0)') ', Score= ' num2str(1/gBestScore1(count1))])
        end
    end
    toc

end
%% 

figure(1)
col=hsv(size(optimunBeta,2));
count=0;
marker = ['o','s','d','*','v'];
for i= DepLevel
    count = count+1;
    ydata =optimunBeta(:,i);
    xdata =BudgetArray;
    p1=plot(xdata,ydata,'Color',col(count,:));
    p1.LineStyle  = '-'; p1.LineWidth = 0.5; p1.Marker = marker(count); 
    p1.MarkerSize = 8; p1.MarkerFaceColor = col(count,:);
    leg1{count} = ['Dep Level  = ' num2str(DepLevel(count) - 1) ]; 
    hold on
end
    
set(gca,'FontSize',14); %set the axis font size
%set(gca,'xtick',0:0.25:1.5)
box on;
grid on;
xlabel('Cost Bound', 'FontSize', 18)
ylabel('Optimal Score', 'FontSize', 18)
legend(leg1,'FontSize',18,'Location','southeast');
hold off


% save ('BetaModel.mat','BudgetArray','optimunBeta')
%% 
matrixSoluction = cell2mat(gBest1(:,count1));
IniMatrix = cell(Max_SubIteration,1);
for i=1:Max_SubIteration
    IniMatrix{i} = cell2mat(model.Initiative(matrixSoluction(i,:)~=0)');
end
disp(['Initiatives: ' unique(IniMatrix)'])