% Budget Scenario is a decision making problem, then each decision varible
% is a binary variable (0,1), then if we use the uniform distribution for
% such in matlab
clc
xDist = randi([0,1],[100,1]); %we have 15 Initiatives
disp(mean(xDist))
disp(std(xDist))

% as we can see above the mean of such distribution is 0.5 and the standart
% deviation is also 0.5, when we are using a uniform distribution.

%% Now we need to know that distribution is inside of another one, the fact
%we can choose different combinations of a set of option, in this
%particular example in N. Order data we have 15 initiatives, then we can
%check the number of way we have to peak
nIn = 15;
NumberWay = ones(nIn,1)*NaN;
for i=1:nIn
    NumberWay(i) = nchoosek(nIn,i);
end
nEle=1:15;

figure()
plot(NumberWay);
xlabel('number of elements in the subset', 'FontSize',16)
ylabel('number of possible combinations', 'FontSize',16)
grid on
box on

%GETTING THE FIRST MOMENT OF THE DISTRIBUTION
DistMean = trapz(nEle,nEle.*NumberWay')/trapz(nEle,NumberWay');

% the valuo of DistMean, shows that at 7.5 as the number of elements on the
% subset has a higher mode, I would expected that it would be more likely
% to have subsets with 7 or 8 elements as answer for the BS if we dont have
% bounds. 

%% I am trying to understand how this would affect the heuristic model and 
% how the statistics of the problem is to extract from the model the
% features to be use to training a Machine Learning algorithm

Budget = 1:30; % boundary on the cost
DPI = 1;    % level of dependence of initiatives on the scenarios
numRun = 10000;
Fitness = ones(numRun,length(Budget))*NaN;

for j =1:length(Budget)
    model=CreateModel_CK(Budget(j),DPI); % here is where the model is stablished 
    noV = length(model.Cost); % number of decision variables, match with the size of cost vector
    CostFunction=@(x) MyCost_CK(x,model); % Modify or replace Mycost.m according to your cost funciton
    for i=1:numRun
        xVar = randi([0,1],[15,1]);
        Fitness(i,j) = CostFunction(xVar');
    end
end

%get the mean value
MeanFit = ones(length(Budget),1)*NaN;
stdVar = ones(length(Budget),1)*NaN;
for i=1:length(Budget)   
    MeanFit(i) = mean(Fitness(Fitness(:,i)~=Inf,i));
    stdVar(i) = std(Fitness(Fitness(:,i)~=Inf,i));
end
figure()
ydata =MeanFit;
ydataErr = stdVar;
xdata =  Budget;
p1=errorbar(xdata,ydata,ydataErr,ydataErr,'Color','b');
p1.LineStyle  = '-'; p1.LineWidth = 0.5; p1.Marker = 'o'; p1.MarkerSize = 8; p1.MarkerFaceColor = 'b';
xlabel('Budget', 'FontSize',16)
ylabel('Fitness mean ', 'FontSize',16)
grid on
box on



