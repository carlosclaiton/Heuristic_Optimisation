%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This functions implements the PSEQEA algorithm proposed by 
% Hossain,Kowsar, et al published on Proceedings of 13th International 
% Conference on Computer and Information Technology (ICCIT 2010)
% 23-25 December, 2010, Dhaka, Bangladesh.
%
%Developer: Carlos Kuhn
%contact: Carlos.Kuhn@dst.defence.gov.au
%Date: 24/03/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gBest,gBestScore,ConvergenceCurve,IterNumb]=PSEQEA(noP,Max_iteration,CostFunction,noV)

% %%\\\\\\\\\\\\\\\\\\\\\\\Entre data for debuging\\\\\\\\\\\\\\\\\\\\\\\\\\
% clear all
% noP = 30; % number of q-bits
% noV = 15; % number of states
% Max_iteration = 1000;
% Budget = 2;
% DPI = 1;
% model=CreateModel_CK(Budget,DPI); % here is where the model is stablished 
% CostFunction=@(x) MyCost_CK(x,model);

%/////////////////weight variable////////////////////////////////////////
c1=2;       % cognitive weight
c2=2;       % social weight

%///////////////////Angle for the Quantum operator///////////////////////
Theta = zeros(noP,noV); %angle will replace Velocity vector when using QEA

%////////////////initial q-bits all with the same propability/////////////
qbit = ones(2,noP,noV)./sqrt(2);

%////////////////initialize Position vector//////////////////////////////
Position=randi([0,1],noP,noV);

%/////////////////////////Cognitive component////////////////////////////
pBestScore=zeros(noP);
pBest=zeros(noP,noV);

%//////////////////////////Social component///////////////////////////////
gBestScore=inf;
gBest=zeros(1,noV);

%/////////////////////////////Q-bit opertor///////////////////////////////
OpQ =@(x)[cos(x), -sin(x); sin(x), cos(x)]; 

%/////////////////////////////Convergency/////////////////////////////////
ConvergenceCurve=zeros(1,Max_iteration); %Convergence vector

%/////////////////////////////////////////////////////////////////////////
%Calculate cost for each particle
NumbConv = 200;
IterNumb=1;
stepEqual = 0;

while IterNumb < Max_iteration % stepEqual < NumbConv && IterNumb < Max_iteration  
    
    %Calculate the cost for each state
    for i=1:size(Position,1)  
        [fitness, ~]=CostFunction(Position(i,:));
        
        if(pBestScore(i)>fitness)
            pBestScore(i)=fitness;
            pBest(i,:)=Position(i,:);
        end
        if(gBestScore>fitness)
            gBestScore=fitness;
            gBest=Position(i,:);
        end
    end

    %Calculating the rotation angle for each q-bit
     for i=1:size(Position,1)
        for j=1:size(Position,2) 
            Theta(i,j)= c1*(pBest(i,j)-Position(i,j))+c2*(gBest(j)-Position(i,j));
        end
     end
     %updating the q-bits
     for i=1:size(Position,1)
        for j=1:size(Position,2) 
            qbit(:,i,j)= OpQ(Theta(i,j))*qbit(:,i,j);
        end
     end
     
     %updating the states
    for i=1:size(Position,1) % For each particle
        for j=1:size(Position,2) % For each variable
            if rand < qbit(2,i,j)^2
                Position(i,j)=1;
            else
                Position(i,j)=0;
            end
        end
    end
    if ConvergenceCurve(IterNumb) - gBestScore < 10^-10
        stepEqual = stepEqual+ 1;
    else
        stepEqual = 0; 
    end
    
    if stepEqual > NumbConv % || IterNumb > Max_iteration-1
        break
    end
   
    ConvergenceCurve(IterNumb)=gBestScore;
    IterNumb=IterNumb+1;

%     stepEqual = sum(ConvergenceCurve == gBestScore);
end
end

