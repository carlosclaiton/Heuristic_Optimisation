%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  BPSO and VPSO source codes version 1.0                           %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%       Homepage: http://www.alimirjalili.com                       %
%   Main paper: S. Mirjalili and A. Lewis, "S-shaped versus         %
%               V-shaped transfer functions for binary Particle     %
%               Swarm Optimization," Swarm and Evolutionary         %
%               Computation, vol. 9, pp. 1-14, 2013.                %
%                                                                   %
% History of modifications:                                         %
% Carlos Kuhn, 24/03/2020                                           %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gBest,gBestScore,ConvergenceCurve,IterNumb]=BPSO(noP,Max_iteration,BPSO_num,CostFunction,noV)

%%Initial Parameters for PSO
% w=2;              % Inirtia weight
wMax=0.9;           % Max inirtia weight
wMin=0.4;           % Min inirtia weight
c1=5;               % cognitive weight
c2=5;               % social weight
Vmax=6;


Velocity=zeros(noP,noV);        %Velocity vector
Position=randi([0,1],noP,noV);  %Position vector, Uniformly random O or 1 

%////////Cognitive component///////// 
pBestScore=inf(noP,1);    %initial conditions is infinity
pBest=zeros(noP,noV);
%////////////////////////////////////

%////////Social component///////////
gBestScore=inf;         %initial condition is infinity
gBest=zeros(1,noV);
%///////////////////////////////////
NumbConv = 200;
ConvergenceCurve=zeros(1,Max_iteration); %Convergence vector
IterNumb=1;
stepEqual = 0;
while IterNumb < Max_iteration 
    
    %Calculate cost for each particle
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

    %update the W of PSO (inertia weight)
    w=wMax-IterNumb*((wMax-wMin)/Max_iteration);
    %Update the Velocity and Position of particles
    for i=1:size(Position,1)
        for j=1:size(Position,2) 
            %Equation (1)
            Velocity(i,j)=w*Velocity(i,j)+c1*rand()*(pBest(i,j)-Position(i,j))+c2*rand()*(gBest(j)-Position(i,j));
%             Velocity(i,j)=w*Velocity(i,j)+c1*(pBest(i,j)-Position(i,j))+c2*(gBest(j)-Position(i,j));
            
            if(Velocity(i,j)>Vmax)
                Velocity(i,j)=Vmax;
            end
            if(Velocity(i,j)<-Vmax)
                Velocity(i,j)=-Vmax;
            end  
            
            if BPSO_num==1
                s=1/(1+exp(-2*Velocity(i,j))); %S1 transfer function
            end
%             if BPSO_num==2
%                 s=1/(1+exp(-Velocity(i,j)));   %S2 transfer function              
%             end
            if BPSO_num==2
%                s=1/(1+exp(-(noV/45)*Velocity(i,j)));  %S2 transfer function
               s=1/(1+exp(-(0.013*noV+0.13)*Velocity(i,j)));  %S2 transfer function
            end 
            if BPSO_num==3
                s=1/(1+exp(-Velocity(i,j)/2)); %S3 transfer function              
            end
            if BPSO_num==4
               s=1/(1+exp(-Velocity(i,j)/3));  %S4 transfer function
            end    

            
            if BPSO_num<=4 %S-shaped transfer functions
                if rand<s % Equation (4) and (8)
                    Position(i,j)=1;
                else
                    Position(i,j)=0;
                end
            end
            
            if BPSO_num==5
                s=abs(erf(((sqrt(pi)/2)*Velocity(i,j)))); %V1 transfer function
            end            
            if BPSO_num==6
                s=abs(tanh(Velocity(i,j))); %V2 transfer function
            end            
            if BPSO_num==7
                s=abs(Velocity(i,j)/sqrt((1+Velocity(i,j)^2))); %V3 transfer function
            end            
            if BPSO_num==8
                s=abs((2/pi)*atan((pi/2)*Velocity(i,j))); %V4 transfer function (VPSO)         
            end
            
            if BPSO_num>4 && BPSO_num<=8 %V-shaped transfer functions
                if rand<s %Equation (10)
                    Position(i,j)=~Position(i,j); 
                end
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



