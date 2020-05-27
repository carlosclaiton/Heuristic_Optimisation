%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exhaustive Search, to help to compare time comsumption            %
%                                                                   %
% Carlos Kuhn, 15/05/2020                                           %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Best,BestScore,IterNumb]=BruteForce(Max_iteration,CostFunction,noV)


%%
noV=4;


Position=zeros(noV,2^noV);  %Position vector, Uniformly random O or 1 

BestScore=inf;    %initial conditions is infinity
Best=zeros(noV);

IterNumb=1;

    for i =1:noV
        Position(i,i) = 1;
        
    end
    Position(:,noV+1:2*noV) = circshift(Position(:,1:noV),1,1);

%%
    for i=1:2^noV
            %Get the score
        [fitness, ~]=CostFunction(Position(:,i));
        %Save if the score is better
        if(BestScore>fitness)
            BestScore=fitness;
            Best=Position(i,:);
        end
    end

    
    
end



