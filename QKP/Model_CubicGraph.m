% this code will generate the data that is simulate a cubic graph

function model=Model_CubicGraph(Budget,NumNodes)

% clc   
% NumNodes=16; %number of initiatives

    Edge = zeros(NumNodes);
    Cost = ones(NumNodes,1);
    Value = ones(NumNodes,1);

    Edge(1,end)=1;
    Edge(end,1)=1;
    
    for i=1:NumNodes-1
        Edge(i,i+1) = 1;
        Edge(i+1,i) = 1;
%         [i,i+1]
%         [i+1,i]        
    end
    
    for i=1:NumNodes/4
        Edge(i,i+NumNodes/4) = 1;
        Edge(i+NumNodes/4,i) = 1;
%         [i,i+NumNodes/4]
%         [i+NumNodes/4,i]        
    end
    
    for i=NumNodes/2+1:3*NumNodes/4
        Edge(i,i+NumNodes/4) = 1;
        Edge(i+NumNodes/4,i) = 1;
%         [i,i+NumNodes/4]
%         [i+NumNodes/4,i]        
    end
    
%%THIS IS WHAT THE MODEL WILL PASS THROUGH
    model.Budget = Budget; 
    model.Cost=Cost;
    model.Value =Value;
    model.Edge=Edge;
end