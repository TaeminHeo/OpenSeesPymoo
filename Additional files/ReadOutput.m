function [nodeIO, nodeLS, nodeCP, DriftIO, DriftLS, DriftCP, IndexStep]=ReadOutput(DispIO, DispLS, DispCP, nStory)

load node.out -ascii; 
load EQDrift.out -ascii;

[m,n]=size(node);

IndexStep = zeros(1,3);

for i=2:m
    if ( (node(i-1,2)<=DispIO) && (node(i,2)>=DispIO) )  
        IndexStep(1,1) = i;
        nodeIO = node(i,:);
        DriftIO = EQDrift(i,:);
        DriftIO(1,1) = max(EQDrift(i,2:(nStory+1)));
    end
    if ( (node(i-1,2)<=DispLS) && (node(i,2)>=DispLS) )  
        IndexStep(1,2) = i;
        nodeLS = node(i,:);
        DriftLS = EQDrift(i,:);
        DriftLS(1,1) = max(EQDrift(i,2:(nStory+1)));    
    end
    if ( (node(i-1,2)<=DispCP) && (node(i,2)>=DispCP) )  
        IndexStep(1,3) = i;
        nodeCP = node(i,:);
        DriftCP = EQDrift(i,:);
        DriftCP(1,1) = max(EQDrift(i,2:(nStory+1)));  
        tmp = DriftCP*1000;
        DriftCP = round(tmp)/1000;
        break;
    end
end

