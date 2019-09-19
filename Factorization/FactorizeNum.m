function [gen,product,CounterNum,CounterDen,found,VisitedNum,VisitedDen]=FactorizeNum(gen,SwitchLimit,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,product,CounterNum,CounterDen,VisitedNum,VisitedDen,kEnd)
%[product,CounterNum,CounterDen,found,VisitedNum,VisitedDen]=FactorizeDen(gen,SwitchLimit,NumOrDen,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,product,CounterNum,CounterDen,VisitedNum,VisitedDen,kStart)
%
%INPUTS: Same as FactorizeHub
%
%OUTPUTS: Same as FactorizeHub
%
%DESCRIPTION: FactorizeNum tries to factor the numerator by going through the
%table recursively. If it doesn't find it (maybe the numerator is 1), 
%it goes back to the hub and asks it to try to factor the denominator instead


if ~any(gen) %If gen=[0,0,0] then we are done
    found=1;
    return
end
    
for i=1:size(MoreIrreducibles,2)
    if isequal(gen,MoreIrreducibles{i})
        found=1;
        CounterNum(i)=CounterNum(i)+1;
        return
    end
end

        

gen_PosInd=PositiveIndexer(gen);
commaSeparatedList=num2cell(gen_PosInd);

if safety
    if VisitedNum(commaSeparatedList{:},kEnd)==1
        found=0;
        return
    else
        VisitedNum(commaSeparatedList{:},kEnd)=1;
    end
end


index=sortedDistance(gen,1,BasicIrreducibles,BasicIrreducibles);
%index=1:size(BasicIrreducibles,2);

for i=index
    RelevantPart=Table{i,commaSeparatedList{:}}; %This can be empty or a cell
    if isempty(RelevantPart)
        continue
    end
    newgen=gen-BasicIrreducibles{i};
    
    for kStart=1:size(RelevantPart,2)
        if RelevantPart{kStart}(kEnd)~=1 
            continue
        end
        %Else
        Counternew=CounterNum;
        Counternew(i)=Counternew(i)+1;
        [genprov,productprov,CounterNumprov,CounterDenprov,found,VisitedNum,VisitedDen]=FactorizeNum(newgen,SwitchLimit,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,product,Counternew,CounterDen,VisitedNum,VisitedDen,kStart);
        if found==1
            gen=genprov;
            product=productprov;
            CounterNum=CounterNumprov;
            CounterDen=CounterDenprov;
            return
        end
         %Else continue to search with the previous gen, Counter
    end
end

%If you've searched through thick and thin but didn't find, go back to the Hub

[productprov,CounterNumprov,CounterDenprov,found,VisitedNum,VisitedDen]=FactorizationHub(gen,SwitchLimit,-1,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,product,CounterNum,CounterDen,VisitedNum,VisitedDen,kEnd);
if found==1
    product=productprov;
    CounterNum=CounterNumprov;
    CounterDen=CounterDenprov;
    return
else
    return
end

end