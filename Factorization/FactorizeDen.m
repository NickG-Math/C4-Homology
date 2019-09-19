function [gen,product,CounterNum,CounterDen,found,VisitedNum,VisitedDen]=FactorizeDen(gen,SwitchLimit,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,product,CounterNum,CounterDen,VisitedNum,VisitedDen,kStart)
%[product,CounterNum,CounterDen,found,VisitedNum,VisitedDen]=FactorizeDen(gen,SwitchLimit,NumOrDen,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,product,CounterNum,CounterDen,VisitedNum,VisitedDen,kStart)
%
%INPUTS: Same as FactorizeHub
%
%OUTPUTS: Same as FactorizeHub
%
%DESCRIPTION: FactorizeDen tries to factor the denominator by going through the
%table recursively. If it doesn't find it (maybe the denominator is 1), 
%it goes back to the hub and asks it to try to factor the nominator instead


%Shortcut: If we get a positive then hand it to FactorizeNum through FactorizeGenerator
if all(gen>=0) 
    [productprov,CounterNumprov,CounterDenprov,found,VisitedNum,VisitedDen]=FactorizationHub(gen,SwitchLimit,1,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,product,CounterNum,CounterDen,VisitedNum,VisitedDen,kStart);
    if found==1
        product=productprov;
        CounterNum=CounterNumprov;
        CounterDen=CounterDenprov;
        return
    end
end


for i=1:size(MoreIrreducibles,2) %If we get an irreducible then it's over
    if isequal(gen,MoreIrreducibles{i})
        found=1; 
        CounterNum(i)=CounterNum(i)+1;
        return
    end
end




gen_PosInd=PositiveIndexer(gen);

if any(gen_PosInd>size(IndexedHomology)) %Hit the walls of the table
   found=0;
   disp('Out of bounds')
   return
end
commaSeparatedList=num2cell(gen_PosInd);

if safety
    if VisitedDen(commaSeparatedList{:},kStart)==1 %Mark our visit
        found=0;
        return
    else
        VisitedDen(commaSeparatedList{:},kStart)=1;
    end
end




Homologystart=IndexedHomology{commaSeparatedList{:}}(kStart); %Get the homology at generator; this is to check if multiplication by x is
%an isomorphism y/x->x


index=sortedDistance(gen,-1,BasicIrreducibles,BasicIrreducibles);

for i=index
    newgen=BasicIrreducibles{i}+gen;
    newgen_PosInd=PositiveIndexer(newgen);
    if any(newgen_PosInd>size(IndexedHomology)) %Within bounds
        continue
    end
    newcommaSeparatedList=num2cell(newgen_PosInd);

    Homologyend=IndexedHomology{newcommaSeparatedList{:}};

    MultEnd=Table{i,newcommaSeparatedList{:}}{kStart};
    
    kEndPossibilities=zeros(1,size(Homologyend,2));
    for k=1:size(Homologyend,2)
        if (Homologystart==1 && Homologyend(k)~=1) || Homologystart>Homologyend(k) || MultEnd(k)==0 %We want an injection for division
            continue
        end
        kEndPossibilities(k)=k;
    end
    kEndPossibilities(kEndPossibilities==0)=[];
    %Now kEndPossibilities has all valid k's that give an injection
    for k=kEndPossibilities
        productnew=product*MultEnd(k);
        Counternew=CounterDen;
        Counternew(i)=CounterDen(i)+1;
        [genprov,productprov,CounterNumprov,CounterDenprov,found,VisitedNum,VisitedDen]=FactorizeDen(newgen,SwitchLimit,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,productnew,CounterNum,Counternew,VisitedNum,VisitedDen,k);
        %Provisional
        if found==1
            gen=genprov;
            product=productprov;
            CounterDen=CounterDenprov;
            CounterNum=CounterNumprov;
            return
        end
        %Else continue to search with the same gen, product and Counter
    end
end


%If you've searched through thick and thin but didn't find, go back to the Hub

[productprov,CounterNumprov,CounterDenprov,found,VisitedNum,VisitedDen]=FactorizationHub(gen,SwitchLimit,1,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,product,CounterNum,CounterDen,VisitedNum,VisitedDen,kStart);
if found==1 
    product=productprov;
    CounterNum=CounterNumprov;
    CounterDen=CounterDenprov;
    return
else
    return
end


end
