function Sorted=sortedDistance(gen,ProdOrDiv,BasicIrreducibles,DistanceAgainst)
%
%Sorted=sortedDistance(gen,ProdOrDiv,BasicIrreducibles,DistanceAgainst)
%
%INPUTS: Array gen, logical ProdOrDiv, cells of arrays BasicIrreducibles,
%DistanceAgainst
%
%OUPUTS: Array sorted
%
%DESCRIPTION: Given the generator, it computes the newgen that will be used
%in FactorizeNum or FactorizeDen (ProdOrDiv determines which one). Then the
%distances of newgen to the BasicIrreducibles are computed and the minimum
%is found. Finally the newgens are sorted in terms of that minimum
%distance, in ascending order (to be precise we order the index of the
%BasicIrreducibles, not the newgens themselves)


minDistance=zeros(1,size(BasicIrreducibles,2));
for i=1:size(BasicIrreducibles,2)
    if ProdOrDiv==1 %product
        newgen=gen-BasicIrreducibles{i};
    else %Division
        newgen=gen+BasicIrreducibles{i};
    end
    minDistance(i)=sum(abs(newgen));
    if minDistance(i)==0
        Sorted=i;
        return
    end
    for j=1:size(DistanceAgainst,2)
        a=sum(abs(newgen-DistanceAgainst{j}));
        if a==0
            Sorted=i;
            return
        elseif a<minDistance(i)
            minDistance(i)=a;
        end
    end
end
[~,Sorted]=sort(minDistance);
end