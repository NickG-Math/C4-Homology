function [basis,normalBasis]=Multiply(order,level,k1,k2,rankC,rankD,C,D,GC,GD,useData,Data) 
%
%[basis,normalBasis]=Multiply(order,level,k1,k2,rankC,rankD,C,D,GC,GD,useData,Data) 
%
%INPUT: ints order,level,k1,k2, cells of cells of arrays rankC,rankD, cells C,D, columns GC,GD, logical useData, struct Data
%
%OUTPUT: arrays basis, normalBasis
%
%DESCRIPTION: Given two chains C,D and generators GC,GD of their homology
%at k1,k2 and level, write their product as a linear combination of the
%generators at k1+k2+1 and store them in basis
%
%normalBasis is basis normalized: if basis(i) generates the i-th homology
%then normalBasis(i)=1
%
%Note: C,D are given on bottom level but rankC,rankD are given on bottom level
%and "level" level (we pretransferred them to avoid repeat calculations). 
%GC,GD are given on "level" level


%In the product we are interested at index k1+k2. 
%So we have two differentials D0,D1 and we need to transfer them so we also need the ranks at k1+k2-1, k1+k2 and k1+k2+1. 
%We call these low, mid, high ranks resp. 

rankLow=rankBox(k1+k2-1,rankC{1},rankD{1}); %Needed to transfer D0 
[rankMid,detailedrankMid]=rankBox(k1+k2,rankC{1},rankD{1}); %Where our generator lives. We need the detailedrank so as to pad appropriately. See below
rankHigh=rankBox(k1+k2+1,rankC{1},rankD{1}); %Needed to transfer D1 

%Restrict our generators to the bottom, so as to be able to multiply them (we can only multiply in the equivariant bases there)
ResGC=restrict(GC,rankC{level}{k1+1},rankC{1}{k1+1});
ResGD=restrict(GD,rankD{level}{k2+1},rankD{1}{k2+1});
    
%Multiply on the bottom. To do that we need to pad left and right by zeros (think about how the tensor of chain complexes works)
padleft=0;
for j=0:k1-1  
    padleft=padleft+sum(detailedrankMid{j+1});
end
padright=0;
for j=k1+1:k1+k2
    padright=padright+sum(detailedrankMid{j+1});
end

%The product on the bottom in the left convenient basis is [GC*GD(1);GC*CD(2);....]
%Vectorized:
leftconvproduct=ResGC.*ResGD(:)'; 
leftconvproduct=leftconvproduct(:);  

%Then we change the basis
[convlefttocanon, ~]=boxchangebasis(rankC{1}{k1+1},rankD{1}{k2+1},useData,Data);
productcanon=convlefttocanon*leftconvproduct;

%Finally we pad it
ResProduct=[zeros(padleft,1);productcanon;zeros(padright,1)];

%This is the product of the restrictions on the bottom level. 
%To get the actual product we invert the restrictions. 

rankMidtop=ranktransfer(rankMid,level,order); %Transfer rankMid to be used below
product=invres(ResProduct,rankMidtop,rankMid);
%So finally we have the product of the generators at our level!

%Now we get the differentials of the box product.
D0_Bottom=BoxDiff(k1+k2,rankC{1},rankD{1},C,D,useData,Data); %Exiting Differential at k1+k2
D1_Bottom=BoxDiff(k1+k2+1,rankC{1},rankD{1},C,D,useData,Data); %Entering Differential at k1+k2

%We transfer the boxed differentials upstairs (no need to transfer the detailedrank, that's only been used for padding). 
D0_Level=transferdifferential(D0_Bottom,order/level,rankMid,rankLow);
D1_Level=transferdifferential(D1_Bottom,order/level,rankHigh,rankMid);

%Finally we compute the basis
[~,Homologyatproduct,SmithVariables]=Homology(D1_Level,D0_Level);
q=Homologyelement({product},SmithVariables);
basis=q{1};

%We normalize the basis
normalBasis=abs(basis);
toBeNormalized=Homologyatproduct~=1 & basis~=0; %We don't want to change Z and the second condition is because gcd(0,?)=0 instead of ?.
normalBasis(toBeNormalized)=gcd(basis(toBeNormalized),Homologyatproduct(toBeNormalized));

end

