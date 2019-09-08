function [basis,normalizedBasis]=C4mult(level,first,second,useData,Data,varargin) 
%Inputs: int level, arrays of first, second, logical useData, struct Data
%Optional inputs: array=varargin{1}
%Outputs: arrays basis and normalizedBasis

%Description: Computes the generator GC in coordinates of first=(k1,n1,m1) on the given level.
%Same for GD and second=(k2,n2,m2). The level is the same for both so we can multiply.
%So then we take the product GC*GD and write it as a linear combination of generators and
%store the coefficients in basis.
%The normalizedBasis is equal to the basis modulo signs.

%If there are multiple generators in the coordinates of GC (non cyclic answer) we specify
%which generator we want through varargin{1}=[x,0] where x is the index of the generator in our list of generators 

%Similarly if we have multiple generators at GD we specify them via
%varargin{2}=[0,y]. Finally if there are multiple generators at both GC
%and GD then we use varargin{1}=[x,y]


if isequal(first,zeros(1,3)) || isequal(second,zeros(1,3)) %We multiply with [0,0,0] i.e. 1. No point in any more computations
    basis=1;
    normalizedBasis=1;
    return
end
    
k1=first(1); n1=first(2); m1=first(3); k2=second(1); n2=second(2); m2=second(3);

%Reindex appropriately
k1=C4kreindex(k1,n1,m1);
k2=C4kreindex(k2,n2,m2);

%Basic check
if k1<0 || k1>abs(n1)+2*abs(m1) || k2<0 || k2>abs(n2)+2*abs(m2)
    basis=[];
    normalizedBasis=[];
    return
end

%Get our chains
[rankC,C]=C4allChains(n1,m1,useData,Data);
[rankD,D]=C4allChains(n2,m2,useData,Data);

%Get generator and homology for both. Return empty basis if generator=0
[GC{level},HC{level}]=Homology(C{level}{k1+2},C{level}{k1+1});
if isequal(HC{level},0)
    basis=[];
    normalizedBasis=[];
    return
elseif size(HC{level},2)>1  %If multiple generators, use varargin to decide which to use
    if isempty(varargin)
        disp(HC{level})
        error('Please enter which generator x of the above you want as [x,0] where x=1 if it is the first option etc. ')
    else
        GC{level}=GC{level}(:,varargin{1}(1));
    end
end
[GD{level},HD{level}]=Homology(D{level}{k2+2},D{level}{k2+1});
if isequal(HD{level},0)
    basis=[];
    normalizedBasis=[];
    return
elseif size(HD{level},2)>1
    if isempty(varargin)
        disp(HD{level})
        error('Please enter which generator of the above you want as [0,y] where y=1 if it is the first option etc.')
    else
        GD{level}=GD{level}(:,varargin{1}(2));
    end
end

%In the product we are interested at index k1+k2. So we have two
%differentials D0,D1 and we need to transfer them so we also need 
%the ranks at k1+k2-1, k1+k2 and k1+k2+1. We call these low, mid, high
%resp. 

rankMid=cell(1,4); %We will only transfer this 

rankLow=rankBox(k1+k2-1,rankC{1},rankD{1}); %Needed to transfer D0 
[rankMid{1},detailedrankMid]=rankBox(k1+k2,rankC{1},rankD{1}); %Where our generator lives. We need the detailedrank so as to pad appropriately. See below
rankHigh=rankBox(k1+k2+1,rankC{1},rankD{1}); %Needed to transfer D1 


D0{1}=BoxDiff(k1+k2,rankC{1},rankD{1},C{1},D{1},useData,Data); %Exiting Differential at k1+k2
D1{1}=BoxDiff(k1+k2+1,rankC{1},rankD{1},C{1},D{1},useData,Data); %Entering Differential at k1+k2

%Now we restrict our generators to the bottom, take the products of the restrictions and
%finally invert restrictions so as to get the product at the correct level

%Restrict our generators to the bottom, so as to be able to multiply them (we can only multiply in the equivariant bases there)
lvl=level;
while lvl>1
    GC{lvl/2}=restrict(GC{lvl},rankC{lvl}{k1+1},rankC{lvl/2}{k1+1});
    GD{lvl/2}=restrict(GD{lvl},rankD{lvl}{k2+1},rankD{lvl/2}{k2+1});
    lvl=lvl/2;
end
    
%Now multiply on the bottom. To do that we need to pad left and right by
%zeros (think about how the tensor of chain complexes works)

padleft=0;
for j=0:k1-1  
    padleft=padleft+sum(detailedrankMid{j+1});
end
padright=0;
for j=k1+1:k1+k2
    padright=padright+sum(detailedrankMid{j+1});
end

%The product on the bottom in the left convenient basis
leftconvproduct=GC{1}.*GD{1}(:)'; 
leftconvproduct=leftconvproduct(:);  
%I.e. multiply the GC{1} with each entry of GD{1} and then concatenate vertically.


%Then we change the basis
[convlefttocanon, ~]=boxchangebasis(rankC{1}{k1+1},rankD{1}{k2+1},useData,Data);
productcanon=convlefttocanon*leftconvproduct;
%Finally we pad it
product=cell(1,4);
product{1}=[zeros(padleft,1);productcanon;zeros(padright,1)];
%This is the product of the restrictions on the bottom level


%Invert restrictions. Remember l was used to go from level to 1. Now we do the opposite

while lvl<level
    rankMid{2*lvl}=ranktransfer(rankMid{lvl},lvl);
    product{2*lvl}=invres(product{lvl},rankMid{lvl},lvl);
    lvl=2*lvl;
end

%Now transfer the chains upstairs (no need to transfer the detailedrank, that's only for padding). 
%We don't transfer the product, rather use the inverse of restriction as we want to get the product, not a multiple of it (transfer of restriction). 

D0{level}=transferdifferential(D0{1},4/level,rankMid{1},rankLow);
D1{level}=transferdifferential(D1{1},4/level,rankHigh,rankMid{1});

[~,Homologyatproduct,SmithVariables]=Homology(D1{level},D0{level});
q=Homologyelement(product(level),SmithVariables);
basis=q{1};

%We normalize the basis
normalizedBasis=abs(basis);
if isequal(Homologyatproduct,4) && basis==3
    normalizedBasis=1;
end
end