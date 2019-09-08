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

%First we get our chains
[rankC,C]=C4allChains(n1,m1,useData,Data);
[rankD,D]=C4allChains(n2,m2,useData,Data);
%Reindex appropriately
k1=C4kreindex(k1,n1,m1);
k2=C4kreindex(k2,n2,m2);

%Basic check
if k1<0 || k1>abs(n1)+2*abs(m1) || k2<0 || k2>abs(n2)+2*abs(m2)
    basis=[];
    normalizedBasis=[];
    return
end

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

%Box at the bottom (you can only box there)
Boxed=cell(1,4); totalrank=cell(1,4);
[totalrank{1},detailedrank{1},Boxed{1}]=Box(rankC{1},rankD{1},C{1},D{1},useData,Data); 
%We need the detailedrank so as to pad appropriately. See below

%Restrict our generators to the bottom, so as to be able to multiply them (we can only multiply in the equivariant bases there)
l=level;
while l>1
    GC{l/2}=restrict(GC{l},rankC{l}{k1+1},rankC{l/2}{k1+1});
    GD{l/2}=restrict(GD{l},rankD{l}{k2+1},rankD{l/2}{k2+1});
    l=l/2;
end
    
%Now multiply on the bottom. To do that we need to pad left and right by
%zeros (think about how the tensor of chain complexes works)

padleft=0;
for j=0:k1-1  
    padleft=padleft+sum(detailedrank{1}{k1+k2+1}{j+1});
end
padright=0;
for j=k1+1:k1+k2
    padright=padright+sum(detailedrank{1}{k1+k2+1}{j+1});
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


%Now transfer the chains upstairs (no need to transfer the detailedrank, that's only for padding). 
%We don't transfer the product, rather use the inverse of restriction as we want to get the product, not a multiple of it (transfer of restriction). 
%Remember l was used to go from level to 1. Now we do the opposite

Boxed{level}=cell(1,k1+k2+2);
for i=[k1+k2+1,k1+k2+2] %Only transfer the differential where we need it
    if i>1
        Boxed{level}{i}=transferdifferential(Boxed{1}{i},4/level,totalrank{1}{i},totalrank{1}{i-1});
    end
end

        
while l<level
    totalrank{2*l}{k1+k2+1}=ranktransfer(totalrank{l}{k1+k2+1},l);
    product{2*l}=invres(product{l},totalrank{l}{k1+k2+1},l);
    l=2*l;
end

[~,Homologyatproduct,SmithVariables]=Homology(Boxed{level}{k1+k2+2},Boxed{level}{k1+k2+1});
q=Homologyelement(product(level),SmithVariables);
basis=q{1};

%We normalize the basis
normalizedBasis=abs(basis);
if isequal(Homologyatproduct,4) && basis==3
    normalizedBasis=1;
end
end