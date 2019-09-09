function [basis,normalizedBasis]=C4mult(level,first,second,useData,Data,varargin) 
%Inputs: int level, arrays of first, second, logical useData, struct Data
%Optional inputs: array varargin{1}=[x,y]
%Outputs: arrays basis and normalizedBasis
%Description: First computes the generator GC in coordinates of first=(k1,n1,m1) on the given level. Same for GD and second=(k2,n2,m2). The level must be the same for both so we can multiply.
%Then takes the product GC*GD and writes it as a linear combination of generators and stores the coefficients in basis. 
%The basis is empty if GC=0 or GD=0. The normalizedBasis is equal to the basis modulo signs.
%If there are multiple generators in the coordinates of GC (non cyclic answer) we specify which generator we want through varargin{1}=[x,0] where x=1,2,... corresponds to the first,second,... generator 
%Similarly if we have multiple generators at GD we specify them via varargin{2}=[0,y]. Finally if there are multiple generators at both GC and GD then we use varargin{1}=[x,y]

if isequal(first,zeros(1,3)) || isequal(second,zeros(1,3)) %We multiply with [0,0,0] i.e. 1. No point in any more computations
    basis=1;
    normalizedBasis=1;
    return
end
    
k1=first(1); n1=first(2); m1=first(3); k2=second(1); n2=second(2); m2=second(3);

%Reindex appropriately
k1=C4kreindex(k1,n1,m1);
k2=C4kreindex(k2,n2,m2);

%We will need these a few times
max1=abs(n1)+2*abs(m1);
max2=abs(n2)+2*abs(m2); 

%Basic check
if k1<0 || k1>max1 || k2<0 || k2>max2
    basis=[];
    normalizedBasis=[];
    return
end

%We now compute the chains C,D that the given generators live in.
%We will only need to transfer them at k1,k2 resp but we preallocate all levels 1,2,4 for indexing convenience
rankC=cell(1,4); rankC{1}=cell(1,max1+2); rankC{2}=rankC{1}; rankC{4}=rankC{1};
C=rankC;
rankD=cell(1,4); rankD{1}=cell(1,max2+2); rankD{2}=rankD{1}; rankD{4}=rankD{1};
D=rankD;

%We get these chains up to k1+k2+1 since we need the differential of the their Box at k1+k2
for i=0:min(k1+k2+1,max1)
        [rankC{1}{i+1},~,C{1}{i+1}]=C4Diff(i,n1,m1,useData,Data);
end
for i=0:min(k1+k2+1,max2)
    [rankD{1}{i+1},~,D{1}{i+1}]=C4Diff(i,n2,m2,useData,Data);
end

%Now we get the generators of C,D at k1,k2. First we need to transfer C,D at k1,k2 up to the given level
lvl=1;
while lvl<level
    rankC{2*lvl}{k1+1}=ranktransfer(rankC{lvl}{k1+1},lvl);
    rankD{2*lvl}{k2+1}=ranktransfer(rankD{lvl}{k2+1},lvl);
    lvl=2*lvl;
end

%We transfer the exiting differential
if k1>0 %No point transferring empties
    C{level}{k1+1}=transferdifferential(C{1}{k1+1},4/level,rankC{1}{k1+1},rankC{1}{k1});
end
if k2>0 
    D{level}{k2+1}=transferdifferential(D{1}{k2+1},4/level,rankD{1}{k2+1},rankD{1}{k2});
end

%We transfer the entering differential
C{level}{k1+2}=transferdifferential(C{1}{k1+2},4/level,rankC{1}{k1+2},rankC{1}{k1+1});
D{level}{k2+2}=transferdifferential(D{1}{k2+2},4/level,rankD{1}{k2+2},rankD{1}{k2+1});

%Get the generator for C at k1 and level. Return empty basis if generator=0
[GC{level},HC{level}]=Homology(C{level}{k1+2},C{level}{k1+1});
if isequal(HC{level},0)
    basis=[];
    normalizedBasis=[];
    return
elseif size(HC{level},2)>1  %If multiple generators exist, use varargin to decide which to use
    if isempty(varargin)
        disp(HC{level})
        error('Multiple generators in the first sphere; please specify which of the above to use')
    else
        GC{level}=GC{level}(:,varargin{1}(1));
    end
end

%Get the generator for D at k2 and level. Return empty basis if generator=0
[GD{level},HD{level}]=Homology(D{level}{k2+2},D{level}{k2+1});
if isequal(HD{level},0)
    basis=[];
    normalizedBasis=[];
    return
elseif size(HD{level},2)>1 %If multiple generators exist, use varargin to decide which to use
    if isempty(varargin)
        disp(HD{level})
        error('Multiple generators in the second sphere; please specify which of the above to use')
    else
        GD{level}=GD{level}(:,varargin{1}(2));
    end
end


%In the product we are interested at index k1+k2. So we have two 
%differentials D0,D1 and we need to transfer them so we also need 
%the ranks at k1+k2-1, k1+k2 and k1+k2+1. We call these low, mid, high ranks resp. 

rankMid=cell(1,4); %We will only transfer this rank

rankLow=rankBox(k1+k2-1,rankC{1},rankD{1}); %Needed to transfer D0 
[rankMid{1},detailedrankMid]=rankBox(k1+k2,rankC{1},rankD{1}); %Where our generator lives. We need the detailedrank so as to pad appropriately. See below
rankHigh=rankBox(k1+k2+1,rankC{1},rankD{1}); %Needed to transfer D1 


D0{1}=BoxDiff(k1+k2,rankC{1},rankD{1},C{1},D{1},useData,Data); %Exiting Differential at k1+k2
D1{1}=BoxDiff(k1+k2+1,rankC{1},rankD{1},C{1},D{1},useData,Data); %Entering Differential at k1+k2

%Restrict our generators to the bottom, so as to be able to multiply them (we can only multiply in the equivariant bases there)
%Remember lvl=level now
while lvl>1
    GC{lvl/2}=restrict(GC{lvl},rankC{lvl}{k1+1},rankC{lvl/2}{k1+1});
    GD{lvl/2}=restrict(GD{lvl},rankD{lvl}{k2+1},rankD{lvl/2}{k2+1});
    lvl=lvl/2;
end
    
%Multiply on the bottom. To do that we need to pad left and right by zeros (think about how the tensor of chain complexes works)
padleft=0;
for j=0:k1-1  
    padleft=padleft+sum(detailedrankMid{j+1});
end
padright=0;
for j=k1+1:k1+k2
    padright=padright+sum(detailedrankMid{j+1});
end

%The product on the bottom in the left convenient basis is [GC{1}*GD{1}(1);GC{1}*CD{1}(2);....]
%Vectorized:
leftconvproduct=GC{1}.*GD{1}(:)'; 
leftconvproduct=leftconvproduct(:);  

%Then we change the basis
[convlefttocanon, ~]=boxchangebasis(rankC{1}{k1+1},rankD{1}{k2+1},useData,Data);
productcanon=convlefttocanon*leftconvproduct;

%Finally we pad it
product=cell(1,4);
product{1}=[zeros(padleft,1);productcanon;zeros(padright,1)];

%This is the product of the restrictions on the bottom level. 
%To get the actual product we invert the restrictions. 

%Remember lvl=1 now
while lvl<level
    rankMid{2*lvl}=ranktransfer(rankMid{lvl},lvl);
    product{2*lvl}=invres(product{lvl},rankMid{lvl},lvl);
    lvl=2*lvl;
end
%So finally we have the product of the generators at our level!

%Now transfer the boxed differentials upstairs (no need to transfer the detailedrank, that's only been used for padding). 

D0{level}=transferdifferential(D0{1},4/level,rankMid{1},rankLow);
D1{level}=transferdifferential(D1{1},4/level,rankHigh,rankMid{1});

%Finally we compute the basis
[~,Homologyatproduct,SmithVariables]=Homology(D1{level},D0{level});
q=Homologyelement(product(level),SmithVariables);
basis=q{1};

%We normalize the basis
normalizedBasis=abs(basis);
if isequal(Homologyatproduct,4) && basis==3
    normalizedBasis=1;
end
end