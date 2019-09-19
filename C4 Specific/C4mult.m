function [basis,normalBasis]=C4mult(level,first,second,useData,Data,varargin) 
%
%[basis,normalBasis]=C4mult(level,first,second,useData,Data,varargin) 
%
%INPUT: int level, arrays first and second, logical useData, struct Data
%
%OPTIONAL INPUT: array varargin{1}=[x,y]
%
%%OUTPUT: arrays basis and normalBasis
%
%DESCRIPTION: C4mult computes the generators GC, GD in coordinates first=(k1,n1,m1) and second=(k2,n2,m2) resp, on the given level.
%
%Then takes the product GC*GD and writes it as a linear combination of generators at (k1+k2,n1+n2,m1+m2) and stores the coefficients in variable basis. 
%
%The basis is empty if GC=0 or GD=0. The normalBasis is equal to the basis modulo signs.
%
%If there are multiple options for GC (non cyclic homology) we specify which generator we want through varargin{1}=[x,0] where x=1,2,... corresponds to the first,second,... generator 
%
%Similarly if we have multiple options GD we specify which we want via varargin{2}=[0,y]. Finally if there are multiple options for both GC and GD then we use varargin{1}=[x,y]

if isequal(first,zeros(1,3)) || isequal(second,zeros(1,3)) %We multiply with [0,0,0] i.e. 1. No point in any more computations
    basis=1;
    normalBasis=1;
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
    normalBasis=[];
    return
end

%We now compute the chains C,D that the given generators live in.
limit1=min(k1+k2+1,max1); limit2=min(k1+k2+1,max2); %How far our chains go

%We preallocate our chains and ranks. We will need to store the transfers
%of our ranks but not the chains, so we preallocate accordingly.

rankC=cell(1,level); C=cell(1,limit1+1); rankD=cell(1,level); D=cell(1,limit2+1); 
%We get the chains up to k1+k2+1 since we need the differential of their Box at k1+k2
for i=0:limit1
        [rankC{1}{i+1},~,C{i+1}]=C4Diff(i,n1,m1,useData,Data);
end
for i=0:limit2
    [rankD{1}{i+1},~,D{i+1}]=C4Diff(i,n2,m2,useData,Data);
end

%Now we get the generators of C,D at k1,k2 and level. First is C
[HC_Levels,GC_Levels,rankC_Levels]=C4Homology(level,k1,n1,m1,useData,Data);
HC=HC_Levels{level}; GC=GC_Levels{level}; rankC{level}{k1+1}=rankC_Levels{level};
if ~any(HC) %If homology is trivial return empty basis
    basis=[];
    normalBasis=[];
    return
elseif size(HC,2)>1  %If multiple generators exist, use varargin to decide which to use
    if isempty(varargin)
        disp(HC)
        error('Multiple generators in the first sphere; please specify which of the above to use as final argument [x,0]')
    else
        GC=GC(:,varargin{1}(1));
    end
end
%Now it's the D's turn
[HD_Levels,GD_Levels,rankD_Levels]=C4Homology(level,k2,n2,m2,useData,Data);
HD=HD_Levels{level}; GD=GD_Levels{level}; rankD{level}{k2+1}=rankD_Levels{level};
if ~any(HD)
    basis=[];
    normalBasis=[];
    return
elseif size(HD,2)>1 %If multiple generators exist, use varargin to decide which to use
    if isempty(varargin)
        disp(HD)
        error('Multiple generators in the second sphere; please specify which of the above to use  as final argument [0,y]')
    else
        GD=GD(:,varargin{1}(2));
    end
end

%We have our generators GC and GD. We now take their product
[basis,normalBasis]=Multiply(4,level,k1,k2,rankC,rankD,C,D,GC,GD,useData,Data);
end