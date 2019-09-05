function [basis,varargout]=C4mult(level,first,second,useData,Data,varargin)  %varargin selects the generator if there are multiple ones
%varargout is simplified basis (see the end), homology, generators resp.
varargout=cell(1,3); varargout{2}=cell(1,3); varargout{3}=cell(1,3);

k1=first(1); n1=first(2); m1=first(3); k2=second(1); n2=second(2); m2=second(3);

%Exception when complex is 0->Z->0 so we don't have to fix the chains all
%the time
if isequal([k1,n2,m1,k2,n2,m2],zeros(1,6)) %1 times 1
    basis=1;
    varargout{1}=1;
    return
end

%First we get our chains
[rankC,C]=C4allChains(n1,m1,useData,Data);
[rankD,D]=C4allChains(n2,m2,useData,Data);
k1=C4kreindex(k1,n1,m1);
k2=C4kreindex(k2,n2,m2);

%Basic check
if k1<0 || k1>abs(n1)+2*abs(m1) || k2<0 || k2>abs(n2)+2*abs(m2)
    basis=[];
    return
end

%Get generator and homology for both. Return empty basis if generator=0
[GC{level},HC{level}]=Homology(C{level}{k1+2},C{level}{k1+1});
if isequal(HC{level},0)
    basis=[];
    return
elseif size(HC{level},2)>1  %If multiple generators, use varargin to decide which to use
    if isempty(varargin) || size(varargin{1},2)==0
        disp(HC{level})
        error('Please enter which generator of the above you want as x or [x,y]')
    else
        GC{level}=GC{level}(:,varargin{1}(1));
    end
end
[GD{level},HD{level}]=Homology(D{level}{k2+2},D{level}{k2+1});
if isequal(HD{level},0)
    basis=[];
    return
elseif size(HD{level},2)>1
    if isempty(varargin) || size(varargin{1},2)==1
        disp(HD{level})
        error('Please enter which generator of the above you want as [0,y] or [x,y]')
    else
        GD{level}=GD{level}(:,varargin{1}(2));
    end
end

%Box at the bottom (you can only box downstairs!!)
Boxed=cell(1,4); totalrank=cell(1,4);

[totalrank{1},detailedrank{1},Boxed{1}]=Box(rankC{1},rankD{1},C{1},D{1},useData,Data); 


%Restrict our generators to the bottom, so as to be able to multiply them (only multiply downstairs!!)
if level==4
    GC{2}=restrict(GC{4},rankC{4}{k1+1},rankC{2}{k1+1});
    GD{2}=restrict(GD{4},rankD{4}{k2+1},rankD{2}{k2+1});
end
if level>=2
    GC{1}=restrict(GC{2},rankC{2}{k1+1},rankC{1}{k1+1});
    GD{1}=restrict(GD{2},rankD{2}{k2+1},rankD{1}{k2+1});
end
    
%Now multiply on the bottom. To do that we need to pad left and right by zeros

padleft=0;
for j=0:k1-1  %If the j corresponds to [] in rank, then sum=0 so we are good!
    padleft=padleft+sum(detailedrank{1}{k1+k2+1}{j+1});
end
padright=0;
for j=k1+1:k1+k2
    padright=padright+sum(detailedrank{1}{k1+k2+1}{j+1});
end

%The product on the bottom in the left convenient basis
leftconvproduct=GC{1}.*GD{1}(:)'; 
leftconvproduct=leftconvproduct(:);  %We multiply the GC{1} with each entry of GD{1} and then concatenate vertically.


%Then we change the basis
[convlefttocanon, ~]=boxchangebasisfast(rankC{1}{k1+1},rankD{1}{k2+1},useData,Data);
productcanon=convlefttocanon*leftconvproduct;
%Finally we pad it
product=cell(1,4);
product{1}=[zeros(padleft,1);productcanon;zeros(padright,1)];



%Now transfer upstairs (no need to transfer the detailedrank, that's only for padding). 
%We don't transfer the product, rather use the inverse of restriction as we want to get the product, not a multiple of it. 
for l=[2,4]
    Boxed{l}=cell(1,k1+k2+2);
    totalrank{l}=cell(1,k1+k2+2); %We only really need to transfer at k1+k2+1 and k1+k2+2, so we'll just stop there
    for i=[k1+k2+1,k1+k2+2]
        domain=totalrank{1}{i};
        if i==1
            range=0;
        else
            range=totalrank{1}{i-1};
        end
        Boxed{l}{i}=transferdifferential(Boxed{1}{i},4/l,domain,range);
        totalrank{l}{i}=ranktransfer(totalrank{l/2}{i},l/2);
    end
    product{l}=invres(product{l/2},totalrank{l/2}{k1+k2+1},l/2);
end

[Generatoratproduct,Homologyatproduct,SmithVariables]=Homology(Boxed{level}{k1+k2+2},Boxed{level}{k1+k2+1});
q=Homologyelement(product(level),SmithVariables);
basis=q{1};
%For varargout{1} we output the simplified basis, where are allowed to
%rechoose the generator so as to get the answer to be 0,1,2 for Z/4 etc.
%So it's true only up to signs
varargout{1}=abs(basis);
if Homologyatproduct==4 && basis==3
    varargout{1}=1;
end
varargout{2}{1}=HC{level};
varargout{2}{2}=HD{level};
varargout{2}{3}=Homologyatproduct;

varargout{3}{1}=GC{level};
varargout{3}{2}=GD{level};
varargout{3}{3}=Generatoratproduct;


end