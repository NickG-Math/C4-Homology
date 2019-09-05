function [convlefttocanon, convrighttocanon]=boxchangebasis(rankC,rankD,useData,Data)
%Computes the canonical and left/right convenient bases for an equivariant
%basis and returns the change of basis matrices.

if isequal(rankC,1)
    convlefttocanon=eye(sum(rankD));
    convrighttocanon=eye(sum(rankD));
    return
end
if isequal(rankD,1)
    convlefttocanon=eye(sum(rankC));
    convrighttocanon=eye(sum(rankC));
    return
end


%down the cases where you are in range
if useData
   % try  Try and catch if you expect to go outsite data range, but will slow
        first=matrixToLinear(rankC,Data);
        second=matrixToLinear(rankD,Data);
        convlefttocanon=Data.ChangeBasis{first,second}{1};
        convrighttocanon=Data.ChangeBasis{first,second}{2}; %A bit faster than [convlefttocanon,convrighttocanon]=Data.ChangeOfBasisData{matrixToLinear(rankC,Data),matrixToLinear(rankD,Data)}{:};
        return
   % catch
        %  fprintf('Change of basis for %s times %s performed without using data ',mat2str(rankC),mat2str(rankD))
   % end
end




s=size(rankC,2); %The basis for C is x_1,gx_1,...,g^{rankC(1)-1}x_1,x_2,....,g^{rankC(s)-1}x_s
t=size(rankD,2); %The basis for D is y_1,gy_1,...,g^{rankD(1)-1}y_1,y_2,....,g^{rankD(t)-1}y_t
%We work with one equivariant basis at a time
canonical=cell(s,t); leftconv=cell(s,t); rightconv=cell(t,s);
for c=1:s
    for d=1:t
        rank=rankmult(rankC(c),rankD(d));
        %The canonical basis for C tensor D is
        %xy,gxgy,...,g^{rank(1)-1}xg^{rank(1)-1}y,
        %xgy,...
        %....
        %xg^{size(rank,2)-1}y,...,
        
        i=0:size(rank,2)-1;
        j=0:rank(1)-1;

        helper=2^c*3^d*5.^mod(j,rankC(c)).*7.^mod(j+i',rankD(d)); %Due to the behavior of bsxfun, we don't need 5.^repmat(mod(j,rankC(c)),1,ilimit) in the presence of the matrix given by 7^.
        helper=helper'; helper=helper(:); 
        canonical{c,d}=helper'; %Reshape into an array so as to avoid inconsistent concatenation
        %The assembly would otherwise be canonical(1,1);canonical(2,1);... so there is no problem with reshaping right now. 
        %We don't have to/shouldn't replace left and right convenient bases right now as
        %1. they are not assembled in this way, and
        %2. their dimensions are consistent in the assembly (no padding needed)

        %The left convenient basis for C tensor D is
        %xy,gxy,...,g^{rankC(1)-1}xy,
        %xgy,...
        %....
        %xg^{rankD(1)-1}y,...,g^{rankC{1}-1}xg^{rankD{1}-1}y
                
        i=0:rankD(d)-1;
        j=0:rankC(c)-1;
        leftconv{c,d}=2^c*3^d*5.^j.*7.^repmat(i',1,rankC(c)); %Due to the behavior of bsxfun, we only need one repmat; the other is then done automatically

        %The right convenient basis for C tensor D is
        %xy,xgy,...,xg^{rankD(1)-1}y,
        %gxy,...
        %....
        %g^{rankC(1)-1}xy,...,g^{rankC{strings1}-1}xg^{rankD{1}-1}y
        rightconv{c,d}=leftconv{c,d}';
    end
end

%Assemble all equivariant bases together

%Vectorized version, barely any faster than the usual with the loop below
helper=reshape(canonical,1,s*t); %Get canonical(1,1),canonical(2,1),...
totalcanon=[helper{:}]; %Make it into a matrix.

%I am not sure hwo to do this faster as I will need to go through a cell of
%matrices row by row and concatenate. Not supported natively due to
%inconsistent dimensions.
totalleft=zeros(sum(rankD),sum(rankC));  
trackhorleft=0;
trackvertleft=0;
for d=1:t
    for c=1:s        
        totalleft(trackvertleft+1:trackvertleft+rankD(d),trackhorleft+1:trackhorleft+rankC(c))=leftconv{c,d}; %It's [left(1,1),left(2,1),..., left(s,1); left(1,2),...]
        trackhorleft=trackhorleft+rankC(c);
    end
    trackvertleft=trackvertleft+rankD(d);
    trackhorleft=0;
end

totalright=zeros(sum(rankC),sum(rankD));
trackhorright=0;
trackvertright=0;

for c=1:s
    for d=1:t
        totalright(trackvertright+1:trackvertright+rankC(c),trackhorright+1:trackhorright+rankD(d))=rightconv{c,d}; %It's [right(1,1),right(1,2),..., right(1,s); right(2,1),...]
        trackhorright=trackhorright+rankD(d);
    end
    trackvertright=trackvertright+rankC(c);
    trackhorright=0;
end

convlefttocanon=changeofbasis(totalleft,totalcanon);
convrighttocanon=changeofbasis(totalright,totalcanon);