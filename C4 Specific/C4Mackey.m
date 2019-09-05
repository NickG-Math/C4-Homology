function [Homol,varargout]=C4Mackey(k,n,m,useData,Data)

Homol=cell(1,4); Gen=cell(1,4); Tr=cell(1,4); Res=cell(1,4); act=cell(1,2);
%Tr(2) is transfer 1->2, Tr(4) is 1->2, Res(2) is 2->1 and Res(4) is 4->2.The others are empty cells

k=C4kreindex(k,n,m);

if k<0 || k>abs(n)+2*abs(m)
    for level=[1,2,4]
        [Homol{level},Gen{level},Tr{level},Res{level},act{level}]=deal(0); %sets all of them to 0
    end
    varargout{1}=Gen;
    varargout{2}=Tr;
    varargout{3}=Res;
    varargout{4}=act;
    varargout{5}="0";
    return
end

[rank,D]=C4allChains(n,m,useData,Data);
SmithVariables=cell(1,4);
for level=[1,2,4]
    [Gen{level},Homol{level},SmithVariables{level}]=Homology(D{level}{k+2},D{level}{k+1});
end

Tr{2}=[]; Tr{4}=[];
Res{2}=[]; Res{4}=[];
act{1}=[]; act{2}=[];

for level=[1,2,4]
    starttransfer=level/2; %We only use it for level=2,4
    startrestriction=2*level; %We only use it for level=1,2
    restricted={}; transferred={}; actioned={};
    if level~=4
        restricted=cell(1,size(Gen{startrestriction},2));
        actioned=cell(1,size(Gen{level},2));
        for i=1:size(Gen{startrestriction},2)
            restricted{i}=restrict(Gen{startrestriction}(:,i),rank{startrestriction}{k+1},rank{level}{k+1});
        end
        for i=1:size(Gen{level},2)
            actioned{i}=action(Gen{level}(:,i),rank{level}{k+1});
        end
    end
    if level~=1
        transferred=cell(1,size(Gen{starttransfer},2));
        for i=1:size(Gen{starttransfer},2)
            transferred{i}=transfergenerator(Gen{starttransfer}(:,i),rank{starttransfer}{k+1},rank{level}{k+1});
        end
    end
    
    element=Homologyelement([restricted,actioned,transferred],SmithVariables{level});
    
    
    if level~=4
        for i=1:size(restricted,2)
            Res{startrestriction}=[Res{startrestriction},element{i}];
        end
        for i=size(restricted,2)+1:size(restricted,2)+size(actioned,2)
            act{level}=[act{level},element{i}];
        end
    end
    if level~=1
        for i=size(restricted,2)+size(actioned,2)+1:size(element,2)
            Tr{level}=[Tr{level},element{i}];
        end
    end
end

varargout{1}=Gen;
varargout{2}=Tr;
varargout{3}=Res;
varargout{4}=act;

%We now write the Mackey functor in our notation
H=Homol; T=Tr; R=Res;


if size(H{4},2)<=1
    for i=1:2
        if  sign(T{2*i})==-1 && sign(R{2*i})==-1 && H{i}==1 && H{2*i}==1
            T{2*i}=-T{2*i}; R{2*i}=-R{2*i};
        end
    end
    OurMackeyFunctor=[act{2},act{1},R{4},R{2},T{4},T{2},H{4},H{2},H{1}];
    for i=1:16
        if isequal(OurMackeyFunctor,Data.MackeyList{i}{1})
            varargout{5}=Data.MackeyList{i}{2};
        end
    end
else
    OurMackeyFunctor={[act{2},act{1},R{2},T{2},H{2},H{1}],H{4},T{4},R{4}};
    for i=17:18
        if isequal(OurMackeyFunctor,Data.MackeyList{i}{1})
            varargout{5}=Data.MackeyList{i}{2};
        end
    end
end


end