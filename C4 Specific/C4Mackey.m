function [Name,Homol,Gen,Tr,Res,act]=C4Mackey(k,n,m,useData,Data)
%
%INPUT: ints k,n,m, logical useData and struct Data
%
%OUTPUT: String Name, cell of arrays Homol,Tr,Res,act and cell of matrices Gen 
%
%DESCRIPTION: Computes the RO(C4) homology of a point at (k,n,m) as a Mackey Functor. 
%
%Tr{2} is transfer 1->2, Tr{4} is 2->4, Res{2} is 2->1 and Res{4} is 4->2 and act{i} is the action on level i for i=1,2
%
%Name is the name of the Mackey functor


Homol=cell(1,4); Gen=cell(1,4); Tr=cell(1,4); Res=cell(1,4); act=cell(1,2);
k=C4kreindex(k,n,m); %Get a consistent index no matter the signs of n,m

if k<0 || k>abs(n)+2*abs(m) %Empty chain complexes in that range
    for level=[1,2,4]
        [Homol{level},Gen{level},Tr{level},Res{level},act{level}]=deal(0); %sets all of them to 0
    end
    Name="0";
    return
end

%rank is the rank of the chains at k, we will use it for transfers/restrictions of generators.
%SmithVariables will be used in Homologyelement to write these transfers/restrictions in terms of the generators

[Homol,Gen,rank,SmithVariables]=C4Homology([1,2,4],k,n,m,useData,Data); %Requesting all levels 1,2,4

%Now we compute transfers/restrictions/actions. First preallocate

Tr{2}=[]; Tr{4}=[];
Res{2}=[]; Res{4}=[];
act{1}=[]; act{2}=[];

for level=[1,2,4]
    domtr=level/2; %The domain of the transfer. We only use it for level=2,4
    domres=2*level; %We only use it for level=1,2
    restricted={}; transferred={}; actioned={};
    if level~=4
        restricted=cell(1,size(Gen{domres},2));
        actioned=cell(1,size(Gen{level},2));
        
        for i=1:size(Gen{domres},2)
            restricted{i}=restrict(Gen{domres}(:,i),rank{domres},rank{level}); %Restrict the i-th generator
        end
        for i=1:size(Gen{level},2)
            actioned{i}=action(Gen{level}(:,i),rank{level}); %Find the action on the i-th generator
        end
    end
    if level~=1
        transferred=cell(1,size(Gen{domtr},2));
        for i=1:size(Gen{domtr},2)
            transferred{i}=transfergenerator(Gen{domtr}(:,i),rank{domtr},rank{level});  %Transfer the i-th generator
        end
    end
    
    %Get the three elements defined above as linear combinations of generators
    element=Homologyelement([restricted,actioned,transferred],SmithVariables{level}); 
    
    %Some final bookkeeping
    if level~=4
        for i=1:size(restricted,2)
            Res{domres}=[Res{domres},element{i}];
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

%Sometimes we get isomorphic Mackey functors to the ones we have in our list because say
%Res{2}=-1 and Tr{2}=-1; -1 is an automorphism so that's the same as
%Res{2}=1 and Tr{2}=1. So now we normalize our Mackey functors
H=Homol; T=Tr; R=Res;

MackeyList=Data.MackeyList; %Load the MackeyList
if size(H{4},2)<=1 %if the level is cyclic
    %Normalize
    for i=1:2
        if  sign(T{2*i})==-1 && sign(R{2*i})==-1 && H{i}==1 && H{2*i}==1
            T{2*i}=-T{2*i}; R{2*i}=-R{2*i};
        end
    end
    %Check against our list
    OurMackeyFunctor=[act{2},act{1},R{4},R{2},T{4},T{2},H{4},H{2},H{1}];
    for i=1:16
        if isequal(OurMackeyFunctor,MackeyList{i}{1})
            Name=MackeyList{i}{2};
        end
    end
else % if we have a noncyclic level we need to be slightly more careful.
    OurMackeyFunctor={[act{2},act{1},R{2},T{2},H{2},H{1}],H{4},T{4},R{4}};
    for i=17:18
        if isequal(OurMackeyFunctor,MackeyList{i}{1})
            Name=MackeyList{i}{2};
        end
    end
end


end