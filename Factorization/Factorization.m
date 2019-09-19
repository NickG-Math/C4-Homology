function [product,CounterNum,CounterDen,MoreIrreducibles]=Factorization(gen,SwitchLimit,NumOrDen,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,varargin)
%
%[product,CounterNum,CounterDen,MoreIrreducibles]=Factorization(gen,SwitchLimit,NumOrDen,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,varargin)
%
%INPUTS: array gen, ints SwitchLimit,NumOrDen,safety, cells Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles
%
%OPTIONAL INPUTS: matrix Visited=varargin{1}, int k=varargin{2}
%
%OUTPUTS: int product,  arrays CounterNum,CounterDen and cell MoreIrreducibles
%
%DESCRIPTION: Given a generator at gen and the multiplication table Table, 
%Factorization tries to write it as x/y where x,y are products of BasicIrreducibles
%(eg Euler+Orientation classes) and MoreIrreducibles(eg BasicIrreducibles and s3,w3,x11 for C4). 
%x,y are determined by the output variables as follows:
%x=product*MoreIrreducibles(1)^CounterNum(1)*...*MoreIrreducibles(end)^CounterNum(end)
%and 
%y=MoreIrreducibles(1)^CounterDen(1)*...*MoreIrreducibles(end)^CounterDen(end)
%
%If a factorization couldn't be found, MoreIrreducibles is extended with
%the new generator
%
%The optional argument k is needed if there are multiple generators at gen 
%(non cyclic homology). The algorithm is recursive and the variables SwitchLimit, 
%NumORDen, safety, Visited determine its behavior; see FactorizationHub.

if size(varargin,2)>0
    Visited=varargin{1};
else
    Visited=zeros(size(IndexedHomology)); %Who 
end
if size(varargin,2)>1
    k=varargin{2};
else
    k=1;
end


product=1; CounterNum=zeros(1,size(MoreIrreducibles,2)); CounterDen=CounterNum; %Visited(gen_Pos(1),gen_Pos(2),gen_Pos(3))=-1;
[product,CounterNum,CounterDen,found]=FactorizationHub(gen,SwitchLimit,NumOrDen,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,product,CounterNum,CounterDen,Visited,Visited,k);
if ~found
    disp('Probably irreducible')
    MoreIrreducibles{end+1}=gen;
end