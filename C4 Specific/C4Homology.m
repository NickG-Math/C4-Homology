function [Homol,Gen,rank,SmithVariables]=C4Homology(levels,k,n,m,useData,Data)
%
%INPUT: array levels, ints k,n,m, logical useData and struct Data
%
%OUTPUT: cell of arrays Homol, cell of matrices Gen, cell of arrays rank, cell of cells of matrices SmithVariables
%
%DESCRIPTION: Computes the homology and generators at the requested levels=[1,2,4] at (k,n,m) for REINDEXED k 
%Also prints the rank at level and the SmithVariables (useful for other calculations)


%Preallocate for all possible levels
Homol=cell(1,4); Gen=cell(1,4); rank=cell(1,4); SmithVariables=cell(1,4); D0=cell(1,4); D1=cell(1,4);

%Do bottom level
[rank{1},rankRan,D0{1}]=C4Diff(k,n,m,useData,Data); %The rank at k, at k-1 and the differential exiting k
[rankDom,~,D1{1}]=C4Diff(k+1,n,m,useData,Data); %The rank at k+1 and the differential exiting k

%We transfer the differentials and the rank at k (no point transferring the other two ranks)
for level=levels
    if level==1 %No need to transfer
        [Gen{level},Homol{level},SmithVariables{level}]=Homology(D1{level},D0{level}); %Compute the homology
    else
        D0{level}=transferdifferential(D0{1},4/level,rank{1},rankRan);
        D1{level}=transferdifferential(D1{1},4/level,rankDom,rank{1});
        [Gen{level},Homol{level},SmithVariables{level}]=Homology(D1{level},D0{level}); %Compute the homology
        rank{level}=ranktransfer(rank{1},level,4); % Transfer the rank to level
    end
end

end