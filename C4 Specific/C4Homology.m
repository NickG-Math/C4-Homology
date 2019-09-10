function [Homol,Gen,rank,SmithVariables]=C4Homology(level,k,n,m,useData,Data)
%
%INPUT: ints level,k,n,m, logical useData and struct Data
%
%OUTPUT: array Homol, matrix Gen, array rank, cell of matrices SmithVariables
%
%DESCRIPTION: Computes the homology and generators at level at (k,n,m) for REINDEXED k 
%Also prints the rank at level and the SmithVariables (useful for other calculations)

%We compute the ranks and differentials, all at the bottom level
[rank_Bottom,rankRan,D0_Bottom]=C4Diff(k,n,m,useData,Data); %The rank at k, at k-1 and the differential exiting k
[rankDom,~,D1_Bottom]=C4Diff(k+1,n,m,useData,Data); %The rank at k+1 and the differential exiting k

%We transfer the differentials and the rank at k (no point transferring the other two ranks)
D0_Level=transferdifferential(D0_Bottom,4/level,rank_Bottom,rankRan);
D1_Level=transferdifferential(D1_Bottom,4/level,rankDom,rank_Bottom);

[Gen,Homol,SmithVariables]=Homology(D1_Level,D0_Level); %Compute the homology
rank=ranktransfer(rank_Bottom,level,4); % Transfer the rank to level
end