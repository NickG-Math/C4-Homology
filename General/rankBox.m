function [rank,detailedrank]=rankBox(i,rankC,rankD)
%
%INPUT: int i and cells of arrays rankC,rankD
%
%OUTPUT: Arrays rank and detailedrank
%
%DESCRIPTION: Compute the rank of the box product of chain complexes given their ranks at the index i. 
%
%detailedrank stores the rank of each summand in the box product separately. It's only used when multiplying generators (see C4mult).

s=size(rankC,2)-2; %C_0 up to C_s and C_{s+1}=0
t=size(rankD,2)-2; %D_0 up to D_t and D_{t+1}=0
%The box product is (C otimes D)_0 up to (C otimes D)_{s+t} and (C otimes D)_{s+t+1}=0

if i>=0
    detailedrank=cell(1,i+1);
    for j=max(0,i-t):min(i,s) %(C otimes D)_i=C_0 otimes D_i+...+C_j otimes D_{i-j}+...+C_i otimes D_i. This range is when C_j otimes D_{i-j} actually exists.
        detailedrank{j+1}=rankmult(rankC{j+1},rankD{i-j+1});
    end
    rank=[detailedrank{:}]; %Append them all to the rank
else  %To avoid extra if/else in other functions
    detailedrank{1}=[];
    rank=[];
end

end