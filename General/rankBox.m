function [rank,detailedrank]=rankBox(i,rankC,rankD)
%
%[rank,detailedrank]=RANKBOX(i,rankC,rankD)
%
%INPUT: int i and cells of arrays rankC,rankD
%
%OUTPUT: Arrays rank and detailedrank
%
%DESCRIPTION: Compute the i-th rank of the box product of chain complexes 
%given their ranks. detailedrank stores the rank of each summand separately.
%detailedrank is only used in the function "Multiply"
%
s=size(rankC,2)-1; %C_0 up to C_s
t=size(rankD,2)-1; %D_0 up to D_t
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