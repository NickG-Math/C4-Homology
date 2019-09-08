function [totalrank,rank,Boxed]=Box(rankC,rankD,C,D,useData,Data)  
%Inputs: Arrays rankC,rankD, matrices C,D, logical useData and struct Data
%Outputs: Cell of arrays totalrank, 2d cell of arrays rank and cell of matrices Boxed. 

%Description: Given chain complexes C,D of ranks rankC,rankD, Boxed is C tensor D. totalrank is the rank at every level
%rank is more detailed in that it records the rank for every summand separately (totalrank just append them all together). 
%This is useful for padding in multiplication

s=size(C,2)-2; %From C_0 to C_s and differentials C{1}:C_0->0 to C{s+1}:C_s->C_{s-1} and empty:0->C_{s+1}
t=size(D,2)-2; %From D_0 to D_t and differentials D{1}:D_0->0 to D{t+1}:D_s->D_{t-1} and empty:0->D_{t+1}

totalrank=cell(1,s+t+2); rank=cell(1,s+t+2); Boxed=cell(1,s+t+2);

%We go from 0 to s+t
for i=0:s+t %We compute the left&right differentials out (C box D)_i and the ranks
    rank{i+1}=cell(1,i+1); LeftDiff=cell(1,i+1); RightDiff=cell(1,i+1);
    
    for j=max(0,i-t):min(i,s) %(C box D)_i=C_0 box D_i + ... + C_j box D_{i-j} + ... + C_i box D_0. We work with one summand j at a time
        
        rank{i+1}{j+1}=rankmult(rankC{j+1},rankD{i-j+1}); %The rank of C_j box D_{i-j} in matrix form
        
        
        %We make special provisions if one of C_j or D_{i-j} is Z. In that case some change of basis matrices below are guaranteed to be identities, 
        %so we are saving time by not multiplying with identity matrices. Ultimately all these if/else statements can be removed (keeping only the if part), at the cost of speed.
        if ~isequal(rankC{j+1},1) && ~isequal(rankD{i-j+1},1)
            [ConvLeftToCanonDomain, ConvRightToCanonDomain]=boxchangebasis(rankC{j+1},rankD{i-j+1},useData,Data);
            flag=1;
        else
            flag=0;
        end
        %We have computed the change of basis matrices in the domains to be used for the computations below
        
        %Now we compute the left differential C_j box D_{i-j}-> C_{j-1} box D_{i-j}
        if j>=1 %Otherwise the left differential is empty
            
            convLeftdiff=blkdiagopt(C{j+1},sum(rankD{i-j+1}));              %The left differential in the convenient basis   
            
            %Another one of these provisions follows
            if ~isequal(rankC{j},1) && ~isequal(rankD{i-j+1},1)
                [ConvLeftToCanonRange, ~]=boxchangebasis(rankC{j},rankD{i-j+1},useData,Data); %The change of basis matrix in the range
                if flag
                    LeftDiff{j+1}=ConvLeftToCanonRange*convLeftdiff*ConvLeftToCanonDomain';  %The left differential in the canonical basis
                else
                    LeftDiff{j+1}=ConvLeftToCanonRange*convLeftdiff;  %The left differential in the canonical basis
                end
            else
                if flag
                    LeftDiff{j+1}=convLeftdiff*ConvLeftToCanonDomain';  %The left differential in the canonical basis
                else
                    LeftDiff{j+1}=convLeftdiff;  %The left differential in the canonical basis
                end
            end
        end
        
        %Now we compute the right differential C_j box D_{i-j}-> C_j box D_{i-j-1}
        if i-j>=1 %Otherwise the right differential is empty
            
            convRightdiff=blkdiagopt(D{i-j+1},sum(rankC{j+1}));
            
            %Again with the special provisions
            if ~isequal(rankC{j+1},1) && ~isequal(rankD{i-j},1)
                [~,ConvRightToCanonRange]=boxchangebasis(rankC{j+1},rankD{i-j},useData,Data); %The change of basis matrix in the range
                if flag
                   RightDiff{j+1}=(-1)^j*ConvRightToCanonRange*convRightdiff*ConvRightToCanonDomain'; %The right differential in the canonical basis
                else
                    RightDiff{j+1}=(-1)^j*ConvRightToCanonRange*convRightdiff; %The right differential in the canonical basis
                end
            else
                if flag
                    RightDiff{j+1}=(-1)^j*convRightdiff*ConvRightToCanonDomain'; %The right differential in the canonical basis
                else
                    RightDiff{j+1}=(-1)^j*convRightdiff; %The right differential in the canonical basis
                end
            end
        end
    end
    
    %We now combine the left and right differentials
    Boxed{i+1}=matrixmixer(LeftDiff,RightDiff);
    totalrank{i+1}=[rank{i+1}{:}]; %We get the total rank by appending
end
end