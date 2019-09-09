function Boxed=BoxDiff(i,rankC,rankD,C,D,useData,Data)
%
%INPUT: int i, cells of arrays rankC,rankD, cells of matrices C,D, logical useData and struct Data
%
%OUTPUT: Matrix Boxed. 
%
%DESCRIPTION: Given chain complexes C,D of ranks rankC,rankD, Boxed is the i-th differential of C tensor D.

s=size(C,2)-2; %C_0 up to C_s and C_{s+1}=0
t=size(D,2)-2; %D_0 up to D_t and D_{t+1}=0
%The box product is (C otimes D)_0 up to (C otimes D)_{s+t} and (C otimes D)_{s+t+1}=0

if i>0
    LeftDiff=cell(1,i+1); RightDiff=cell(1,i+1);
    %Initialized left and right differentials. Now we compute them bit by bit
    
    for j=max(0,i-t):min(i,s) %(C otimes D)_i=C_0 otimes D_i+...+C_j otimes D_{i-j}+...+C_i otimes D_i. This range is when C_j otimes D_{i-j} actually exists.
        
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
        if j>=1  %Otherwise the left differential is empty
            
            convLeftdiff=blkdiagopt(C{j+1},sum(rankD{i-j+1}));   %The left differential in the convenient basis   
            
            %Another one of these provisions follows
            if ~isequal(rankC{j},1) && ~isequal(rankD{i-j+1},1)
                [ConvLeftToCanonRange, ~]=boxchangebasis(rankC{j},rankD{i-j+1},useData,Data); %The change of basis matrix in the range
                if flag
                    LeftDiff{j+1}=ConvLeftToCanonRange*convLeftdiff*ConvLeftToCanonDomain'; %The left differential in the canonical basis
                else
                    LeftDiff{j+1}=ConvLeftToCanonRange*convLeftdiff; %The left differential in the canonical basis
                end
            else
                if flag
                    LeftDiff{j+1}=convLeftdiff*ConvLeftToCanonDomain'; %The left differential in the canonical basis
                else
                    LeftDiff{j+1}=convLeftdiff; %The left differential in the canonical basis
                end
            end
        end
        
        %Now we do the right differential. Same story so we won't repeat the comments
        if i-j>=1
            convRightdiff=blkdiagopt(D{i-j+1},sum(rankC{j+1}));
            if ~isequal(rankC{j+1},1) && ~isequal(rankD{i-j},1)
                [~,ConvRightToCanonRange]=boxchangebasis(rankC{j+1},rankD{i-j},useData,Data);
                if flag
                    RightDiff{j+1}=(-1)^j*ConvRightToCanonRange*convRightdiff*ConvRightToCanonDomain';
                else
                    RightDiff{j+1}=(-1)^j*ConvRightToCanonRange*convRightdiff;
                end
            else
                if flag
                    RightDiff{j+1}=(-1)^j*convRightdiff*ConvRightToCanonDomain';
                else
                    RightDiff{j+1}=(-1)^j*convRightdiff;
                end
            end
        end
    end
    
    %Left and right differentials have been computed and we can assemble them
    Boxed=matrixmixer(LeftDiff,RightDiff);
else %i==0
    Boxed=[]; %No point computing an empty differential
end

end