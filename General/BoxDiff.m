function Boxed=BoxDiff(i,rankC,rankD,C,D,useData,Data)
%
%Boxed=BoxDiff(i,rankC,rankD,C,D,useData,Data)
%
%INPUT: int i, cells of arrays rankC,rankD, cells of matrices C,D, logical useData and struct Data
%
%OUTPUT: Matrix Boxed. 
%
%DESCRIPTION: Given chain complexes C,D of ranks rankC,rankD, Boxed is the i-th differential of C tensor D.

s=size(C,2)-1;  %C_0 up to C_s 
t=size(D,2)-1;  %D_0 up to D_t
%The box product is (C otimes D)_0 up to (C otimes D)_{s+t}

if i>0 %The first differential i=0 always maps to 0
    LeftDiff=cell(1,i+1);  RightDiff=cell(1,i+1);
    %Initialized left and right differentials. Now we compute them bit by bit
    
    for j=max(0,i-t):min(i,s) %(C otimes D)_i=C_0 otimes D_i+...+C_j otimes D_{i-j}+...+C_i otimes D_i. This range is when C_j otimes D_{i-j} actually exists.
        
        [Left_to_Canon_Domain,Right_to_Canon_Domain]=boxchangebasis(rankC{j+1},rankD{i-j+1},useData,Data); %Get the change of basis permutations in the domain
                
        %Now we compute the left differential C_j box D_{i-j}-> C_{j-1} box D_{i-j}
        if j>=1  %Otherwise the left differential is empty
            convLeftdiff=blkdiagopt(C{j+1},sum(rankD{i-j+1}));   %The left differential in the convenient basis
            [Left_to_Canon_Range,~]=boxchangebasis(rankC{j},rankD{i-j+1},useData,Data); %change of basis permutation in the range
            LeftDiff{j+1}=convLeftdiff(Left_to_Canon_Range,Left_to_Canon_Domain); %applying the permutations to the matrix to get the leftdiff in the canonical basis
        end
        
        %Now we do the right differential. Same story.
        if i-j>=1
            convRightdiff=blkdiagopt(D{i-j+1},sum(rankC{j+1}));
            [~,Right_to_Canon_Range]=boxchangebasis(rankC{j+1},rankD{i-j},useData,Data);
            RightDiff{j+1}=(-1)^j*convRightdiff(Right_to_Canon_Range,Right_to_Canon_Domain);
        end
    end
    %Left and right differentials have been computed and we can assemble them
    Boxed=matrixmixer(LeftDiff,RightDiff);
else %i==0
    Boxed=[];
end

end