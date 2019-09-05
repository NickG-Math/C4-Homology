function [totalrank,rank,Boxed]=Box(rankC,rankD,C,D,useData,Data)  
%Forms C box D. rank will give the rank for every summand separately (useful for padding in multiplication), while totalrank will just append them.

s=size(C,2)-2; %From C_0 to C_s and differentials C{1}:C_0->0 to C{s+1}:C_s->C_{s-1} and empty:0->C_{s+1}
t=size(D,2)-2; %From D_0 to D_t and differentials D{1}:D_0->0 to D{t+1}:D_s->D_{t-1} and empty:0->D_{t+1}

totalrank=cell(1,s+t+2); rank=cell(1,s+t+2); Boxed=cell(1,s+t+2);
%At s+t+2 we will always be empty, so no point in writing Left and Right Diffs

for i=0:s+t %We compute the left&right differentials out (C box D)_i and its rank
    rank{i+1}=cell(1,i+1); LeftDiff=cell(1,i+1); RightDiff=cell(1,i+1);
    
    for j=max(0,i-t):min(i,s) %(C box D)_i=C_0 box D_i + ... + C_j box D_{i-j} + ... + C_i box D_0. We work with one summand j at a time
        
        rank{i+1}{j+1}=rankmult(rankC{j+1},rankD{i-j+1}); %The rank of C_j box D_{i-j} in matrix form
        
        if rankC{j+1}~=1 && rankD{i-j+1}~=1 %We make special provisions if one of C_j or D_{i-j} is Z. In that case some change of basis matrices below are guaranteed to be identities, 
                                            %so we are saving time by not multiplying with them. Ultimately all these if statements can be removed, at the cost of some speed.
            [ConvLeftToCanonDomain, ConvRightToCanonDomain]=boxchangebasis(rankC{j+1},rankD{i-j+1},useData,Data);
            flag=1;
        else
            flag=0;
        end
        %The change of basis matrices in the domains to be used for the computations below
        
        %Now we compute the left differential C_j box D_{i-j}-> C_{j-1} box D_{i-j}
        if j>=1 %I.e. we have a left differential
            
            
            convLeftdiff=blkdiagopt(C{j+1},sum(rankD{i-j+1}));              %The left differential in the convenient basis   
            
            if rankC{j}~=1 && rankD{i-j+1}~=1
                [ConvLeftToCanonRange, ~]=boxchangebasis(rankC{j},rankD{i-j+1},useData,Data); %The change of basis matrix in the range
                if flag
                    LeftDiff{j+1}=ConvLeftToCanonRange*convLeftdiff*ConvLeftToCanonDomain';             %The left differential in the canonical basis
                else
                    LeftDiff{j+1}=ConvLeftToCanonRange*convLeftdiff;
                end
            else
                if flag
                    LeftDiff{j+1}=convLeftdiff*ConvLeftToCanonDomain';
                else
                    LeftDiff{j+1}=convLeftdiff;
                end
            end
        end
        
        %Now we compute the right differential C_j box D_{i-j}-> C_j box D_{i-j-1}
        if i-j>=1 %I.e. we have a right differential
            
            convRightdiff=blkdiagopt(D{i-j+1},sum(rankC{j+1}));
            if rankC{j+1}~=1 && rankD{i-j}~=1
                [~,ConvRightToCanonRange]=boxchangebasis(rankC{j+1},rankD{i-j},useData,Data); %The change of basis matrix in the range
                if flag
                   RightDiff{j+1}=(-1)^j*ConvRightToCanonRange*convRightdiff*ConvRightToCanonDomain'; %The signed right differential in the canonical basis
                else
                    RightDiff{j+1}=(-1)^j*ConvRightToCanonRange*convRightdiff; %The signed right differential in the canonical basis
                end
            else
                if flag
                    RightDiff{j+1}=(-1)^j*convRightdiff*ConvRightToCanonDomain'; %The signed right differential in the canonical basis
                else
                    RightDiff{j+1}=(-1)^j*convRightdiff; %The signed right differential in the canonical basis
                end
            end
        end
        
    end
    
    %We now combine the left and right differentials
    Boxed{i+1}=matrixmixer(LeftDiff,RightDiff);
    totalrank{i+1}=[rank{i+1}{:}]; %We get the total rank by appending
end
end