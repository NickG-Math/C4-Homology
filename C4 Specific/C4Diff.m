function [rankDom,rankRan,Diff]=C4Diff(i,n,m,useData,Data)
%
%INPUT: ints i,n,m, logical useData, struct Data
%
%OUTPUT: arrays rankDom, rankRan and matrix Diff
%
%DESCRIPTION: Computes the i-th differential for S^{nsigma+mlambda} and the ranks of its domain and range. The latter is required for transferring

if n*m>=0 %No need to box
        [rankDom,Diff]=C4StandardDiff(i,n,m);
        rankRan=size(Diff,1);
end

if n*m<0 %We need to box        
        rankCsigma=cell(1,min(i+1,abs(n)+1)); Csigma=rankCsigma; rankClambda=cell(1,min(i+1,2*abs(m)+1)); Clambda=rankClambda;
        for j=0:size(rankCsigma,2)-1
            [rankCsigma{j+1},Csigma{j+1}]=C4StandardDiff(j,n,0);
        end
        for j=0:size(rankClambda,2)-1
            [rankClambda{j+1},Clambda{j+1}]=C4StandardDiff(j,0,m);
        end
 
%Now deprecated due to high memory usage. If ever used again, make sure NOT to use the whole chains but to remove the last [] !!
%rankCsigma=Data.rankStandard{sign(n)+2,2,abs(n)+1,1}; %We need all (or at least most) of the Chains to box
%         Csigma=Data.ChainsStandard{sign(n)+2,2,abs(n)+1,1};
%         rankClambda=Data.rankStandard{2,sign(m)+2,1,abs(m)+1}; %We need all (or at least most) of the Chains to box
%         Clambda=Data.ChainsStandard{2,sign(m)+2,1,abs(m)+1};

%Time to box. As before, the rank of the range is needed for transfering
    
    %We box lambda first sigma second
    
    rankDom=rankBox(i,rankClambda,rankCsigma);
    rankRan=rankBox(i-1,rankClambda,rankCsigma);
    Diff=BoxDiff(i,rankClambda,rankCsigma,Clambda,Csigma,useData,Data);
    
    %We can also box the other way
    
%     rankDom{1}=rankBox(i,rankCsigma,rankClambda);
%     rankRan{1}=rankBox(i-1,rankCsigma,rankClambda);
%     Diff{1}=BoxDiff(i,rankCsigma,rankClambda,Csigma,Clambda,useData,Data);

end
end