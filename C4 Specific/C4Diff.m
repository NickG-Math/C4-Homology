function [rankDom,Diff]=C4Diff(i,n,m,useData,Data)
%Inputs: ints i,n,m, logical useData, struct Data
%Outputs: Cell of arrays rankDom, Cell of matrices Diff
%Description: Computes the i-th differiantial for S^{nsigma+mlambda} and
%the rank of its domain at all levels. Eg Diff{2} is the diff. at level 2

rankDom=cell(1,4); Diff=cell(1,4);

if n*m>=0 %No need to box
    if useData %Get the chains from Data
        A=Data.rankStandard{sign(n)+2,sign(m)+2,abs(n)+1,abs(m)+1};
        B=Data.ChainsStandard{sign(n)+2,sign(m)+2,abs(n)+1,abs(m)+1};
        rankDom{1}=A{i+1};
        
        %Now we compute the rank of the range. This is needed to transfer the differential
        if i>=1 %Otherwise the rank at the range is empty.
            rankRan{1}=A{i};
        end
        Diff{1}=B{i+1};
    else %Don't get them from Data
        [A,B]=C4standard(n,m);
        rankDom{1}=A{i+1};
        %Now we compute the rank of the range. This is needed to transfer the differential
        if i>=1 %Otherwise the rank at the range is empty.
            rankRan{1}=A{i};
        end
        Diff{1}=B{i+1};
    end
end


if n*m<0 %We need to box
    if useData %Get the basic chains from Data
        rankCsigma=Data.rankStandard{sign(n)+2,2,abs(n)+1,1}; %We need all (or at least most) of the Chains to box
        Csigma=Data.ChainsStandard{sign(n)+2,2,abs(n)+1,1};
        rankClambda=Data.rankStandard{2,sign(m)+2,1,abs(m)+1}; %We need all (or at least most) of the Chains to box
        Clambda=Data.ChainsStandard{2,sign(m)+2,1,abs(m)+1};
    else %Don't get them from Data
        [rankCsigma,Csigma]=C4standard(n,0);
        [rankClambda,Clambda]=C4standard(0,m);
    end
    
    %Time to box. As before, the rank of the range is needed for transfering
    
    %We box lambda first sigma second
    
    rankDom{1}=rankBox(i,rankClambda,rankCsigma);
    rankRan{1}=rankBox(i-1,rankClambda,rankCsigma);
    Diff{1}=BoxDiff(i,rankClambda,rankCsigma,Clambda,Csigma,useData,Data);
    
    
    %We can also box the other way
    
%     rankDom{1}=rankBox(i,rankCsigma,rankClambda);
%     rankRan{1}=rankBox(i-1,rankCsigma,rankClambda);
%     Diff{1}=BoxDiff(i,rankCsigma,rankClambda,Csigma,Clambda,useData,Data);

end

%Now we transfer

for level=[2,4]
    rankDom{level}=ranktransfer(rankDom{level/2},level/2); %No point transferring the rank of the range
    if i>=1
        Diff{level}=transferdifferential(Diff{1},4/level,rankDom{1},rankRan{1});
    end
    if n*m>=0 && i==0 %The one correction to account for 0->Z->0 
        Diff{level}=0;
    end
end
end