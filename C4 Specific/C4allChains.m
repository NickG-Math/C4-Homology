function [rank,D]=C4allChains(n,m,useData,Data)
%Inputs: ints n,m, logical useData, struct Data
%Outputs: cell of cells of arrays rank and cell of cells of matrices D
%Description: rank{1},rank{2},rank{4} and D{1},D{2},D{4} are the ranks
%and differentials on each level, with 1 being bottom. 


rank=cell(1,4); D=cell(1,4);
%First compute the bottom level. This might need boxing
if useData %Use the standard chains stored in Data
    if n*m>=0
        rank{1}=Data.rankStandard{sign(n)+2,sign(m)+2,abs(n)+1,abs(m)+1}; %Remember our indexing changes
        D{1}=Data.ChainsStandard{sign(n)+2,sign(m)+2,abs(n)+1,abs(m)+1};
    else
        rankCsigma=Data.rankStandard{sign(n)+2,2,abs(n)+1,1}; %Load S^{nsigma}
        Csigma=Data.ChainsStandard{sign(n)+2,2,abs(n)+1,1};
        rankClambda=Data.rankStandard{2,sign(m)+2,1,abs(m)+1}; %Load S^{mlambda}
        Clambda=Data.ChainsStandard{2,sign(m)+2,1,abs(m)+1};
      %  [rank{1},~,D{1}]=Box(rankClambda,rankCsigma,Clambda,Csigma,useData,Data); %One way of doing it
        [rank{1},~,D{1}]=Box(rankCsigma,rankClambda,Csigma,Clambda,useData,Data); %The other way of doing it
    end
else  %Don't use the standard chains stored in Data
    if n*m>=0
        [rank{1},D{1}]=C4standard(n,m);
    else
        [rankCsigma,Csigma]=C4standard(n,0);
        [rankClambda,Clambda]=C4standard(0,m);        
        [rank{1},~,D{1}]=Box(rankClambda,rankCsigma,Clambda,Csigma,useData,Data); %One way of doing it
       % [rank{1},~,D{1}]=Box(rankCsigma,rankClambda,Csigma,Clambda,useData,Data); %The other way of doing it
    end
end

%Now transfer

for level=[2,4]
    D{level}=cell(1,abs(n)+2*abs(m)+2);
    rank{level}=cell(1,abs(n)+2*abs(m)+2);
    for i=1:abs(n)+2*abs(m)+2    
        rank{level}{i}=ranktransfer(rank{level/2}{i},level/2);
        if n*m>=0
            D{level}{1}=0; %Only useful in the extreme case that we get 0->Z->0 and that needs to transfer correctly to get H_0S^0=Z.
        end
        if i>1
            D{level}{i}=transferdifferential(D{1}{i},4/level,rank{1}{i},rank{1}{i-1});
        end
    end
end
end