function [rank,D]=C4allChains(n,m,useData,Data)

rank=cell(1,4); D=cell(1,4);

if n*m>=0        
        rank{1}=Data.rankStandard{sign(n)+2,sign(m)+2,abs(n)+1,abs(m)+1};
        D{1}=Data.ChainsStandard{sign(n)+2,sign(m)+2,abs(n)+1,abs(m)+1};
else
    rankCsigma=Data.rankStandard{sign(n)+2,2,abs(n)+1,1}; %Load S^{nsigma}
    Csigma=Data.ChainsStandard{sign(n)+2,2,abs(n)+1,1};
    rankClambda=Data.rankStandard{2,sign(m)+2,1,abs(m)+1}; %Load S^{mlambda}
    Clambda=Data.ChainsStandard{2,sign(m)+2,1,abs(m)+1};
   
  %[rank{1},~,D{1}]=Box(rankClambda,rankCsigma,Clambda,Csigma); %One way of doing it
   [rank{1},~,D{1}]=Box(rankCsigma,rankClambda,Csigma,Clambda,useData,Data); %The other way of doing it
end


for level=[2,4]
    D{level}=cell(1,abs(n)+2*abs(m)+2);
    rank{level}=cell(1,abs(n)+2*abs(m)+2);
    for i=1:abs(n)+2*abs(m)+2    
        rank{level}{i}=ranktransfer(rank{level/2}{i},level/2);
        D{level}{1}=0;
        if i>1
            D{level}{i}=transferdifferential(D{1}{i},4/level,rank{1}{i},rank{1}{i-1});
        end
    end
end
end