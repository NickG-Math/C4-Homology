function write_C4_Standard_Chains(rangeN,rangeM)
rankStandard=cell(3,3,rangeN+1,rangeM+1);
ChainsStandard=cell(3,3,rangeN+1,rangeM+1);
for n=-rangeN:rangeN
    for m=-rangeM:rangeM
        if n*m>=0 %If they have the same sign then we have the standard resolution
            [rankStandard{sign(n)+2,sign(m)+2,abs(n)+1,abs(m)+1},ChainsStandard{sign(n)+2,sign(m)+2,abs(n)+1,abs(m)+1}]=C4standard(n,m);
        end
    end
end
save('C4_Standard_Chains.mat');