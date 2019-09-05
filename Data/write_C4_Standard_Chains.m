function write_C4_Standard_Chains(range_n,range_m)
rankStandard=cell(3,3,range_n+1,range_m+1);
ChainsStandard=cell(3,3,range_n+1,range_m+1);
for n=-range_n:range_n
    for m=-range_m:range_m
        if n*m>=0
            [rankStandard{sign(n)+2,sign(m)+2,abs(n)+1,abs(m)+1},ChainsStandard{sign(n)+2,sign(m)+2,abs(n)+1,abs(m)+1}]=C4standard(n,m);
        end
    end
end
save('Data\C4_Standard_Chains.mat');