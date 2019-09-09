function write_C4_Standard_Chains(rangeN,rangeM)
%
%INPUT: positive ints rangeN, rangeM
%
%DESCRIPTION: Writes cell of arrays rankStandard and cell of matrices ChainsStandard to C4_Standard_Chains.mat
%
%These are the ranks and Chains of S^{nsigma+mlambda} for -rangeN<=n<=rangeN and -rangeM<=m<=rangeM

%For C4 we only really need S^nsigma and S^mlambda for quick access but whatever

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