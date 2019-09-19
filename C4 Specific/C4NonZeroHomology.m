function [SerialHomology, IndexedHomology]=C4NonZeroHomology(rangeN,rangeM,useData,Data)
SerialHomology={};
IndexedHomology{1,1,1}=1; %[k,n,m]=[0,0,0]
for m=-rangeM:rangeM
    for n=-rangeN:rangeN
        max=abs(n)+2*abs(m);
        for k=-max:max
            l=C4kreindex(k,n,m);
            if l<0 || l>max || (n==0 && m==0)
                continue
            end
            H=C4Homology(4,l,n,m,useData,Data);
            if any(H{4})
                SerialHomology{end+1}=[k,n,m];
                IndexedHomology{PositiveIndexer(k),PositiveIndexer(n),PositiveIndexer(m)}=H{4};
            end
        end
    end
end