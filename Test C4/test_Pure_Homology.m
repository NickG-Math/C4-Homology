function Test_Pure_Homology(rangeN,rangeM,useData,Data)

for n=0:rangeN
    for m=0:rangeM
        for k=0:n+2*m
            Answer1=C4Mackey(k,n,m,useData,Data);
            Answer2=C4MackeyAnswer(k,n,m);
            if isequal(Answer1,Answer2)
                sprintf('The %d homology of the (%d,%d) sphere is %s',k,n,m,Answer1)
            else
                error('The %d homology of the (%d,%d) sphere is not %s',k,n,m,Answer2)
            end
        end
    end
end