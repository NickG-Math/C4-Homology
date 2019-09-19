function test_Parallel_Code_Distributed(rangeN,rangeM,useData,Data)


[rangeN,rangeM,distributedInputs]=distribute_Inputs_To_Workers(ceil((rangeN+1)/2),rangeM);
sprintf('Adjusted range to %d,%d for (n,m)',2*rangeN-1,rangeM)
distributedInputs(1,:)=2*distributedInputs(1,:)-1; %Get n to be odds from 1 and above

total=rangeN*rangeM;

opts =parforOptions(gcp,'RangePartitionMethod','fixed','SubrangeSize',total/12);

x=distributedInputs(1,:);
y=distributedInputs(2,:);
parfor (linear=1:total,opts)
    %for linear=1:total
    n=x(linear); %All odds from 0 up to the 2*range-1
    m=y(linear); %All from 0 up to the range
    for k=-n:2*m
        Answer1=C4Mackey(k,-n,m,useData,Data);
        Answer2=C4MackeyAnswer(k,-n,m);
        if k==-n+2*m && k>=-1
            if isequal(Answer1,Answer2)
                sprintf('The %d homology of (%d,%d) sphere is %s',k,-n,m,Answer1)
            else
                error('The %d homology of (%d,%d) sphere is not %s!',k,-n,Answer2)
            end
        end
    end
end

