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
        Answer=C4Mackey(k,-n,m,useData,Data);
        found=0;
        if k==-n+2*m && k>=-1
            if isequal(Answer,"Z_-")
                sprintf('The %d homology of (%d,%d) sphere is Z_-',k,-n,m)
                found=found+1;
            else
                disp(Answer)
                error('The %d homology of (%d,%d) sphere is not Z_-!',k,-n,m)
            end
        end
        if (-n+1<=k && k<2*m-n &&  mod(k,2)==0) || (2*m-n<k && k<=-3 && mod(k,2)==1 && k~=2*m-n)
            if isequal(Answer,"Z/2")
                sprintf('The %d homology of (%d,%d) sphere is Z/2',k,-n,m)
                found=found+1;
            else
                disp(Answer)
                error('The %d homology of (%d,%d) sphere is not Z/2!',k,-n,m)
            end
        end
        if -1<=k &&  k<-n+2*m && mod(k,2)==1
            if isequal(Answer,"V")
                sprintf('The %d homology of (%d,%d) sphere is V',k,-n,m)
                found=found+1;
            else
                disp(Answer)
                error('The %d homology of (%d,%d) sphere is not  V!',k,-n,m)
            end
            
        end
        
        if k==-n && k<=-3
            if isequal(Answer,"Q")
                sprintf('The %d homology of (%d,%d) sphere is Q',k,-n,m)
                found=found+1;
            else
                disp(Answer)
                error('The %d homology of (%d,%d) sphere is not Q!',k,-n,m)
            end
        end
        if k>=-n+2 && k<-n+2*m && k<=-3 && mod(k,2)==1
            if isequal(Answer,"Z/2+V")
                sprintf('The %d homology of (%d,%d) sphere is Z/2 plus V',k,-n,m)
                found=found+1;
            else
                disp(Answer)
                error('The %d homology of (%d,%d) sphere is not Z/2 plus V!',k,-n,m)
            end
        end
        
        if  k==-n+2*m && k<=-3
            if isequal(Answer,"Z/2+Z_-")
                sprintf('The %d homology of (%d,%d) sphere is Z_- plus Z/2',k,-n,m)
                found=found+1;
            else
                disp(Answer)
                error('The %d homology of (%d,%d) sphere is not Z_- plus Z/2!',k,-n,m)
            end
        end
        
        if isequal(Answer,"0")
            sprintf('The %d homology of (%d,%d) sphere is 0',k,-n,m)
            found=found+1;
        end
        if found~=1
            disp(Answer)
            error('The %d homology of (%d,%d) sphere has not been accounted for',k,-n,m)
        end
    end
end 

