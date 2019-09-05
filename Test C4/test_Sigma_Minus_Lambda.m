function test_Sigma_Minus_Lambda(rangeN,rangeM,useData,Data)

for n=1:2:rangeN
    for m=1:rangeM
        for k=-2*m:n
            [~,~,~,~,~,Answer]=C4Mackey(k,n,-m,useData,Data);
            found=0;
            if  (n-2*m<k && k<=n-4 && mod(k,2)==1) || (0<=k && k<n-2*m && mod(k,2)==0)
                if isequal(Answer,"Z/2")
                    sprintf('The %d homology of (%d,%d) sphere is Z/2',k,n,-m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of (%d,%d) sphere is not Z/2',k,n,-m)
                end
            end
            if k==n-2*m && m>=2
                if isequal(Answer,"L_-") 
                    sprintf('The %d homology of (%d,%d) sphere is L_-',k,n,-m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of (%d,%d) sphere is not L_-',k,n,-m)
                end
            end
            if  k==n-2 && m==1
                if isequal(Answer,"Z^o_-")
                    sprintf('The %d homology of (%d,%d) sphere is Z_-^o',k,n,-m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of (%d,%d) sphere is not Z^{op} minus',k,n,-m)
                end
            end
            if  k==n-3 && n>=3 && m>=2
                if isequal(Answer,"Q^o")
                    sprintf('The %d homology of (%d,%d) sphere is Q^o',k,n,-m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of (%d,%d) sphere is not Q^o',k,n,-m)
                end
            end
            if (n-2*m<k && k<=n-5 && k<0 && mod(k,2)==0) || (k==-2 && n==1 && m>=2)
                if isequal(Answer,"V")
                    sprintf('The %d homology of (%d,%d) sphere is overline Z/2',k,n,-m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of (%d,%d) sphere is not  overline Z/2',k,n,-m)
                end
            end
            if n-2*m<k && k<=n-5 && 0<=k && mod(k,2)==0
                if isequal(Answer,"Z/2+V")
                    sprintf('The %d homology of (%d,%d) sphere is Z/2 plus overline Z/2',k,n,-m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of (%d,%d) sphere is not Z/2 plus overline Z/2',k,n,-m)
                end
            end
            if isequal(Answer,"0")
                sprintf('The %d homology of (%d,%d) sphere is 0',k,n,-m)
                found=found+1;
            end
            if found~=1
                disp(Answer)
                disp(found)
                error('The %d homology of (%d,%d) sphere is not accounted for',k,n,-m)
            end
        end
    end
end

    






for n=2:2:rangeN
    for m=1:rangeM
        for k=-2*m:n
            [~,~,~,~,~,Answer]=C4Mackey(k,n,-m,useData,Data);
            found=0;
            if 0<=k && k<=n-4 && mod(k,2)==0 && k~=n-2*m
                if isequal(Answer,"Z/2")
                    sprintf('The %d homology of (%d,%d) sphere is Z/2',k,n,-m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of (%d,%d) sphere is not Z/2',k,n,-m)
                end
            end
            if n-2*m<k && k<n-3 && mod(k,2)==1
                if isequal(Answer,"Z/4")
                    sprintf('The %d homology of (%d,%d) sphere is Z/4',k,n,-m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of (%d,%d) sphere is not Z/4',k,n,-m)
                end
            end
            if k==n-2*m && n-2*m<0 && m>=2
                if isequal(Answer,"L")
                    sprintf('The %d homology of (%d,%d) sphere is L',k,n,-m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of (%d,%d) sphere is not L',k,n,-m)
                end
            end
            if k==n-2 && m==1
                if isequal(Answer,"L^o")
                    sprintf('The %d homology of (%d,%d) sphere is L^o',k,n,-m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of (%d,%d) sphere is not L^o',k,n,-m)
                end
            end
            if k==n-3  && m>=2
                if isequal(Answer,"Q^o")
                    sprintf('The %d homology of (%d,%d) sphere is Q^o',k,n,-m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of (%d,%d) sphere is not Q^o',k,n,-m)
                end
            end
            if k==n-2*m && n-2*m>=0 && m>=2
                if isequal(Answer, "L+Z/2")
                    sprintf('The %d homology of (%d,%d) sphere is L+Z/2',k,n,-m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of (%d,%d) sphere is not L+Z/2',k,n,-m)
                end
            end
            if isequal(Answer,"0")
                sprintf('The %d homology of (%d,%d) sphere is 0',k,n,-m)
                found=found+1;
            end
            if found~=1
                disp(Answer)
                disp(found)
                error('The %d homology of (%d,%d) sphere is not accounted for',k,n,-m)
            end
        end
    end
end
