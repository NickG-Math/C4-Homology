function test_Pure_Homology(rangeN,rangeM,useData,Data)


for n=0:rangeN
    for m=0:rangeM
        for k=0:n+2*m
            [~,~,~,~,~,Answer]=C4Mackey(k,n,m,useData,Data);
            found=0;
            if k==n+2*m && mod(n,2)==0
                if isequal(Answer,"Z")
                    sprintf('The %d homology of the (%d,%d) sphere is Z',k,n,m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of the (%d,%d) sphere is not Z',k,n,m)
                end
            end
            if k==n+2*m &&  mod(n,2)==1
                if isequal(Answer,"Z_-")
                    sprintf('The %d homology of the (%d,%d) sphere is Z_-',k,n,m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of the (%d,%d) sphere is not Z_-',k,n,m)
                end
            end
            if ((mod(n,2)==0 && k<n) || (mod(n,2)==1 && k<n+2*m ))&& mod(k,2)==0
                if isequal(Answer,"Z/2")
                    sprintf('The %d homology of the (%d,%d) sphere is Z/2',k,n,m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of the (%d,%d) sphere is not Z/2',k,n,m)
                end
            end
            if mod(n,2)==0 && k>=n && k<n+2*m && mod(k,2)==0
                if isequal(Answer,"Z/4")
                    sprintf('The %d homology of the (%d,%d) sphere is Z/4',k,n,m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of the (%d,%d) sphere is not Z/4',k,n,m)
                end
            end
            if mod(n,2)==1 && k>=n && k<n+2*m && mod(k,2)==1
                if isequal(Answer,"V")
                    sprintf('The %d homology of the (%d,%d) sphere is overline Z/2',k,n,m)
                    found=found+1;
                    
                else
                    disp(Answer)
                    error('The %d homology of the (%d,%d) sphere is not overline Z/2',k,n,m)
                end
            end
            if isequal(Answer,"0")
                sprintf('The %d homology of the (%d,%d) sphere is 0',k,n,m)
                found=found+1;
                
            end
            if found~=1
                disp(Answer)
                disp(found)
                error('The %d homology of the (%d,%d) sphere is not accounted for',k,n,m)
            end
            
        end
    end
end