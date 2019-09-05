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





for n=-rangeN:0
    for m=-rangeM:0
        for k=-abs(n)-abs(2*m):0
            [~,~,~,~,~,Answer]=C4Mackey(k,n,m,useData,Data);
            found=0;
            
            if k==0 && n==0 && m==0
                if isequal(Answer,"Z")
                    sprintf('The %d homology of the (%d,%d) sphere is Z',k,n,m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of the (%d,%d) sphere is not Z',k,n,m)
                end
            end

            
            
            if k==n+2*m && n==-1 && m==0
                if isequal(Answer,"Z_-")
                    sprintf('The %d homology of the (%d,%d) sphere is Z_-',k,n,m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of the (%d,%d) sphere is not Z_-',k,n,m)
                end
            end
            if k==n+2*m && m==0 && n~=0 && mod(n,2)==0
                if isequal(Answer,"p*L")
                    sprintf('The %d homology of the (%d,%d) sphere is p^*L',k,n,m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of the (%d,%d) sphere is not p^*L',k,n,m)
                end
            end
            
            if k==n+2*m && n<=-3 && m==0 && mod(n,2)==1
                if isequal(Answer,"p*L_-")
                    sprintf('The %d homology of the (%d,%d) sphere is p^*L_-',k,n,m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of the (%d,%d) sphere is not p^*L_-',k,n,m)
                end
            end
            
            
            if k==n+2*m && m~=0 && mod(n,2)==0
                if isequal(Answer,"L")
                    sprintf('The %d homology of the (%d,%d) sphere is L',k,n,m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of the (%d,%d) sphere is not L',k,n,m)
                end
            end
            
            
            if k==n+2*m && m~=0 && mod(n,2)==1
                if isequal(Answer,"L_-")
                    sprintf('The %d homology of the (%d,%d) sphere is L_-',k,n,m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of the (%d,%d) sphere is not L_-',k,n,m)
                end
            end
            
            if n<0 && abs(k)>=3 && mod(k,2)==1 && ((m==0 && abs(k)<abs(n)) || (m~=0 && mod(n,2)==0 &&  abs(k)<=abs(n)+1 ) || (m~=0 && mod(n,2)==1 &&  abs(k)<abs(n)+abs(2*m)))
                if  isequal(Answer,"Z/2")
                    sprintf('The %d homology of the (%d,%d) sphere is Z/2',k,n,m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of the (%d,%d) sphere is not Z/2',k,n,m)
                end
            end
            
            if mod(n,2)==0 && abs(k)>=abs(n)+3 && abs(k)<abs(n+2*m) && mod(k,2)==1
                if isequal(Answer,"Z/4")
                    sprintf('The %d homology of the (%d,%d) sphere is Z/4',k,n,m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of the (%d,%d) sphere is not Z/4',k,n,m)
                end
            end
            
            
            if mod(n,2)==1 && abs(k)>=abs(n)+3 && abs(k)<abs(n+2*m) && mod(k,2)==0
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