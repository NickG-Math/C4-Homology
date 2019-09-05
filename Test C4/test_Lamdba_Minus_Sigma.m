function test_Lamdba_Minus_Sigma(rangeN,rangeM,useData,Data)

% %For Even n>=2 and any m>=1 means that I only get Z,Z/4,Z/2,Q,0

for n=2:2:rangeN
    for m=1:rangeM
        for k=-n:2*m
            [~,~,~,~,~,Answer]=C4Mackey(k,-n,m,useData,Data);
            found=0;
            if k==-n+2*m
                if isequal(Answer,"Z")
                    sprintf('The %d homology of (%d,%d) sphere is Z',k,-n,m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of (%d,%d) sphere is not Z',k,-n,m)
                end
            end
            if -n+1<=k && k<=-3 && mod(k,2)==1
                if isequal(Answer,"Z/2")
                    sprintf('The %d homology of (%d,%d) sphere is Z/2',k,-n,m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of (%d,%d) sphere is not Z/2',k,-n,m)
                end
            end
            if -n+2<=k && k<-n+2*m && mod(k,2)==0
                if  isequal(Answer,"Z/4")
                    sprintf('The %d homology of (%d,%d) sphere is Z/4',k,-n,m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of (%d,%d) sphere is not Z/4',k,-n,m)
                end
            end
            if  k==-n
                if isequal(Answer,"Q")
                    sprintf('The %d homology of (%d,%d) sphere is Q',k,-n,m)
                    found=found+1;
                else
                    disp(Answer)
                    error('The %d homology of (%d,%d) sphere is not Q',k,-n,m)
                end
            end
            if isequal(Answer,"0")
                sprintf('The %d homology of (%d,%d) sphere is 0',k,-n,m)
                found=found+1;
            end
            if found~=1
                disp(Answer)
                error('The %d homology of (%d,%d) sphere is not accounted for',k,-n,m)
            end
        end
    end
end









%%%%%Odd n means Z_-, Z/2 plus Z_-, overlineZ/2, Z/2, Z/2 plus overlineZ/2, Q

for n=1:2:rangeN
    for m=1:rangeM
        for k=-n:2*m
            [~,~,~,~,~,Answer]=C4Mackey(k,-n,m,useData,Data);
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
end

