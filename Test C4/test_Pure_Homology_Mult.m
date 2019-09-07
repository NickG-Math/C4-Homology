function test_Pure_Homology_Mult(rangeN1,rangeN2,rangeM1,rangeM2,useData,Data)


asigma=[0,1,0];
u2sigma=[2,2,0];
usigma=[1,1,0];
alambda=[0,0,1];
ulambda=[2,0,1];

%Top level
for n1=0:rangeN1
    for n2=0:rangeN2
        for m1=0:rangeM1
            for m2=0:rangeM2
                if n1<=1 || m2==0 %Otherwise we get the Gold
                    elementnow=n1*asigma+n2*u2sigma+m1*alambda+m2*ulambda;
                    if n1==0 && n2==0 && m1==0 && m2==0 
                        sprintf('Top Generator verified: 1')
                        continue
                    end
                    
                    if n1>0
                        [~,a]=C4mult(4,asigma,elementnow-asigma,useData,Data);
                    elseif n2>0
                        [~,a]=C4mult(4,u2sigma,elementnow-u2sigma,useData,Data);
                    elseif m1>0
                        [~,a]=C4mult(4,alambda,elementnow-alambda,useData,Data);
                    elseif m2>0
                        [~,a]=C4mult(4,ulambda,elementnow-ulambda,useData,Data);
                    end
                    if a==1
                        sprintf('Top Generator verified: asigma^%d*u2sigma^%d*alambda^%d*ulambda^%d',n1,n2,m1,m2)
                    else
                        error('Top Generator not verified')
                    end
                end
            end
            
        end
    end
end


%Mid level
for n=0:rangeN1
    for m1=0:rangeM1
        for m2=0:rangeM2
            elementnow=n*usigma+m1*alambda+m2*ulambda;
            if n==0 && m1==0 && m2==0 
                sprintf('Mid Generator verified: 1')
                continue
            end
            if n>0
                [~,a]=C4mult(2,usigma,elementnow-usigma,useData,Data);
            elseif m1>0
                [~,a]=C4mult(2,alambda,elementnow-alambda,useData,Data);
            elseif m2>0
                [~,a]=C4mult(2,ulambda,elementnow-ulambda,useData,Data);
            end
            if a==1
                sprintf('Mid Generator verified: usigma^%d*bar(alambda)^%d*bar(ulambda)^%d',n,m1,m2)
            else
                error('Mid Generator not verified')
            end
        end
    end
end


% Gold verification
[~,a]=C4mult(4,2*asigma,ulambda,useData,Data);
if a==2
    sprintf('Gold Verified!')
else
    error('Gold not verified!')
end


