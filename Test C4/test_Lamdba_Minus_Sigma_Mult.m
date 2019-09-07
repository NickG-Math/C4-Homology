function test_Lamdba_Minus_Sigma_Mult(rangeN1,rangeN2,rangeM1,rangeM2,useData,Data)

asigma=[0,1,0];
u2sigma=[2,2,0];
usigma=[1,1,0];
alambda=[0,0,1];
ulambda=[2,0,1];
w3=[-3,-3,0];


%Top level

%%(ulambda)/ (u2sigma)
for n2=-rangeN2:-1
    for m2=1:rangeM2
        element=n2*u2sigma+m2*ulambda;
        if n2<0
            [~,a]=C4mult(4,u2sigma,element,useData,Data);
        end
        if a==1
            sprintf('Top Generator verified: ulambda^%d/u2sigma^%d',m2,-n2)
        else
            error('Top Generator is not verified')
        end
    end
end



%%%(alambda ulambda)/ (asigma u2sigma)
for n1=[-1,0]
    for n2=-rangeN2:0
        for m1=1:rangeM1
            for m2=0:rangeM2
                element=n1*asigma+n2*u2sigma+m1*alambda+m2*ulambda;
                if n1<0
                    [~,a]=C4mult(4,asigma,element,useData,Data);  
                elseif n2<0    
                    [~,a]=C4mult(4,u2sigma,element,useData,Data);  
                end    
                if (n1==-1 && n2==0 && m2==0) || (n1==-1 && m2~=0) || (n1==0 && n2==-1 && m2==0) 
                    %asigma*(2 alambda)/(asigma)=2*alambda
                    %asigma*(2 alambda ulambda)/(asigma u2sigma)=2*(alambda ulambda)/(u2sigma) 
                    %u2sigma*(2 alambda/u2sigma)=2 alambda
                    if a==2
                        if m2==0 || n1~=0
                            sprintf('Top Generator verified: (2*alambda^%d*ulambda^%d)/(asigma^%d*u2sigma^%d)',m1,m2,-n1,-n2)
                        else
                            sprintf('Top Generator verified: (alambda^%d*ulambda^%d)/u2sigma^%d',m1,m2,-n2)
                        end
                    else
                        error('Top Generator not verified')
                    end
                else
                    if a==1
                        if m2==0 || n1~=0
                            sprintf('Top Generator verified: (2 alambda^%d*ulambda^%d)/(asigma^%d*u2sigma^%d)',m1,m2,-n1,-n2)
                        else
                            sprintf('Top Generator verified: (alambda^%d*ulambda^%d)/u2sigma^%d',m1,m2,-n2)
                        end
                    else
                        error('Top Generator not verified')
                    end
                end
            end
        end
    end
end




%%(alambda w3)/ (asigma u2sigma)
for n1=-rangeN1:0
    for n2=-rangeN2:0
        for m1=1:rangeM1
            element=n1*asigma+n2*u2sigma+m1*alambda+w3;
            if n1<0
                [~,a]=C4mult(4,asigma,element,useData,Data);
            elseif n2<0
                [~,a]=C4mult(4,u2sigma,element,useData,Data);
            end
            if a==1
                sprintf('Top Generator verified: (alambda^%d*w3)/ (asigma^%d*u2sigma^%d)',m1,-n1,-n2)
            else
                error('Top Generator not verified')
            end
        end
    end
end


%Mid level

for n2=-rangeN2:0
    for m1=0:rangeM1
        for m2=0:rangeM2
            if n2==0 && m1==0 && m2==0
                continue
            end
            element=n2*usigma+m1*alambda+m2*ulambda;
            if n2<0
                [~,a]=C4mult(2,usigma,element,useData,Data);
            end
            if a==1
                    sprintf('Mid Generator verified: (alambda^%d*ulambda^%d)/usigma^%d',m1,m2,-n2)
            else
                error('Mid Generator not verified')
            end
        end
    end
end

