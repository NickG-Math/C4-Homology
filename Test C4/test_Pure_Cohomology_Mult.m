function test_Pure_Cohomology_Mult(rangeN1,rangeN2,rangeM1,rangeM2,useData,Data)

asigma=[0,1,0];
u2sigma=[2,2,0];
usigma=[1,1,0];
alambda=[0,0,1];
ulambda=[2,0,1];
w3=[-3,-3,0]; %Not needed due to the Frobenius Relations, which we assume.
x11=[-3,-1,-1];
s3=[-3,0,-2];


%Top level

%u2sigma and ulambda
for n2=-rangeN2:0
    for m2=-rangeM2:0
        element=n2*u2sigma+m2*ulambda;
        if n2==0 && m2==0
            continue
        end
        if n2<0
            [~,a]=C4mult(4,u2sigma,element,useData,Data);  %Multiply with what you have
        elseif m2<0
            [~,a]=C4mult(4,ulambda,element,useData,Data);
        end
        if n2==-1 && m2==0 %2/u2sigma times u2sigma is 2
            if a==2
                sprintf('Top Generator verified: 2/u2sigma^%d',-n2)
            else
                error('Top Generator not verified: 2/u2sigma^%d',-n2)
            end
        elseif n2==0 && m2==-1 %4/ulambda times ulambda is 4
            if a==4
                sprintf('Top Generator verified: 4/ulambda^%d',-m2)
            else
                error('Top Generator not verified: 4/ulambda^%d',-m2)
            end
        else
            if a==1
                if m2==0
                    sprintf('Top Generator verified: 2/u2sigma^%d',-n2)
                else
                    sprintf('Top Generator verified: 4/(u2sigma^%d*ulambda^%d)',-n2,-m2)
                end
            else
                error('Top Generator not verified')
            end
        end
    end
end


%u2sigma and alambda and ulambda and s3
for n2=-rangeN2:0
    for m1=-rangeM1:0
        for m2=-rangeM2:0
            element=n2*u2sigma+m1*alambda+m2*ulambda+s3;
            if n2<0
                [~,a]=C4mult(4,u2sigma,element,useData,Data);  
            elseif m1<0
                [~,a]=C4mult(4,alambda,element,useData,Data);
            elseif m2<0
                [~,a]=C4mult(4,ulambda,element,useData,Data);
            end
            if a==1
                sprintf('Top Generator verified: s_3/(u2sigma^%d*alambda^%d*ulambda^%d)',-n2,-m1,-m2)
            else
                error('Top Generator not verified')
            end
        end
    end
end


%asigma and u2sigma and alambda and ulambda and x11
for n1=-rangeN1:0
    for n2=-rangeN2:0
        for m1=-rangeM1:0
            for m2=-rangeM2:0
                if n1<0 && m2<0
                    continue; %No proposed generator with asigma and ulambda in denominator
                end
                element=n1*asigma+n2*u2sigma+m1*alambda+m2*ulambda+x11;
                if n1<0
                    [~,a]=C4mult(4,asigma,element,useData,Data);  
                elseif n2<0
                    [~,a]=C4mult(4,u2sigma,element,useData,Data);  
                elseif m1<0
                    [~,a]=C4mult(4,alambda,element,useData,Data);
                elseif m2<0
                    [~,a]=C4mult(4,ulambda,element,useData,Data);
                end
                if a==1
                    sprintf('Top Generator verified: x11/(asigma^%d*u2sigma^%d*alambda^%d*ulambda^%d)',-n1,-n2,-m1,-m2)
                else
                    error('Top Generator not verified')
                end
            end
        end
    end
end


% Mid Level
for n2=-rangeN2:0
    for m2=-rangeM2:0
        element=n2*usigma+m2*ulambda;
        if n2==0 && m2==0
            continue
        end
        if n2<0
            [~,a]=C4mult(2,usigma,element,useData,Data);  %Multiply with what you have
        elseif m2<0
            [~,a]=C4mult(2,ulambda,element,useData,Data);
        end
        if n2==0 && m2==-1 
            if a==2
                sprintf('Mid Generator verified: 2/bar(ulambda^%d)',-m2)
            else
                error('Mid Generator not verified: 2/bar(ulambda^%d)',-m2)
            end
        else
            if a==1
                sprintf('Mid Generator verified: 2/(usigma^%d*bar(ulambda^%d)',-n2,-m2)
            else
                error('Mid Generator not verified')
            end
        end
    end
end

for n2=-rangeN2:0
    for m1=-rangeM1:0
        for m2=-rangeM2:0
            element=n2*usigma+m1*alambda+m2*ulambda+s3;
            if n2<0
                [~,a]=C4mult(2,usigma,element,useData,Data);  
            elseif m1<0
                [~,a]=C4mult(2,alambda,element,useData,Data);
            elseif m2<0
                [~,a]=C4mult(2,ulambda,element,useData,Data);
            end
            if a==1
                sprintf('Mid Generator verified: s_3/(usigma^%d*bar(alambda^%d*ulambda^%d))',-n2,-m1,-m2)
            else
                error('Mid Generator not verified')
            end
        end
    end
end
