function test_Sigma_Minus_Lambda_Mult(rangeN1,rangeN2,rangeM1,rangeM2,useData,Data)

asigma=[0,1,0];
u2sigma=[2,2,0];
usigma=[1,1,0];
alambda=[0,0,1];
ulambda=[2,0,1];
s3=[-3,0,-2];


%%%u2sigma/ulambda
%%Remember we might have Z/2+L so we should choose the second generator for L.
for n2=1:rangeN2
    for m2=-rangeM2:-1
        element=n2*u2sigma+m2*ulambda;
        [~,a]=C4mult(4,ulambda,element,useData,Data,[0,2]);
        if m2==-2 || m2==-1
            if isequal(a,[0,1]) || a==2
                if m2==-1
                    sprintf('Top Generator verified: (2*u2sigma^%d)/ulambda',n2)
                else
                    sprintf('Top Generator verified: (4*u2sigma^%d)/ulambda^%d',n2,-m2)
                end
            else
                error('Top Generator not verified')
            end
        else
            if isequal(a,[0,1]) || isequal(a,[1,1]) || a==1 
%%Remember that we might get [1,1] instead of [0,1] depending on the order we box, but the two mackey functors are isomorphic and the isomorphism exchanges those generators.
                sprintf('Top Generator verified: (4*u2sigma^%d)/ulambda^%d',n2,-m2)
            else
                error('Top Generator not verified')
            end
        end
    end
end




%%% (asigma u2sgima s3)/(alambda ulambda)
for n1=[0,1]
    for n2=1:rangeN2
        for m1=-rangeM1:-1
            for m2=-rangeM2:-1
                a=1;
                element=n1*asigma+n2*u2sigma+m1*alambda+m2*ulambda+s3;
                if m1<0
                    [~,a]=C4mult(4,alambda,element,useData,Data);
                elseif m2<0
                    [~,a]=C4mult(4,ulambda,element,useData,Data);
                end
                if a==1
                    sprintf('Top Generator verified: (asigma^%d*u2sigma^%d*s_3)/(alambda^%d*ulambda^%d)',n1,n2,-m1,-m2)
                else
                    error('Top Generator not verified')
                end
            end
        end
    end
end


%%%% (asigma u2sigma)/(alambda)
% Remember we might have Z/2+L so we should choose the first generator for Z/2.
for n1=3:rangeN1
    for n2=1:rangeN2
        for m1=-rangeM1:-1
            a=1;
            element=n1*asigma+n2*u2sigma+m1*alambda;
            if m1<0
                [~,a]=C4mult(4,alambda,element,useData,Data,[0,1]);
            end
            if isequal(a,[1,0]) || a==1
                sprintf('Top Generator is verified: (asigma^%d*u2sigma^%d)/alambda^%d',n1,n2,-m1)
            else
                error('Top Generator is not verified')
            end
        end
    end
end



% Mid Level


for n2=1:rangeN2
    for m2=-rangeM2:-1
        element=n2*usigma+m2*ulambda;
        [~,a]=C4mult(2,ulambda,element,useData,Data,[0,2]);
        if m2==-1
            if a==2
                sprintf('Mid Generator verified: (2*usigma^%d)/bar(ulambda)',n2)
            else
                error('Mid Generator not verified')
            end
        else
            if a==1 
                sprintf('Mid Generator verified: (2*usigma^%d)/bar(ulambda^%d)',n2,-m2)
            else
                error('Mid Generator not verified')
            end
        end
    end
end



%%% (usgima s3)/(alambda ulambda)
for n2=1:rangeN2
    for m1=-rangeM1:-1
        for m2=-rangeM2:-1
            a=1;
            element=n2*usigma+m1*alambda+m2*ulambda+s3;
            if m1<0
                [~,a]=C4mult(2,alambda,element,useData,Data);
            elseif m2<0
                [~,a]=C4mult(2,ulambda,element,useData,Data);
            end
            if a==1
                sprintf('Mid Generator verified: (usigma^%d*bar(s_3))/bar(alambda^%d*ulambda^%d)',n2,-m1,-m2)
            else
                error('Mid Generator not verified')
            end
        end
    end
end
