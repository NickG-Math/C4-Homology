function test_Factorization(rangeN,rangeM,IndexedHomology,Table)

asigma=[0,1,0];
u2sigma=[2,2,0];
alambda=[0,0,1];
ulambda=[2,0,1];
w3=[-3,-3,0];
x11=[-3,-1,-1];
s3=[-3,0,-2];


BasicIrreducibles={asigma,u2sigma,alambda,ulambda};
MoreIrreducibles={asigma,u2sigma,alambda,ulambda,w3,x11,s3};
Visited=zeros(2*(rangeN+2*rangeM)+1,2*rangeN+1,2*rangeM+1,2,'logical');
nameCounter={'asigma','u2sigma','alambda','ulambda','w3','x11','s3'};

%Positives

safety=1;

SwitchLimit=0; NumOrDen=1;
for n1=0:rangeN
    for n2=0:rangeN
        for m1=0:rangeM
            for m2=0:rangeM
                gen=n1*asigma+n2*u2sigma+m1*alambda+m2*ulambda;
                if (n1==0 && n2==0 && m1==0 && m2==0) || (n1+2*n2>=rangeN) ||  (m1+m2>=rangeM)
                    continue
                end
                [product,CounterNum,CounterDen,MoreIrreducibles]=Factorization(gen,SwitchLimit,NumOrDen,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,Visited,1);
                Counter=CounterNum-CounterDen;
                [~,CounterResult]=C4FactorAnswer(gen);
                if product~=1
                    error("1")
                end
                if  isequal(Counter,CounterResult{1})
                    Name=ExpressionParser(product,Counter,nameCounter);
                    sprintf("Generator verified: %s",Name)
                else
                    error("2");
                end
            end
        end
    end
end
%
%
% %
% %
% %
% % %
% % %
% % %
% % % Top level
% % %
%
SwitchLimit=1; NumOrDen=-1;


%%u2sigma and ulambda
for n2=-rangeN:0
    for m2=-rangeM:0
        gen=n2*u2sigma+m2*ulambda;
        if (n2==0 && m2==0)  || (2*n2<=-rangeN) ||  (m2<=-rangeM)
            continue
        end
        [product,CounterNom,CounterDen]=Factorization(gen,SwitchLimit,NumOrDen,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,Visited,1);
        Counter=CounterNom-CounterDen;
        [productResult,CounterResult]=C4FactorAnswer(gen);
        if product==productResult{1} && isequal(Counter,CounterResult{1})
            Name=ExpressionParser(product,Counter,nameCounter);
            sprintf("Generator verified: %s",Name)
            
        else
            error("2")
        end
    end
end
%
%
% %%%u2sigma and alambda and ulambda and s3
for n2=-rangeN:0
    for m1=-rangeM:0
        for m2=-rangeM:0
            if (2*n2<=-rangeN) ||  (m1+m2-2<=-rangeM)
                continue
            end
            gen=n2*u2sigma+m1*alambda+m2*ulambda+s3;
            [product,CounterNom,CounterDen]=Factorization(gen,SwitchLimit,NumOrDen,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,Visited,1);
            Counter=CounterNom-CounterDen;
            [productResult,CounterResult]=C4FactorAnswer(gen);
            if product==productResult{1} && isequal(Counter,CounterResult{1})
                Name=ExpressionParser(product,Counter,nameCounter);
                sprintf("Generator verified: %s",Name)
                
            else
                error("2")
            end
        end
    end
end
%
%%%asigma and u2sigma and alambda and ulambda and x11
for n1=-rangeN:0
    for n2=-rangeN:0
        for m1=-rangeM:0
            for m2=-rangeM:0
                if n1<0 && m2<0
                    continue; %No proposed generator with asigma and ulambda in denominator
                end
                if (n1+2*n2-2<=-rangeN) ||  (m1+m2-1<=-rangeM)
                    continue
                end
                gen=n1*asigma+n2*u2sigma+m1*alambda+m2*ulambda+x11;
                [product,CounterNom,CounterDen,~]=Factorization(gen,SwitchLimit,NumOrDen,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,Visited,1);
                Counter=CounterNom-CounterDen;
                [productResult,CounterResult]=C4FactorAnswer(gen);
                %NOTE: x11/alambda=(2s3)/asigma so there are two ways this
                %can happen!
                if (product==productResult{1} && isequal(Counter,CounterResult{1})) || (product==productResult{2} && isequal(Counter,CounterResult{2}) || product==productResult{3} && isequal(Counter,CounterResult{3}))
                    Name=ExpressionParser(product,Counter,nameCounter);
                    sprintf("Generator verified: %s",Name)
                    
                else
                    error("2")
                end
            end
        end
    end
end
%
%
%
% %
% %
%
%
%%%%%(ulambda)/ (u2sigma)
for n2=-rangeN:-1
    for m2=1:rangeM
        if (2*n2<=-rangeN) || (m2>=rangeM)
            continue
        end
        gen=n2*u2sigma+m2*ulambda;
        [product,CounterNom,CounterDen]=Factorization(gen,SwitchLimit,NumOrDen,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,Visited,1);
        Counter=CounterNom-CounterDen;
        [productResult,CounterResult]=C4FactorAnswer(gen);
        if product==productResult{1} && isequal(Counter,CounterResult{1})
            Name=ExpressionParser(product,Counter,nameCounter);
            sprintf("Generator verified: %s",Name)
            
        else
            error("2")
        end
        
    end
end

% %%(alambda ulambda)/ (asigma u2sigma)
for n1=[-1,0]
    for n2=-rangeN:0
        for m1=1:rangeM
            for m2=0:rangeM
                if (n1+2*n2<=-rangeN) || (m1+m2>=rangeM) || (n1+2*n2==0)
                    continue
                end
                gen=n1*asigma+n2*u2sigma+m1*alambda+m2*ulambda;
                [product,CounterNom,CounterDen]=Factorization(gen,SwitchLimit,NumOrDen,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,Visited,1);
                Counter=CounterNom-CounterDen;
                [productResult,CounterResult]=C4FactorAnswer(gen);
                if product==productResult{1} && isequal(Counter,CounterResult{1})
                    Name=ExpressionParser(product,Counter,nameCounter);
                    sprintf("Generator verified: %s",Name)
                    
                else
                    error("2")
                end
                
            end
        end
    end
end
%

%
%(alambda w3)/ (asigma)
for n1=-rangeN:0
    for n2=-rangeN:0
        for m1=1:rangeM
            gen=n1*asigma+n2*u2sigma+m1*alambda+w3;
            
            if abs(gen(2))>=rangeN || abs(gen(3))>=rangeM
                continue
            end
            [product,CounterNom,CounterDen,]=Factorization(gen,SwitchLimit,NumOrDen,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,Visited,1);
            Counter=CounterNom-CounterDen;
            [productResult,CounterResult]=C4FactorAnswer(gen);
            if product==productResult{1} && isequal(Counter,CounterResult{1})
                Name=ExpressionParser(product,Counter,nameCounter);
                sprintf("Generator verified: %s",Name)
            else
                error("2")
            end
            
            
        end
    end
end
%

%
% % %
% % %
% % %
% % %
% % %
% %
%u2sigma/ulambda
for n2=1:rangeN
    for m2=-rangeM:-1
        gen=n2*u2sigma+m2*ulambda;
        if 2*n2>=rangeN || m2<=-rangeM
            continue
        end
        gen_Pos=PositiveIndexer(gen);
        if size(IndexedHomology{gen_Pos(1),gen_Pos(2),gen_Pos(3)},2)==1
            [product,CounterNom,CounterDen]=Factorization(gen,SwitchLimit,NumOrDen,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,Visited,1);
        else
            [product,CounterNom,CounterDen]=Factorization(gen,SwitchLimit,NumOrDen,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,Visited,2);
        end
        Counter=CounterNom-CounterDen;
        [productResult,CounterResult]=C4FactorAnswer(gen);
        if size(CounterResult,2)==1 && product==productResult{1} && isequal(Counter,CounterResult{1})
            Name=ExpressionParser(product,Counter,nameCounter);
            sprintf("Generator verified: %s",Name)
        elseif size(CounterResult,2)>1 && product==productResult{2} && isequal(Counter,CounterResult{2})
            Name=ExpressionParser(product,Counter,nameCounter);
            sprintf("Generator verified: %s",Name)
        else
            error('3')
        end
        
    end
end


%
% %
%
%
% %
% %
% %
%
%(asigma u2sigma)/(alambda)
for n1=3:rangeN
    for n2=1:rangeN
        for m1=-rangeM:-1
            if (n1+2*n2>=rangeN) || m2<=-rangeM
                continue
            end
            gen=n1*asigma+n2*u2sigma+m1*alambda;
            [product,CounterNom,CounterDen]=Factorization(gen,SwitchLimit,NumOrDen,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,Visited,1); %Multiple Gens
            Counter=CounterNom-CounterDen;
            [productResult,CounterResult]=C4FactorAnswer(gen);
            if product==productResult{1} && isequal(Counter,CounterResult{1})
                Name=ExpressionParser(product,Counter,nameCounter);
                sprintf("Generator verified: %s",Name)
                
            else
                error('3')
            end
        end
    end
end
% %
% % %
% %
% % %
safety=1; %Needed Here!
NumOrDen=1; %Much faster
%%%u2sgima (asigma  s3)/(alambda ulambda)
for n1=[0,1]
    for n2=1:rangeN
        for m1=-rangeM:-1
            for m2=-rangeM:-1
                gen=n1*asigma+n2*u2sigma+m1*alambda+m2*ulambda+s3;
                if abs(gen(2))>=rangeN || abs(gen(3))>=rangeM
                    continue
                end
                [product,CounterNom,CounterDen,~]=Factorization(gen,SwitchLimit,NumOrDen,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,Visited,1);
                Counter=CounterNom-CounterDen;
                [productResult,CounterResult]=C4FactorAnswer(gen);
                if size(CounterResult,2)==1 && product==productResult{1} && isequal(Counter,CounterResult{1})
                    Name=ExpressionParser(product,Counter,nameCounter);
                    sprintf("Generator verified: %s",Name)
                    
                elseif size(CounterResult,2)>1 && product==productResult{2} && isequal(Counter,CounterResult{2})
                    disp('Non cyclic:')
                    Name=ExpressionParser(product,Counter,nameCounter);
                    sprintf("Generator verified: %s",Name)
                    
                else
                    error('3')
                end
                
                
                
            end
        end
    end
end