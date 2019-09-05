clc
rangen=15;
useData=0;

%Test asigma^n*asigma^m=asigma^{n+m}
for n1=1:2*rangen
    for n2=1:2*rangen
        Answer=C4mult(4,[0,n1,0],[0,n2,0],useData,Data);
        if Answer==1
            sprintf('a_sigma^{%d}*a_sigma^{%d} generates',n1,n2)
        else
            error('n1=%d and n2=%d',n1,n2)
        end
    end
end
%Test usigma^n*usigma^m=usigma^{n+m}
for n1=-rangen:rangen
    for n2=-rangen:-rangen
        Answer=C4mult(2,[n1,n1,0],[n2,n2,0],useData,Data);
        if abs(Answer)==1
            sprintf('u_sigma^{%d}*u_sigma^{%d} generates',n1,n2)
        else
            error('n1=%d and n2=%d',n1,n2)
        end
    end
end
%Test u2sigma^n*u2sigma^m=u2sigma^{n+m}
for n1=-rangen:rangen
    for n2=-rangen:rangen
        Answer=C4mult(4,[2*n1,2*n1,0],[2*n2,2*n2,0],useData,Data);
        if abs(Answer)==1 || (n1==0 && n2==0)
            sprintf('u_2sigma^{%d}*u_2sigma^{%d} generates',n1,n2)
        elseif abs(Answer)==2 && n1<0 && n2>=0
            sprintf('2u_2sigma^{%d}*u_2sigma^{%d} is twice the generator',n1,n2)
        elseif abs(Answer)==2 && n1>=0 && n2<0
            sprintf('u_2sigma^{%d}*2u_2sigma^{%d} is twice the generator',n1,n2)
        elseif abs(Answer)==2 && n1<0 && n2<0
            sprintf('2*u_2sigma^{%d}*2u_2sigma^{%d} is twice the generator',n1,n2)
        else
            error('n1=%d and n2=%d',n1,n2)
        end
    end
end


%Test alambda^n*alambda^m=alambda^{n+m}
for n1=1:rangen
    for n2=1:rangen
        Answer=C4mult(4,[0,0,n1],[0,0,n2],useData,Data);
        if Answer==1
            sprintf('a_lambda^{%d}*a_lambda^{%d} generates',n1,n2)
        else
            error('n1=%d and n2=%d',n1,n2)
        end
    end
end
%Test ulambda^n*ulambda^m=ulambda^{n+m}
for n1=-rangen:rangen
    for n2=-rangen:rangen
        Answer=C4mult(4,[2*n1,0,n1],[2*n2,0,n2],useData,Data);
        if abs(Answer)==1 || (n1==0 && n2==0)
            sprintf('u_lambda^{%d}*u_lambda^{%d} generates',n1,n2)
        elseif abs(Answer)==4 && n1<0 && n2>=0
            sprintf('4u_lambda^{%d}*u_lambda^{%d} is 4 times the generator',n1,n2)
        elseif abs(Answer)==4 && n1>=0 && n2<0
            sprintf('u_lambda^{%d}*4u_lambda^{%d} is 4 times the generator',n1,n2)
        elseif abs(Answer)==4 && n1<0 && n2<0
            sprintf('4*u_lambda^{%d}*4u_lambda^{%d} is 4 times the generator',n1,n2)
        else
            error('n1=%d and n2=%d',n1,n2)
        end
    end
end





%Verify the gold
Answer1=C4mult(4,[0,2,0],[2,0,1],useData,Data);
Answer2=C4mult(4,[2,2,0],[0,0,1],useData,Data);
if abs(Answer1)==2 && abs(Answer2)==1
    sprintf('Gold verified')
else
    sprintf('Gold false!')
end
