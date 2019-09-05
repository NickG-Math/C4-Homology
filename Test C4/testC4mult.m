clc 
clearvars -except Data

range=3;

name="whatever"
%Generators for positive pure part
for a=0:range
    for b=0:range
        for c=0:range
            for d=0:range
                if a<=1
                    [howmany,generator]=C4genmult(4,Data,a,b,c,d);
                    %[~,~,~,~,~,name]=C4Mackeyfast(generator(1),generator(2),generator(3),Data);
                    if howmany==1 
                        strcat(sprintf('asigma^{%d}*u2sigma^{%d}*alambda^{%d}*ulambda^{%d} generates ',a,b,c,d)," ",name)
                    else
                        error('asigma^{%d}*u2sigma^{%d}*alambda^{%d}*ulambda^{%d} does not generate',a,b,c,d)
                    end
                elseif d==0
                    [howmany,generator]=C4genmult(4,Data,a,b,c);
                    %[~,~,~,~,~,name]=C4Mackeyfast(generator(1),generator(2),generator(3),Data);
                    if howmany==1 
                        strcat(sprintf('u2sigma^{%d}*alambda^{%d}*ulambda^{%d} generates',a,b,c)," ",name)
                    else
                        error('u2sigma^{%d}*alambda^{%d}*ulambda^{%d} does not generate',a,b,c)
                    end
                end
            end
        end
    end
end

%Check Gold relations
for a=2:range
    for b=0:range
        for c=0:range
            for d=1:range
                if a==2 %Only asigma^2 means gold
                    [howmany,generator]=C4genmult(4,Data,a,b,c,d);
                    %[~,~,~,~,~,name]=C4Mackeyfast(generator(1),generator(2),generator(3),Data);
                    if howmany==2 
                         strcat(sprintf('Gold Relation verfied: asigma^{%d}*u2sigma^{%d}*alambda^{%d}ulambda^{%d} is twice the generator',a,b,c,d)," ",name)
                    else
                        error('Gold Relation not verified for asigma^{%d}*u2sigma^{%d}*alambda^{%d}*ulambda^{%d}',a,b,c,d)
                    end
                else %asigma^3 will leave one asigma to combine with 2, hence 0
                    [howmany,generator]=C4genmult(4,Data,a,b,c,d);
                    %[~,~,~,~,~,name]=C4Mackeyfast(generator(1),generator(2),generator(3),Data);
                    if howmany==0 
                         strcat(sprintf('Gold Relation verfied: asigma^{%d}*u2sigma^{%d}*alambda^{%d}ulambda^{%d} is 0',a,b,c,d)," ",name)
                    else
                        error('Gold Relation not verified for asigma^{%d}*u2sigma^{%d}*alambda^{%d}*ulambda^{%d}',a,b,c,d)
                    end
                end
            end
        end
    end
end

%Check w3 and x11 are given by the usual transfer
if C4mult(2,[-3,-3,0],[3,3,0],Data)
    sprintf('w3 is given by the usual transfer')
else
    error('w3 is not given by the usual transfer')
end
if C4mult(1,[-3,-1,-1],[3,1,1],Data)
    sprintf('x11 is given by the usual transfer')
else
    error('w3 is not given by the usual transfer')
end



%Pure negative part
for a=-range:0
    for b=-range:0
        if b==0 && a~=0
            [howmany,generator]=C4genmult(4,Data,0,a);
            %[~,~,~,~,~,name]=C4Mackeyfast(generator(1),generator(2),generator(3),Data);
            if howmany==2
                 strcat(sprintf('2u2sigma^{%d} generates',a)," ",name)
            else
                error('2u2sigma^{%d} does not generate',a)
            end
        elseif b~=0
            [howmany,generator]=C4genmult(4,Data,0,a,0,b);
            %[~,~,~,~,~,name]=C4Mackeyfast(generator(1),generator(2),generator(3),Data);
            if howmany==4
                 strcat(sprintf('4u2sigma^{%d}*ulambda^{%d} generates',a,b)," ",name)
            else
                error('4u2sigma^{%d}*ulambda^{%d} does not generate',a,b)
            end
        end
    end
end


for a=-range:0
    for b=-range:0
        for c=-range:0
            [howmany,generator]=C4genmult(4,Data,0,a,b,c,0,0,1);
            %[~,~,~,~,~,name]=C4Mackeyfast(generator(1),generator(2),generator(3),Data);
            if isequal(howmany,[1,1]) || howmany==1 %the second is when we have a single s3
                 strcat(sprintf('u2sigma^{%d}*alambda^{%d}*ulambda^{%d}*s3 generates',a,b,c)," ",name)
            else
                error('u2sigma^{%d}*alambda^{%d}*ulambda^{%d}*s3 does not generate',a,b,c)
            end
        end
    end
end

for a=-range:0
    for b=-range:0
        for c=-range:0
            for d=-range:0
                if a==0 || d==0 %Don't want asigma and ulambda
                    [howmany,generator]=C4genmult(4,Data,a,b,c,d,0,1,0);
                    %[~,~,~,~,~,name]=C4Mackeyfast(generator(1),generator(2),generator(3),Data);
                    if isequal(howmany,[1,1]) || howmany==1 %the second is when we have a single x11
                         strcat(sprintf('asigma^{%d}*u2sigma^{%d}*alambda^{%d}*ulambda^{%d}*x11 generates',a,b,c,d)," ",name)
                    else
                        error('asigma^{%d}*u2sigma^{%d}*alambda^{%d}*ulambda^{%d}*x11 does not generate',a,b,c,d)
                    end
                end
            end
        end
    end
end


%%sigma<0, lambda>0

for a=-range:-1
    for b=0:range
        for c=0:range
            if c==0 && b~=0
                [howmany,generator]=C4genmult(4,Data,0,a,b);
                %[~,~,~,~,~,name]=C4Mackeyfast(generator(1),generator(2),generator(3),Data);
                if isequal(howmany,[2,1]) || howmany==1 
                     strcat(sprintf('2u2sigma^{%d}*alambda^{%d} generates',a,b)," ",name)
                else
                    error('2u2sigma^{%d}*alambda^{%d} does not generate',a,b)
                end
            elseif c~=0
                [howmany,generator]=C4genmult(4,Data,0,a,b,c);
                %[~,~,~,~,~,name]=C4Mackeyfast(generator(1),generator(2),generator(3),Data);
                if isequal(howmany,[1,1]) || howmany==1 
                     strcat(sprintf('u2sigma^{%d}*alambda^{%d}*ulambda^{%d} generates',a,b,c)," ",name)
                else
                    error('u2sigma^{%d}*alambda^{%d}*ulambda^{%d} does not generate',a,b,c)
                end
            end
        end
    end
end

for a=-range:0
    for b=-range:0
        for c=0:range
            [howmany,generator]=C4genmult(4,Data,a,b,c,0,1);
            %[~,~,~,~,~,name]=C4Mackeyfast(generator(1),generator(2),generator(3),Data);
            if isequal(howmany,[1,1]) || howmany==1  %second for a single w3
                 strcat(sprintf('asigma^{%d}*u2sigma^{%d}*alambda^{%d}*w3 generates',a,b,c)," ",name)
            else
                error('asigma^{%d}*u2sigma^{%d}*alambda^{%d}*w3 does not generate',a,b,c)
            end
        end
    end
end

for a=-range:0
    for b=1:range+1
        for c=0:range
            [howmany,generator]=C4genmult(4,Data,-1,a,b,c);
            %[~,~,~,~,~,name]=C4Mackeyfast(generator(1),generator(2),generator(3),Data);
            if isequal(howmany,[2,1])
                strcat(sprintf('u2sigma^{%d}*asigma^{-1}*(2alambda^{%d})*ulambda^{%d} generates',a,b,c)," ",name)
            else
                error('u2sigma^{%d}*asigma^{-1}*(2alambda^{%d})*ulambda^{%d} does not generate',a,b,c)
            end
        end
    end
end

%sigma>0, lambda<0

for a=0:1
    for b=1:range+1
        for c=-range:0
            for d=-range:0
                if (c~=0 || d~=0)
                    [howmany,generator]=C4genmult(4,Data,a,b,c,d,0,0,1);
                    %[~,~,~,~,~,name]=C4Mackeyfast(generator(1),generator(2),generator(3),Data);

                    if isequal(howmany,[1,1])
                        strcat(sprintf('asigma^{%d}*u2sigma^{%d}*alambda^{%d}*ulambda^{%d}*s3 generates',a,b,c,d)," ",name)
                    else
                        error('asigma^{%d}*u2sigma^{%d}*alambda^{%d}*ulambda^{%d}*s3 does not generate',a,b,c,d)
                    end
                end
            end
        end
    end
end
%
for a=1:range+1
    [howmany,generator]=C4genmult(4,Data,0,a,0,-1);
    %[~,~,~,~,~,name]=C4Mackeyfast(generator(1),generator(2),generator(3),Data);
    if isequal(howmany,[2,1])
        strcat(sprintf('2u2sigma^{%d}**ulambda^{-1} generates',a)," ",name)
    else
        error('2u2sigma^{%d}**ulambda^{-1} does not generate',a)
    end
end


for a=0:range
    for b=-range:0
        [howmany,generator]=C4genmult(4,Data,1,a,0,b,0,0,1);
        %[~,~,~,~,~,name]=C4Mackeyfast(generator(1),generator(2),generator(3),Data);
        if isequal(howmany,[1,1]) || howmany==1
            strcat(sprintf('asigma*u2sigma^{%d}*ulambda^{%d}*s3 generates',a,b)," ",name)
        else
            error('asigma*u2sigma^{%d}*ulambda^{%d}*s3 does not generate',a,b)
        end
    end
end


for a=3:2:range+3  % YOU MUST START THIS AT 3. This is in accordance to the tables! We only use odd a to avoid the noncyclic groups
    for b=0:range
        for c=-range:-1
            [howmany,generator]=C4genmult(4,Data,a,b,c);
            %[~,~,~,~,~,name]=C4Mackeyfast(generator(1),generator(2),generator(3),Data);
            if isequal(howmany,[1,1])
                strcat(sprintf('asigma^{%d}*u2sigma^{%d}*alambda^{%d} generates',a,b,c)," ",name)
            else
                error('asigma^{%d}*u2sigma^{%d}*alambda^{%d} does not generate',a,b,c)
            end
        end
    end
end

for a=4:2:range+4
    for b=0:range
        for c=-range-1:-1
            [howmany,generator]=C4genmult(4,Data,a,b,c,0,0,0,0,1);
            %[~,~,~,~,~,name]=C4Mackeyfast(generator(1),generator(2),generator(3),Data);
            if isequal(howmany,[1,1])
                strcat(sprintf('asigma^{%d}*u2sigma^{%d}*alambda^{%d} generates the first factor',a,b,c)," ",name)
            else
                error('asigma^{%d}*u2sigma^{%d}*alambda^{%d} does not generate the first factor',a,b,c)
            end
        end
    end
end


for a=1:range
    for b=-range-2:-2
        [howmany,generator]=C4genmult(4,Data,0,a,0,b,0,0,0,2);
        %[~,~,~,~,~,name]=C4Mackeyfast(generator(1),generator(2),generator(3),Data);
        if isequal(howmany,[4,1])
            strcat(sprintf('4u2sigma^{%d}*ulambda^{%d} generates the second factor',a,b)," ",name)
        else
            error('4u2sigma^{%d}*ulambda^{%d} does not generate the second factor',a,b)
        end
    end
end
