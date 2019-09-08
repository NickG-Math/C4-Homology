function [A,varargout] = Smithoptimal(A,wantP,wantQ)

[M,N] = size(A);
L=min(M,N);
if wantP
    P=eye(M);
end
if wantQ
    Q=eye(N);
end

start=1;
while start<=L
    %It's actually faster NOT to check if A is 0 in the range here
    %Find the minimum element of A in absolute value
    
    minim=abs(A(start,start));
    s=start;
    t=start;
    %Using min to find the index of the minimum element is actually
    %significantly slower and uglier than double looping
    for i=start:M
        for j=start:N
            if A(i,j)~=0 && (minim==0 || abs(A(i,j))<minim)
                minim=abs(A(i,j));
                s=i;
                t=j;
            end
        end
    end
    
    if minim==0 %If A is 0 in our range then we are done
        break
    end
    
    %Bring minimum to the start
    if s~=start
        A([start,s],:)=A([s,start],:);
        if wantP
            P([start,s],:)=P([s,start],:);
        end
    end
    if t~=start
        A(:,[start,t])=A(:,[t,start]);
        if wantQ
            Q(:,[start,t])=Q(:,[t,start]);
        end
    end
    
    %See if A(start,start) divides everything, otherwise add to the first
    %column
    
    %Eliminate using A(start,start)
    
    for i=start+1:M
        if A(i,start)~=0
            thequotient=floor(A(i,start)/A(start,start));
            A(i,:)=A(i,:)-thequotient*A(start,:);
            if wantP
                P(i,:)=P(i,:)-thequotient*P(start,:);
            end
        end
    end
    for j=start+1:N
        if A(start,j)~=0
            thequotient=floor(A(start,j)/A(start,start));
            A(:,j)=A(:,j)-thequotient*A(:,start);
            if wantQ
                Q(:,j)=Q(:,j)-thequotient*Q(:,start);
            end
        end
    end
    
    %It's actually faster not switching with a remainder as soon as we see it, as there may be a smaller remainder
    %later.
    
    %If rows and columns are 0 proceed to next iteration
    if ~any(A(start+1:end,start)) && ~any(A(start,start+1:end))
        start=start+1;
    end
end

if wantP && wantQ
    varargout{1}=P;
    varargout{2}=Q;
elseif wantP
    varargout{1}=P;
elseif wantQ
    varargout{1}=Q;
end

end