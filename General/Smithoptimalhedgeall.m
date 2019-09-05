function [A,P,Q,Pi,Qi]= Smithoptimalhedgeall(A,wantP,wantQ,wantPi,wantQi) 
%This Smith is faster than the Smithoptimal if we are lucky to find the min
%value as soon as we start. Hence the "hedge". 
%For random matrices it's much slower, but for our matrices it's actually
%quite a bit faster!
[M,N] = size(A);
L=min(M,N);
if wantP
    P=eye(M);
else
    P=[];
end
if wantQ
    Q=eye(N);
else
    Q=[];
end 
if wantPi
    Pi=eye(M);
else
    Pi=[];
end
if wantQi
    Qi=eye(N);
else
    Qi=[];
end 





for start=1:L
    i=start+1;
    j=start+1;    
    while (i<=M || j<=N) && (any(A(start+1:end,start)) || any(A(start,start+1:end)))  %If the first row or column are not identically zero
        while i<=M && any(A(start+1:end,start)) %First we do the column
            if A(start,start)==0 %If the first term is zero, exchange it with the first nonzero
                q=find(A(start+1:end,start),1)+start;
                A([start,q],:)=A([q,start],:);
                if wantP
                    P([start,q],:)=P([q,start],:);
                end
                if wantPi
                    Pi(:,[start,q])=Pi(:,[q,start]);
                end

                j=start+1;
                i=q+1;
            else %If first term is nonzero
                if A(i,start)==0
                    i=i+1;
                else
                    if mod(A(i,start),A(start,start))==0 %hedge our bets that this will be first
                        thequotient=A(i,start)/A(start,start);
                        A(i,:)=A(i,:)-thequotient*A(start,:);
                        if wantP
                            P(i,:)=P(i,:)-thequotient*P(start,:);
                        end
                        if wantPi
                            Pi(:,start)=Pi(:,start)+thequotient*Pi(:,i);
                        end

                        i=i+1;
                    elseif mod(A(start,start),A(i,start))==0  %hedge our bets that this will be second
                        A([start,i],:)=A([i,start],:);
                        thequotient=A(i,start)/A(start,start);
                        A(i,:)=A(i,:)-thequotient*A(start,:);
                        if wantP
                            P([start,i],:)=P([i,start],:);
                            P(i,:)=P(i,:)-thequotient*P(start,:);
                        end
                        if wantPi
                            Pi(:,[start,i])=Pi(:,[i,start]);
                            Pi(:,start)=Pi(:,start)+thequotient*Pi(:,i);
                        end

                        i=i+1;
                        j=start+1;
                    else %If none divides the other we get a remainder, which is bad.
                        thequotient=floor(A(i,start)/A(start,start));
                        A(i,:)=A(i,:)-thequotient*A(start,:);
                        A([start,i],:)=A([i,start],:);
                        if wantP
                            P(i,:)=P(i,:)-thequotient*P(start,:);
                            P([start,i],:)=P([i,start],:);
                        end
                        if wantPi
                            Pi(:,start)=Pi(:,start)+thequotient*Pi(:,i);
                            Pi(:,[start,i])=Pi(:,[i,start]);
                        end

                        j=start+1;
                    end
                end
            end
        end
        %We do the row.
        while j<=N && any(A(start,start+1:end))
            if A(start,start)==0
                q=find(A(start,start+1:end),1)+start;
                A(:,[start,q])=A(:,[q,start]);
                if wantQ
                    Q(:,[start,q])=Q(:,[q,start]);
                end
                if wantQi
                    Qi([start,q],:)=Qi([q,start],:);
                end

                i=start+1;
                j=q+1;
            else
                if A(start,j)==0
                    j=j+1;
                else
                    if mod(A(start,j),A(start,start))==0 %If a(start,start) divides a(start,j) eliminate but don't swap (or will infinite loop)
                        thequotient=A(start,j)/A(start,start);
                        A(:,j)=A(:,j)-thequotient*A(:,start);
                        if wantQ
                            Q(:,j)=Q(:,j)-thequotient*Q(:,start);
                        end
                        if wantQi
                            Qi(start,:)=Qi(start,:)+thequotient*Qi(j,:);
                        end

                        j=j+1;
                    elseif   mod(A(start,start),A(start,j))==0 %If a(start,j) divides a(start,start) but not the opposite use a(start,j) instead
                        A(:,[start,j])=A(:,[j,start]);
                        thequotient=A(start,j)/A(start,start);
                        A(:,j)=A(:,j)-thequotient*A(:,start);
                        if wantQ
                            Q(:,[start,j])=Q(:,[j,start]);
                            Q(:,j)=Q(:,j)-thequotient*Q(:,start);
                        end
                        if wantQi
                            Qi([start,j],:)=Qi([j,start],:);
                            Qi(start,:)=Qi(start,:)+thequotient*Qi(j,:);
                        end

                        i=start+1;
                        j=j+1;
                     else %I.e. if a(start,j) does not divide a(start,start) and vice versa eliminate a(start,j) and swap
                        thequotient=floor(A(start,j)/A(start,start));
                        A(:,j)=A(:,j)-thequotient*A(:,start);
                        A(:,[start,j])=A(:,[j,start]);
                        if wantQ
                            Q(:,j)=Q(:,j)-thequotient*Q(:,start);
                            Q(:,[start,j])=Q(:,[j,start]);
                        end
                        if wantQi
                            Qi(start,:)=Qi(start,:)+thequotient*Qi(j,:);
                            Qi([start,j],:)=Qi([j,start],:);
                        end
                        i=start+1;
                    end
                end
            end
        end
    end
end


end