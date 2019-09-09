function [A,P,Q,Pi,Qi]= Smithoptimalhedgeall(A,wantP,wantQ,wantPi,wantQi) 
%
%INPUT: Matrix A and logical wantP,wantQ,wantPi,wantQi
%
%Outputs: Matrices A,P,Q,Pi,Qi
%
%DESCRIPTION: The outputA is the Smith normal form of the inputA. 
%
%If all want? variables are 1, we have outputA=P*inputA*Q and inputA=Pi*outputA*Qi and P,Q,Pi,Qi are invertible matrices with Pi=inv(P) and Qi=inv(Q) 
%
%If a want? variable is 0 then the corresponding ? is [] and should NOT be used in any calculations.

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




%start indicates which row or column we are eliminating
for start=1:L
    i=start+1; %i goes through the elements on the column after start and only increases once an element is eliminated to 0. So if i>M then the whole column has been eliminated
    j=start+1; %j is the same but for rows   
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

                j=start+1;%Exchanging columns means we mess up our rows, so we need to set j back to the start. Ideally this won't happen much
                i=q+1;
            else %If first term is nonzero
                if A(i,start)==0
                    i=i+1;
                else
                    if mod(A(i,start),A(start,start))==0 %Eliminate
                        thequotient=A(i,start)/A(start,start);
                        A(i,:)=A(i,:)-thequotient*A(start,:);
                        if wantP
                            P(i,:)=P(i,:)-thequotient*P(start,:);
                        end
                        if wantPi
                            Pi(:,start)=Pi(:,start)+thequotient*Pi(:,i);
                        end

                        i=i+1;
                    elseif mod(A(start,start),A(i,start))==0  %Swap then eliminate
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
                    else %If none divides the other we get a remainder after eliminating, so we need to swap again.
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
        %The same logic applies to the working the row
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