function [rank,diff]=C4StandardDiff(i,n,m)
%INPUT: ints i,n,m
%
%OUTPUT: array rank and matrix Diff
%
%DESCRIPTION: Diff is the i-th differential in the chain complex for the (n,m) sphere in thebottom level. 
%
%rank is the rank of the domain


%Annoying case work but that's how it is!
if n>0 || m>0 %One is positive
    if ~mod(i,2) %even i
        if i==0
            diff=[];
            rank=1; 
            return
        elseif i<=n
            diff=[1,-1;-1,1]; %x->x-gx
            rank=2; 
            return
        elseif i==n+1 && m %nonzero m
            diff=[1,-1,1,-1;-1,1,-1,1]; %x->x-gx
            rank=4; 
            return
        elseif i<=n+2*m && ~mod(n,2) %n even
            diff=[1,0,0,-1;-1,1,0,0;0,-1,1,0;0,0,-1,1]; %x->x-gx
            rank=4; 
            return
        elseif i<=n+2*m %n odd
            diff=[1,-1,1,-1;-1,1,-1,1;1,-1,1,-1;-1,1,-1,1]; %x->x-gx+g^2x-g^3x
            rank=4; 
            return
        else
            diff=[];
            rank=[]; 
            return
        end
    elseif i==1 
        if n %nonzero n
            diff=ones(1,2); %x->1
            rank=2; 
            return
        elseif m %n=0 but m~=0
            diff=ones(1,4);
            rank=4; 
            return
        else %n=m=0
            diff=[];
            rank=[]; 
            return
        end
    else % i is odd >=3
        if i<=n
            diff=ones(2); %x->x+gx
            rank=2; 
            return
        elseif i==n+1 && m %nonzero m
            diff=ones(2,4);
            rank=4; 
            return
        elseif i<=n+2*m && ~mod(n,2) %even n
            diff=ones(4); %x->x+gx+g^2x+g^3x
            rank=4; 
            return
        elseif i<=n+2*m %odd n
            diff=[1,0,0,1;1,1,0,0;0,1,1,0;0,0,1,1]; %x->x+gx
            rank=4; 
            return        
        else
            diff=[];
            rank=[]; 
            return
        end
    end
elseif n==0 && m==0
    if i==0
        diff=0;
        rank=1;
        return
    else
        diff=[];
        rank=[];
        return
    end
else %One negative i.e. cohomology
    [~,diffpre]=C4StandardDiff(-n-2*m+1-i,-n,-m);
    diff=diffpre'; %Dualize chains to get cochains
    if i==0
        if m %m is nonzero
            rank=4;
        else
            rank=2;
        end
    else
        rank=size(diff,2);
        if rank==0
            rank=[];
        end
    end
end
end