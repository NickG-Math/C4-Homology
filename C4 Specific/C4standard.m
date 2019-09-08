function [rank,Chains]=C4standard(n,m)
%Inputs: ints n,m
%Outputs: Cell of arrays rank and cell of matrices Chains of equal length
%Descriptions: Chains is the chain complex in the bottom level for the (n,m) sphere, with
%Chains{i} being the i-1 differential. rank collects the ranks


rank=cell(1,abs(n)+2*abs(m)+2);
Chains=cell(1,abs(n)+2*abs(m)+2);
Chains{1}=0;   %We start with Z->0 always
if n==0 && m==0
    Chains{2}=[]; %Chains{2}=0 does not make sense
    rank{1}=1;
    return
elseif n<=0 && m<=0  %For cohomology
    [~,Chainspre]=C4standard(abs(n),abs(m));
    for i=1:abs(n+2*m)+1
        Chains{i}=Chainspre{abs(n+2*m)+3-i}';
    end
    Chains{abs(n+2*m)+2}=[];
    if m==0
        Chains{1}=zeros(1,2); %Fix Chains{1}
    else
        Chains{1}=zeros(1,4);
    end
    for i=1:abs(n+2*m)+1
    rank{i}=size(Chains{i},2);
    end
    return
end

if n>0
    Chains{2}=ones(1,2);
elseif m>0
    Chains{2}=ones(1,4);
end
for i=2:2:n %Even at most n
    Chains{i+1}=altmatrix([1,-1],2,2);
end
for i=3:2:n %odd at most n
    Chains{i+1}=ones(2);
end
if mod(n,2)==0
    if n>0
        Chains{n+2}=ones(2,4); %Bridge differential for even n>0
    end
    for i=n+2:2:n+2*m
        Chains{i+1}=altmatrix([1,0,0,-1],4,4); %Higher even differentials for even n
    end
    for i=n+3:2:n+2*m
        Chains{i+1}=ones(4,4); %Higher odd differentials for even n
    end
else
    Chains{n+2}=altmatrix([1,-1],2,4); %Bridge differential for odd n
    for i=n+2:2:n+2*m
        Chains{i+1}=altmatrix([1,0,0,1],4,4);   %Higher odd differentials for odd n
    end
    for i=n+3:2:n+2*m
        Chains{i+1}=altmatrix([1,-1],4,4);  %Higher even differentials for odd n
    end
end
Chains{n+2*m+2}=[];  %EXTEREMELLY IMPORTANT
for i=1:n+2*m+2
    rank{i}=size(Chains{i},2);
end
end
