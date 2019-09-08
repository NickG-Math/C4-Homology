function [Generator,Homology,SmithVariables]=Homology(D1,D0)
%Inputs: Matrices D1,D0
%Outputs: A column Generator, a row Homology and cell of matrices SmithVariables
%Description: Homology is Ker(D0)/Im(D1) in "array form" eg [2,1] means Z/2+Z
%For each entry in Homology the corresponding column of Generator is the
%generator of that group
%SmithVariables is Q0i,zerovectors,P1,modoutcompletely,Homology
%To be fed in Homologyelement


    
SmithVariables=cell(1,5);

if isempty(D0) && isempty(D1)  %If both are empty then it's 0->0->0
    Generator=0;
    Homology=0;
    return
elseif isempty(D0) %If D0=0 it's easier to sometimes give it as D0=[]. Here we fix it back to 0 (of the correct dimension)
    D0=zeros(1,size(D1,1));
elseif isempty(D1) %Same with D1
    D1=zeros(size(D0,2),1);
end

M=size(D1,1);
K=size(D0,1);

[S0,~,Q0,~,Q0i]=Smithoptimalhedgeall(D0,0,1,0,1); 
SmithVariables{1}=Q0i;

%We compute the kernel of S0
kernelS0=zeros(M); 
for i=1:M
    if i>K || S0(i,i)==0         %Kernel of diagonal matrix computation
        kernelS0(i,i)=1;
    end
end

if ~any(any(kernelS0)) %If the kernel is zero don't even bother
    Generator=zeros(M,1);
    Homology=0;
    return
end

Kernel=Q0*kernelS0;%i.e. Ker(A)=Q*Ker(S),  %The columns of Kernel are a basis for Ker(D0) plus some zero vectors. 
%Let's call the nonzero ones x,y,z,... Our generators will be in terms of these. So Kernel=(x,y,z,...) with zeros too
D1=Q0i*D1;

%D1 now in the basis above, including the zero vectors. Its columns impose relations on x,y,z,...


zerovectors=~any(Kernel); %This tells us which columns of the Kernel are 0 and can thus be discarded
SmithVariables{2}=zerovectors;


Kernel(:,zerovectors)=[]; %Discard the zero vectors. Only the x,y,z,... are left
D1(zerovectors,:)=[]; %Discard the values corresponding to zero vectors. Now imposes relations on x,y,z,...   

[S1,P1,~,P1i,~]=Smithoptimalhedgeall(D1,1,0,1,0); %P1i changes from a new basis we denote by u,v,w,... to x,y,z,...
SmithVariables{3}=P1;

Kernel=Kernel*P1i;
%The columns of the Kernel are a new basis u,v,w,... for Ker(D0). So Kernel=(u,v,w,...)

if size(S1,1)==1 || size(S1,2)==1
    Image=abs(S1(1,1)); %diag(A) is overloaded: if A is a row or column it does NOT return the first element of A but rather a square matrix with that diagonal
else
    Image=abs(diag(S1))'; %The absolute value of the diagonal of S1 which tells us which u,v,w,.. to mod out. 
end
modoutcompletely=(Image==1); %This records the diagonal values that are pm 1 i.e. which u,v,w to completely kill.

SmithVariables{4}=modoutcompletely;


Generator=Kernel; %Initialization
Generator(:,modoutcompletely)=[]; %The Generator becomes the kernel mod the image

Homology=Image;%Initialization. 
if size(Image,2)<size(Kernel,2)
    Homology=[Homology,zeros(1,size(Kernel,2)-size(Image,2))]; %Homology must have as many elements as there are columns in the kernel so as to match with the generator after the killing below
end
Homology(modoutcompletely)=[]; %We don't record the Z/Z in homology. 
Homology(Homology==0)=1; %Set the Z/0 into 1. The Z/n are set to n by default. 

if isempty(Generator)  %If everything has been killed, return 0 not empty.
    Generator=zeros(M,1);
    Homology=0;
end

SmithVariables{5}=Homology;