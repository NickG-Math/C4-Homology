function D = blkdiagopt(A,k) 
%
%D = BLKDIAGOPT(A,k) 
%
%INPUT: Matrix A and int k
%
%OUTPUT: Matrix D
%
%DESCRIPTION: D is block diagonal with k blocks all equal to A

n=size(A,1);
m=size(A,2);
D = zeros(n*k,m*k);
for i=0:k-1
    rowstart=n*i;
    columnstart=m*i;
    D(rowstart+1:rowstart+n,columnstart+1:columnstart+m)=A;
end