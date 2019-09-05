function D = blkdiagopt(A,k) %Makes block diagonal with A repeated k times
n=size(A,1);
m=size(A,2);
D = zeros(n*k,m*k);
for i=0:k-1
    rowstart=n*i;
    columnstart=m*i;
    D(rowstart+1:rowstart+n,columnstart+1:columnstart+m)=A;
end