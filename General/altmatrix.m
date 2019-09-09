function alt=altmatrix(A,m,n)
%Not worth to document this.
row=repmat(A,1,n/size(A,2));
alt=zeros(m,size(row,2));
for i=1:m
    alt(i,:)=circshift(row,i-1);
end