function A=PositiveIndexer(A)
A(A>0)=2*A(A>0);
A(A<=0)=1-2*A(A<=0);