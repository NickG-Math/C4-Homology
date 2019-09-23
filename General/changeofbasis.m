function A_to_B=changeofbasis(A,B) 
%
%change=CHANGEOFBASIS(A,B) 
%
%INPUT: Matrices A,B
%
%OUTPUT: Array A_to_B
%
%DESCRIPTION: Given two bases A,B in matrix form, 
%CHANGEOFBASIS turns them into arrays (the actual bases) and 
%produces the permutation defining the change of basis matrix from A to B. 
%
%This only works if one matrix is a permutation of the other.

%We first convert the bases into arrays
if size(A,1)>1
    A=A';
    A=A(:);
    A=A';
end
if size(B,1)>1
    B=B';
    B=B(:);
    B=B';
end
%We put 1 at position (i,j) when B(j)=A(i)
A_to_B=zeros(1,size(A,2));
for i=1:size(A_to_B,2)
    A_to_B(i)=find(B(i)==A,1,'first');
end
end
%If the permutation matrix is desired instead, change=double(B'==A); This is slightly faster
%Comparison of the two methods:
%change*U=U(A_to_B,:)
%U*change'=U(:,A_to_B)
%change1*U*change2'=U(A_to_B_1,A_to_B_2)