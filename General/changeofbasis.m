function change=changeofbasis(A,B) 
%Inputs: Matrices A,B
%Outputs: Logical matrix change
%Description: Given two bases A,B in matrix form turn then into arrays (the actual bases) and get the change of basis matrix from A to B
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
change=(B'==A); %Very fast vectorized version
end
