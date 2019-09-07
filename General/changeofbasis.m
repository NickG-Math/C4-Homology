function change=changeofbasis(A,B) %Given two bases A,B in matrix form, one being a permutation of the other, get the change of basis matrix from A to B
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
change=(B'==A);
end
