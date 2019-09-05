function linearIndex=matrixToLinear(A,Data) %Takes the matrix for rank and produces its linear index to be used in ChangeOfBasisData


A=log2(A)+1; %Converts the A=[2,2,4,1] to [1,1,2,0]+1 which is more easily stored with our prime indexing. We need to add 1 as 0 does not lead to unique primefactorization
primelist=Data.longprimelist(1:size(A,2));
primeIndex=prod(primelist.^A);
linearIndex=Data.primeToLinear(primeIndex);
% %verify
% if ~isequal(Data.linearToMatrixData{linearIndex},2.^(A-1)) %Remmeber we have replaced A by log2(A)+1
%     disp(Data.linearToMatrixData{linearIndex})
%     error('Index Conversion Error')
% end