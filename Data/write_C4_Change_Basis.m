function write_C4_Change_Basis(maxlength)


linearindexmax=(3^(maxlength+1)-3)/2; % 3 options 1,2,4 and at most maxlength places [-,-,-,-] hence geometric series
longprimelist=primes(1000); %Should be more than enough
primeindexmax=prod(longprimelist(1:maxlength).^2);


linearToMatrix=cell(1,linearindexmax); %Takes a linear index and outputs the matrix it came from. 
if primeindexmax>1
    primeToLinear=sparse(1,primeindexmax); %A sparse matrix taking a prime index and outputing the linear index.
else
    primeToLinear=zeros(1,primeindexmax); %No need to make it sparse
end
linearindex=1;
for length=1:maxlength
    tempindex=cell(1,length);
    [tempindex{:}]=ind2sub(3*ones(1,length),1:3^length);
    primelist=longprimelist(1:length);
    for i=1:3^length
        log2matrix=cellfun(@(x)x(i)-1,tempindex); %the actual [2,2,4,1] is [1,1,2,0] here
        primeindex=prod(primelist.^(log2matrix+1)); %we can't store [1,1,2,0] with primes hence we add 1
   %     if primeindex~=2 %Primeindex 2 means rank was 1; no point to store
   %     it as it always gives identity. But there is no difference in
   %     memory
            linearToMatrix{linearindex}=2.^log2matrix;
            primeToLinear(primeindex)=linearindex;
            linearindex=linearindex+1;
   %     end
    end
end

n=size(linearToMatrix,2);
ChangeBasis=cell(n); %Given a linear index, outputs the change of basis matrices.
for i=1:n
    for j=1:n
        ChangeBasis{i,j}=cell(1,2);
        [ChangeBasis{i,j}{:}]=boxchangebasis(linearToMatrix{i},linearToMatrix{j},0,[]);
    end
end

save('C4_Change_Basis.mat','longprimelist','linearToMatrix','primeToLinear','ChangeBasis');
