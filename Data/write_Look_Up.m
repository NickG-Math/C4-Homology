function write_Look_Up(maxpower,maxlengthA,maxlengthB)
%Inputs: ints maxpower, maxlengthA, maxlengthB
%Saves variable ChangeBasis to C4_Change_Basis.mat
%Description: This is the lookup table that allows us to precompute boxchangebasis(A,B) for arrays A,B of length maxlengthA, maxlengthB resp 
%The arrays must be of the form [2^?,2^maxpower,...,2^maxpower] or [2^maxpower,...,2^maxpower,2^?] where ?<=maxpower 
%Changebsais is then a 7d cell of sparse logical matrices (or doubles for maxlengthA=maxlengthB=1)

rankA=RankConstructor(maxpower,maxlengthA); %Create all our arrays
rankB=RankConstructor(maxpower,maxlengthB);


for i=1:size(rankA,2)
    A=rankA{i}; %Grab first array
    firstA=A(1);
    lastA=A(end);
    lengthA=size(A,2);
    for j=1:size(rankB,2)
        B=rankB{j}; %Grab second array
        firstB=B(1);
        lastB=B(end);
        lengthB=size(B,2);
        if size(A,2)==1 && size(B,2)==1 %Store as doubles in this super easy case for no overhead
            [C,D]=boxchangebasis(A,B,0,[]);
            ChangeBasis{1,lengthA,lengthB,firstA,lastA,firstB,lastB}=double(C);
            ChangeBasis{2,lengthA,lengthB,firstA,lastA,firstB,lastB}=double(D);
        else %store sparse logical for memory AND speed
            [C,D]=boxchangebasis(A,B,0,[]);
            ChangeBasis{1,lengthA,lengthB,firstA,lastA,firstB,lastB}=sparse(C);
            ChangeBasis{2,lengthA,lengthB,firstA,lastA,firstB,lastB}=sparse(D);
        end
    end
end

save('C4_Change_Basis.mat','ChangeBasis');