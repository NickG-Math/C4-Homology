function write_Look_Up(maxpower,maxlengthA,maxlengthB)

rankA=RankConstructor(maxpower,maxlengthA);
rankB=RankConstructor(maxpower,maxlengthB);


for i=1:size(rankA,2)
    A=rankA{i};
    firstA=A(1);
    lastA=A(end);
    lengthA=size(A,2);
    for j=1:size(rankB,2)
        B=rankB{j};
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