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
        [ChangeBasis{1,lengthA,lengthB,firstA,lastA,firstB,lastB},ChangeBasis{2,lengthA,lengthB,firstA,lastA,firstB,lastB}]=boxchangebasis(A,B,0,[]);
    end
end

save('C4_Change_Basis.mat','ChangeBasis');