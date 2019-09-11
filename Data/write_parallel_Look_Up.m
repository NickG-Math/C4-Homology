function write_parallel_Look_Up(maxpower,maxlengthA,maxlengthB)
%
%Parallelized version of write_Look_Up
%Not optimized by distributing arrays

rankA=RankConstructor(maxpower,maxlengthA); %Create all our arrays
rankB=RankConstructor(maxpower,maxlengthB);

ChangeBasis=cell(2,maxlengthA,maxlengthB,maxpower,maxpower,maxpower,maxpower);

V=cell(1,size(rankB,2));
U=cell(1,size(rankB,2));
parfor j=1:size(rankB,2)
    Q1=cell(1,size(rankA,2));
    Q2=cell(1,size(rankA,2));
    B=rankB{j};
    for i=1:size(rankA,2)
        A=rankA{i}; %Grab second array
        if isequal(A,1) %No point storing identities
            continue
        end
        if size(A,2)==1 && size(B,2)==1 %Store as doubles in this super easy case for no overhead
            [C,D]=boxchangebasis(A,B,0,[]);
            Q1{i}=C;
            Q2{i}=D;
        else %store sparse for memory AND speed
            [C,D]=boxchangebasis(A,B,0,[]);
            Q1{i}=sparse(logical(C));
            Q2{i}=sparse(logical(D));
        end
    end
    V{j}=Q1;
    U{j}=Q2;
end


for i=1:size(rankA,2)
    A=rankA{i}; %Grab first array
    if isequal(A,1) %No point storing identities
        continue
    end
    firstA=A(1);
    lastA=A(end);
    lengthA=size(A,2);
    for j=1:size(rankB,2)
        B=rankB{j}; %Grab second array
        firstB=B(1);
        lastB=B(end);
        lengthB=size(B,2);
        ChangeBasis{1,lengthA,lengthB,firstA,lastA,firstB,lastB}=V{j}{i};
        ChangeBasis{2,lengthA,lengthB,firstA,lastA,firstB,lastB}=U{j}{i};
    end
end


save('C4_Change_Basis.mat','ChangeBasis','-v7.3');