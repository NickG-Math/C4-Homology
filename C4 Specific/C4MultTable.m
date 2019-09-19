function Table=C4MultTable(level,SerialHomology,IndexedHomology,BasicIrreducibles,useData,Data)
%Table=C4MultTable(level,SerialHomology,IndexedHomology,BasicIrreducibles,useData,Data)
%
%INPUTS: int level, cells SerialHomology,IndexedHomology,BasicIrreducibles, logical useData, struct Data
%
%OUTPUT: 4d cell Table that contains a cell of arrays
%
%Description: Creates the multiplication table of all generators in
%SerialHomology multiplied with the BasicIrreducibles (the Euler and
%orientation classes). It also accounts for multiple generators present in 
%the same degree (non cyclic homology).

for i=1:size(BasicIrreducibles,2)
    for j=1:size(SerialHomology,2)
        index=PositiveIndexer(SerialHomology{j});
        homology=IndexedHomology{index(1),index(2),index(3)};
        for k=1:size(homology,2) %If there are multiple generators
            Sumindex=PositiveIndexer(SerialHomology{j}+BasicIrreducibles{i});
            [~,Table{i,Sumindex(1),Sumindex(2),Sumindex(3)}{k}]=C4mult(level,BasicIrreducibles{i},SerialHomology{j},useData,Data,[0,k]);
        end
    end
end