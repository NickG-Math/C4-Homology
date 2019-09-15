function unrestricted=invres(A,ranktop,rankbottom)
%
%unrestricted=INVRES(A,ranktop,rankbottom)
%
%INPUT: Column A, arrays ranktop and rankbottom
%
%OUTPUT: Column unrestricted
%
%DESCRIPTION: Given A at rankbottom, get unrestricted at ranktop s.t. Res(unrestricted)=A
%
%This assumes that A is in the image of the restriction. 
%Also note that restrictions for free Mackey functors are injective.

unrestricted=zeros(sum(ranktop),1);
tracktop=0;
trackbottom=0;
for i=1:size(ranktop,2)
    unrestricted(tracktop+1:tracktop+ranktop(i))=A(trackbottom+1:trackbottom+ranktop(i));
    tracktop=tracktop+ranktop(i);
    trackbottom=trackbottom+rankbottom(i);
end