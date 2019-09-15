function mixed=matrixmixer(L,R)  
%
%mixed=MATRIXMIXER(L,R)  
%
%INPUT: Cells of matrices L,R
%
%OUTPUT: Matrix mixed
%
%DESCRIPTION: Mixes L,R into mixed as follows: 
%First we have L{1}, then R{1} directly below L{1}
%then L{2} directly to the right of R{1} then R{2} directly below L{2}...

emptyL=cellfun('isempty',L);
emptyR=cellfun('isempty',R);
emptyLandR=emptyL & emptyR;
L(emptyLandR)=[];
R(emptyLandR)=[];

%This clears the common empty values that are really not needed for
%anything and can cause problems


if isempty(L) 
    mixed=[];
    return
end
    
length=size(L,2);

%We track where to put the matrices L,R. This is used for preallocation.


trackvert=zeros(1,length+2);
trackhor=zeros(1,length+1);


trackvert(2)=size(L{1},1);

for i=1:length
    trackvert(i+2)=trackvert(i+1)+size(R{i},1);
    trackhor(i+1)=trackhor(i)+size(R{i},2);
end
if isempty(R{length})
    trackhor(length+1)=trackhor(length)+size(L{length},2);
end

mixed=zeros(trackvert(length+2),trackhor(length+1));
for i=1:length
    if trackvert(i+2)>=trackvert(i)+1 && trackhor(i+1)>=trackhor(i)+1
        mixed(trackvert(i)+1:trackvert(i+2),trackhor(i)+1:trackhor(i+1))=[L{i};R{i}];
    end
end
end