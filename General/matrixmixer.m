function mixed=matrixmixer(L,R)   %Makes L{1} then R{1} directly below, then L{2} directly to the right of R{1} then R{2}...
% lengthL=size(L,2);  %L,R must be cell arrays
% lengthR=size(R,2);
% if lengthL~=lengthR
%     error('the Left and Right cells are not equal in length')
% end
emptyL=cellfun('isempty',L);
emptyR=cellfun('isempty',R);
emptyLandR=emptyL & emptyR;
L(emptyLandR)=[];
R(emptyLandR)=[];

if isempty(L) 
    mixed=[];
    return
end
    
length=size(L,2);

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