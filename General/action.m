function gx=action(x,rank) %Takes the element x living in a group with rank and produces the element gx
gx=zeros(size(x,1),1);
track=0;
for i=1:size(rank,2)
    gx(track+1:track+rank(i))=circshift(x(track+1:track+rank(i)),1);
    track=track+rank(i);
end
end