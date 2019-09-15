function gx=action(x,rank) 
%
%gx=ACTION(x,rank) 
%
%INPUT: column x and row rank
%
%OUTPUT: column gx of the same size as x
%
%DESCRIPTION: Takes the element x living in a level of a free Mackey
%functor with given rank and produces the element g*x where g is the
%generator of our group

gx=zeros(size(x,1),1);
track=0;
for i=1:size(rank,2)
    gx(track+1:track+rank(i))=circshift(x(track+1:track+rank(i)),1);
    track=track+rank(i);
end