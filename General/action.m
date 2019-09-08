function gx=action(x,rank) 
%Inputs: column x and row rank
%Outputs: column gx of the same size as x
%Description: Takes the element x living in a level of a free Mackey functor with given rank and
%produces the element g*x where g is the generator of the group C2n

gx=zeros(size(x,1),1);
track=0;
for i=1:size(rank,2)
    gx(track+1:track+rank(i))=circshift(x(track+1:track+rank(i)),1);
    track=track+rank(i);
end
end