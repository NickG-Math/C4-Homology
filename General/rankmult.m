function rank=rankmult(a,b)
%
%rank=RANKMULT(a,b)
%
%INPUT: Arrays a,b
%
%OUTPUT: Array rank
%
%DESCRIPTION: Produces the rank of the tensor product of two modules given
%their ranks. Not preallocated or vectorized. Code could probably be improved.

if size(a,2)==1 && size(b,2)==1 %Quicker than going through the rest of the code
   rank=max(a,b)*ones(1,min(a,b)); %Just math
else
    rank=[];
    for i=1:size(a,2)
        for j=1:size(b,2) %rank=[rank(a(1),b(1)), rank(a(1),b(2)),...]
            rank=[rank,max(a(i),b(j))*ones(1,min(a(i),b(j)))];
        end
    end
end


% This preallocated version is actrually slower:
%  helpcell=cell(size(a,2),size(b,2));
%     for i=1:size(b,2)
%         for j=1:size(a,2)
%            helpcell{i,j}=rankmult(a(j),b(i));
%         end
%     end
%     rank=[helpcell{:,:}];