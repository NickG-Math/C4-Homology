function rank=rankmult(a,b) %Finds the rank of (Z_{C_{a1}}+...+Z_{C_{an}}) box (Z_{C_{b1}}+...+Z_{C_{bm}})
if size(a,2)==1 && size(b,2)==1
    if a>b
        rank=a*ones(1,b);
    else
        rank=b*ones(1,a);
    end
else
    rank=[];
    for i=1:size(a,2)
        for j=1:size(b,2)
            rank=[rank,rankmult(a(i),b(j))];
        end
    end
end


% The nonvectorized version is actrually faster than
%  helpcell=cell(size(a,2),size(b,2));
%     for i=1:size(b,2)
%         for j=1:size(a,2)
%            helpcell{i,j}=rankmult(a(j),b(i));
%         end
%     end
%     rank=[helpcell{:,:}];