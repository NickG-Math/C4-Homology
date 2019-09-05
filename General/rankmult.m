function rank=rankmult(a,b) %Finds the rank of (Z_{C_{a1}}+...+Z_{C_{an}}) box (Z_{C_{b1}}+...+Z_{C_{bm}})
if size(a,2)==1 && size(b,2)==1
    if a>b
        rank=a*ones(1,b);
    else
        rank=b*ones(1,a);
    end
else
    rank=bsxfun(@rankmult,a,b); %bsxfun will perform rankmult (a(1),b) then a(2),b etc. and then append them together.
end
end