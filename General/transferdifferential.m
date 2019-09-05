function transfer=transferdifferential(A,n,dom,ran)  %We transfer up n levels (1 is the bottom level).

if isempty(A)
    transfer=A; 
    return
end    

[M,~]=size(A);


%We need to discard rows n+1..., up to r(1), then n+r(1),... up to r(2). So we track them
trackvert=1;
discard=cell(1,size(ran,2));
for i=1:size(ran,2)
    limit=min(ran(i)+trackvert-1,M);
    discard{i}=n+trackvert:limit;
    trackvert=trackvert+ran(i);
end
A([discard{:}],:)=[];

transfer=[];
trackhor=1;
for i=1:size(dom,2)
    limit=dom(i)+trackhor-1;
    for j=trackhor:min(trackhor+n-1,limit-n)
        transfer=[transfer,sum(A(:,j:n:limit),2)];   %Add the j, j+n,..., in the range d(1)+...+d(i-1) to d(1)+...+d(i) adjoin the result to what has already been done
    end
    if limit-n<trackhor
        transfer=[transfer,A(:,trackhor:limit)]; %If n was too large nothing was added
    end
    trackhor=trackhor+dom(i);
end
%This makes the columns of the transfer. 




end
% if isempty(varargin)
%     dom=N;   %An array of the number of elements in each equivariant basis of the domain. By default there's just one basis of N elements
%     ran=M;   %An array of the number of elements in each equivariant basis of the range. By default there's just one basis of M elements
% else
%     if isempty(varargin{1})
%         dom=N;
%     else
%         dom=varargin{1};
%     end
%     
%     if isempty(varargin{2})
%         ran=M;
%     else
%         ran=varargin{2};
%     end
% end

