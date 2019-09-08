function transfer=transferdifferential(A,n,dom,ran)  
%Inputs are matrix A, arrays dom,ran and int n
%Outputs: Matrix transfer
%Description: Takes the differential A at bottom level (level=1) and transfers it up n levels. 
%dom and ran are the ranks of the domain and range of the differential and they are needed to transfer correctly for box products
%This function does not follow the usual rules of preallocation+vectorization
%However I had to do extra calculations for the preallocation, so it was actually slower!
%Anyway this code could probably be improved.

if isempty(A)
    transfer=A; 
    return
end    

M=size(A,1);

%We need to discard rows n+1..., up to ran(1), then n+ran(1),... up to ran(2) and so on (this is math). So we need to track what we discard
trackvert=1;
discard=cell(1,size(ran,2));
for i=1:size(ran,2)
    limit=min(ran(i)+trackvert-1,M);
    discard{i}=n+trackvert:limit;
    trackvert=trackvert+ran(i);
end
A([discard{:}],:)=[]; %We discard the aforementioned rows.

%We now build the tranfer by adding certain columns of A indicated below
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