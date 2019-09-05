function [res]=restrict(A,dom,ran)  %Restricts down one level always. 

if size(A,2)>1
    error('You can only restrict column vectors')
end

res=zeros(sum(ran),1);

trackdom=0;
trackran=0;
for i=1:size(dom,2)    
    if ran(i)==2*dom(i) 
        res(trackran+1:trackran+dom(i))=A(trackdom+1:trackdom+dom(i));
        res(trackran+dom(i)+1:trackran+ran(i))=A(trackdom+1:trackdom+dom(i));
    elseif ran(i)==dom(i)
        res(trackran+1:trackran+dom(i))=A(trackdom+1:trackdom+dom(i));
    else %The range should be either equal to the domain or twice that
        error('Restriction only transfers one level')
    end
    trackdom=trackdom+dom(i);
    trackran=trackran+ran(i);
end