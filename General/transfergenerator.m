function transfer=transfergenerator(A,dom,ran) %Tramsfers one level higher

if size(A,2)>1
    error('You can only transfer column vectors')
end

transfer=zeros(sum(ran),1);

trackdom=0;
trackran=0;
for i=1:size(dom,2)    
    if dom(i)==2*ran(i) 
        transfer(trackran+1:trackran+ran(i))=A(trackdom+1:trackdom+ran(i))+A(trackdom+ran(i)+1:trackdom+dom(i));
    elseif ran(i)==dom(i)
        transfer(trackran+1:trackran+ran(i))=2*A(trackdom+1:trackdom+ran(i));
    else
        error('We only transfer one level')
    end
    trackdom=trackdom+dom(i);
    trackran=trackran+ran(i);
end
end