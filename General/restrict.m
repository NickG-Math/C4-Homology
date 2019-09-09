function [res]=restrict(A,dom,ran) 
%
%INPUT: Column A, arrays dom and ran
%
%OUTPUT: Column res
%
%DESCRIPTION: Restricts element A one level down. The dom and ran are the ranks of the two levels. 

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
        error('Domain or range provided for transfer is wrong')
    end
    trackdom=trackdom+dom(i);
    trackran=trackran+ran(i);
end