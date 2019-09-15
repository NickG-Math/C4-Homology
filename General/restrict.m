function restricted=restrict(A,dom,ran) 
%
%restricted=RESTRICT(A,dom,ran) 
%
%INPUT: Column A, arrays dom and ran
%
%OUTPUT: Column res
%
%DESCRIPTION: Restricts A at rank dom to restricted at rank ran. 
%We don't need to specify the levels, only the ranks.

restricted=zeros(sum(ran),1);

trackdom=0;
trackran=0;
for i=1:size(dom,2)    
       restricted(trackran+1:trackran+ran(i))=repmat(A(trackdom+1:trackdom+dom(i)),1,ran(i)/dom(i));
       trackdom=trackdom+dom(i);
       trackran=trackran+ran(i);
end
end