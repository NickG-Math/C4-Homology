function Name=ExpressionParser(product,Counter,nameCounter)
%
%Name=ExpressionParser(product,Counter,nameCounter)
%
%INPUTS: int product, array Counter, cell of strings nameCounter
%
%OUTPUT: string Name
%
%DESCRIPTION: Give the names of the MoreIrreducibles and the product,
%Counter specified by the Factorization function, return the name of the
%generator in the usual form x/y.


num=cell(1,size(Counter,2));
den=cell(1,size(Counter,2));
for i=1:size(Counter,2)
    if Counter(i)>1
        num{i}=strcat(nameCounter{i},'^',num2str(Counter(i)),'*');
    elseif Counter(i)==1
        num{i}=strcat(nameCounter{i},'*');
    elseif Counter(i)==-1
        den{i}=strcat(nameCounter{i},'*');
    elseif Counter(i)<0
        den{i}=strcat(nameCounter{i},'^',num2str(-Counter(i)),'*');
    end
end
if product~=1
    Nominator=strcat(num2str(product),'*',num{:});
else
    Nominator=strcat(num{:});
end
if ~isempty(Nominator)
    Nominator(end)='';
end
Denominator=strcat(den{:});
if ~isempty(Denominator)
    Denominator(end)='';
    Denominator=['/',Denominator];
end

Name=strcat(Nominator,Denominator);
end