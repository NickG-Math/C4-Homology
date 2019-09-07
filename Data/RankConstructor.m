function rank=RankConstructor(maxpower,maxlength) %Constructs all  arrays of length<=maxlength and maximum entry maxentry that are of the form
%[2^first,2^maxpower,...,2^maxpower] or [2^maxpower,...,2^maxpower,2^last] where
%first<=maxpower and last <=maxpower. We hash them by recording first, last, and the total length.

i=1;
for j=0:maxpower-1
    rank{i}=2.^j;
    i=i+1;
end
for length=1:maxlength-1
    rank{i}=2.^repmat(maxpower,1,length);
    i=i+1;
    for extra=0:maxpower-1
        rank{i}=2.^[extra, repmat(maxpower,1,length)];
        rank{i+1}=2.^[repmat(maxpower,1,length),extra];
        i=i+2;
    end
end
rank{end+1}=2.^repmat(maxpower,1,maxlength);