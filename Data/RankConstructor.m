function rank=RankConstructor(maxpower,maxlength) 
%Inputs: ints maxpower, maxlength
%Outputs: cell of arrays rank
%Description: rank contains all arrays of length<=maxlength of the form [2^?,2^maxpower,...,2^maxpower] or [2^maxpower,...,2^maxpower,2^?] where ?<=maxpower
%Could be made faster by preallocating rank but whatever.
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