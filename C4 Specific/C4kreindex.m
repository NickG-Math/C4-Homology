function k=C4kreindex(k,n,m)
if n<=0 && m<=0
    k=k-n-2*m; %Reindexes so that 0<=k and k<=abs(n)+2*abs(m). This gives a consistent index in all cases
elseif n<0 && m>0
    k=k-n;
elseif n>0 && m<0
    k=k-2*m;
end
