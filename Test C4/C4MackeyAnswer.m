function Mackey=C4MackeyAnswer(k,n,m)
%
%Mackey=C4MackeyAnswer(gen)
%
%INPUTS: ints k,n,m
%
%OUTPUTS: string Mackey
%
%DESCRIPTION: Produces the name of the Mackey for the factorization of gen as we
%specify it in our tables. Lots of case work.

if n>=0 && m>=0
    if k==n+2*m && mod(n,2)==0
        Mackey="Z";
    elseif k==n+2*m &&  mod(n,2)==1
        Mackey="Z_-";
    elseif ((mod(n,2)==0 && k<n) || (mod(n,2)==1 && k<n+2*m ))&& mod(k,2)==0 && k>=0
        Mackey="Z/2";
    elseif mod(n,2)==0 && k>=n && k<n+2*m && mod(k,2)==0
        Mackey="Z/4";
    elseif mod(n,2)==1 && k>=n && k<n+2*m && mod(k,2)==1
        Mackey="overline Z/2";
    else
        Mackey="0";
    end
elseif n<=0 && m<=0 && k<=0 && k>=n+2*m
    if k==n+2*m && n==-1 && m==0
        Mackey="Z_-";
    elseif k==n+2*m && m==0 && n~=0 && mod(n,2)==0
        Mackey="p^*L";
    elseif k==n+2*m && n<=-3 && m==0 && mod(n,2)==1
        Mackey="p^*L_-";
    elseif k==n+2*m && m~=0 && mod(n,2)==0
        Mackey="L";
    elseif k==n+2*m && m~=0 && mod(n,2)==1
        Mackey="L_-";
    elseif n<0 && k<=0 && abs(k)>=3 && mod(k,2)==1 && (m==0 && abs(k)<abs(n) || (m~=0 && mod(n,2)==0 &&  abs(k)<=abs(n)+1 ) || (m~=0 && mod(n,2)==1 &&  abs(k)<abs(n)+abs(2*m)))
        Mackey="Z/2";
    elseif mod(n,2)==0 && k<=0 && abs(k)>=abs(n)+3 && abs(k)<abs(n+2*m) && mod(k,2)==1
        Mackey="Z/4";
    elseif mod(n,2)==1&& k<=0  && abs(k)>=abs(n)+3 && abs(k)<abs(n+2*m) && mod(k,2)==0
        Mackey="overline Z/2";
    else
        Mackey="0";
    end
elseif n>0 && m<0
    m=-m;
    if mod(n,2)==1
        if  (n-2*m<k && k<=n-4 && mod(k,2)==1) || (0<=k && k<n-2*m && mod(k,2)==0)
            Mackey="Z/2";
        elseif k==n-2*m && m>=2
            Mackey="L_-";
        elseif  k==n-2 && m==1
            Mackey="Z_-^{flat}";
        elseif  k==n-3 && n>=3 && m>=2
            Mackey="Q^{sharp}";
        elseif (n-2*m<k && k<=n-5 && k<0 && mod(k,2)==0) || (k==-2 && n==1 && m>=2)
            Mackey="overline Z/2";
        elseif n-2*m<k && k<=n-5 && 0<=k && mod(k,2)==0
            Mackey="Z/2+overline Z/2";
        else
            Mackey="0";
        end
    else
        if 0<=k && k<=n-4 && mod(k,2)==0 && k~=n-2*m
            Mackey="Z/2";
        elseif n-2*m<k && k<n-3 && mod(k,2)==1
            Mackey="Z/4";
        elseif k==n-2*m && n-2*m<0 && m>=2
            Mackey="L";
        elseif k==n-2 && m==1
            Mackey="L^{sharp}";
        elseif k==n-3  && m>=2
            Mackey="Q^{sharp}";
        elseif k==n-2*m && n-2*m>=0 && m>=2
            Mackey= "L+Z/2";
        else
            Mackey="0";
        end
    end
else
    n=-n;
    if mod(n,2)==0
        if k==-n+2*m
            Mackey="Z";
        elseif -n+1<=k && k<=-3 && mod(k,2)==1
            Mackey="Z/2";
        elseif -n+2<=k && k<-n+2*m && mod(k,2)==0
            Mackey="Z/4";
        elseif  k==-n
            Mackey="Q";
        else
            Mackey="0";
        end
    else
        if k==-n+2*m && k>=-1
            Mackey="Z_-";
        elseif (-n+1<=k && k<2*m-n &&  mod(k,2)==0) || (2*m-n<k && k<=-3 && mod(k,2)==1 && k~=2*m-n)
            Mackey="Z/2";
        elseif -1<=k &&  k<-n+2*m && mod(k,2)==1
            Mackey="overline Z/2";
            
        elseif k==-n && k<=-3
            Mackey="Q";
        elseif k>=-n+2 && k<-n+2*m && k<=-3 && mod(k,2)==1
            Mackey="Z/2+overline Z/2";
            
        elseif  k==-n+2*m && k<=-3
            Mackey="Z/2+Z_-";
        else
            Mackey="0";
        end
    end
end
    
    
    