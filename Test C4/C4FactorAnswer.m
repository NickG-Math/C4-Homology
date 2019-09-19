function [product,Counter]=C4FactorAnswer(gen)
%
%[product,Counter]=C4FactorResult(gen)
%
%INPUTS: Array gen
%
%OUTPUTS: Cells product, Counter;
%
%DESCRIPTION: produces the product and Counter for the factorization of gen as we 
%specify it in our tables. Lots of case work. Cells are needed due to the
%existence of multiple equivalent ways to write a generator

product{1}=1;
if all(gen>=0)
    k=gen(1); n=gen(2); m=gen(3);
    if mod(n,2)==0
        if k==n+2*m
            Counter{1}=[0,n/2,0,m,0,0,0];
        elseif n<=k && mod(k,2)==0
            i=(k-n)/2;
            Counter{1}=[0,n/2,m-i,i,0,0,0];
        elseif mod(k,2)==0
            i=k/2;
            Counter{1}=[n-2*i,i,m,0,0,0,0];
        end
    else
        if mod(k,2)==0 && k<=n
            i=k/2;
            Counter{1}=[n-2*i,i,m,0,0,0,0];
        elseif mod(k,2)==0
            i=(k-n+1)/2;
            Counter{1}=[1,(n-1)/2,m-i,i,0,0,0];
        end
    end
end

if all(gen<=0)
    k=gen(1); n=-gen(2); m=-gen(3);
    if mod(n,2)==0
        if k==-n-2*m && m~=0
            Counter{1}=[0,-n/2,0,-m,0,0,0];
            product{1}=4;
        elseif k==-n  && m==0
            Counter{1}=[0,-n/2,0,0,0,0,0];
            product{1}=2;
        elseif -n-2*m<k && k<-n-1 && mod(k,2)==1
            i=(n+2*m+3+ k)/2;
            Counter{1}=[0,-n/2,-(i-2),-(m-i),0,0,1];
            product{1}=1;
        elseif -n-1<=k && k<-1
            i=(-k-3)/2;
            %x11/alambda=(2s3)/asigma so there are two ways to do this
            %So even though there are not two generators, we will still use
            %2 ideals to exproductess the two different forms
            product{1}=2;
            Counter{1}=[-(n-2*i),-i,-(m-2),0,0,0,1];
            Counter{2}=[-(n-2*i-1),-i,-(m-1),0,0,1,0];
            product{2}=1;
        end
    else
        if k==-n-2*m && m~=0
            %x12=(asigma s3)/u2sigma so multiple options
            product{1}=1;
            Counter{1}=[0,-(n-1)/2,0,-(m-1),0,1,0];
            product{2}=1;
            Counter{2}=[1,-(n+1)/2,0,-(m-2),0,0,1];
        elseif -n-2*m<k && k<-n-1 && mod(k,2)==1
            %x12=(asigma s3)/u2sigma so two ways
            i=(n+2*m+2+k)/2;
            product{1}=2;
            Counter{1}=[-1,-(n-1)/2,-(i-2),-(m-i),0,0,1];
            product{2}=1;
            Counter{2}=[1,-(n+1)/2,-(i-1),-(m-i-1),0,0,1];
        elseif -n-1<=k && k<-1 && mod(k,2)==1
            i=(-k-3)/2;
            %x11/alambda=(2s3)/asigma so two ways
            product{1}=2;
            Counter{1}=[-(n-2*i),-i,-(m-2),0,0,0,1];
            Counter{2}=[-(n-2*i-1),-i,-(m-1),0,0,1,0];
            product{2}=1;
        end

    end
end

if gen(2)<0 && gen(3)>0 %lambda-sigma
    k=gen(1); n=-gen(2); m=gen(3);
    if mod(n,2)==0
        if k==2*m-n
            Counter{1}=[0,-n/2,0,m,0,0,0];
            product{1}=1;
        elseif -n+2<=k && mod(k,2)==0
            i=(2*m-n-k)/2;
            Counter{1}=[0,-n/2,i,m-i,0,0,0];
            product{1}=1;
        elseif k==-n
            product{1}=2;
            Counter{1}=[0,-n/2,m,0,0,0,0];
        elseif -n+1<=k && k<=-3 && mod(k,2)==1
            i=(-k-1)/2;
            product{1}=1;
            Counter{1}=[-(n-2*i-1),-(i-1),m,0,1,0,0];

        end
    else
        if -n+1<=k && k<2*m-n && mod(k,2)==0
            product{1}=2;
            i=(2*m-n+1-k)/2;
            Counter{1}=[-1,-(n-1)/2,i,m-i,0,0,0];
        elseif -n<=k && k<=-3 && mod(k,2)==1
            i=(-k-1)/2;
            product{1}=1;
            Counter{1}=[-(n-2*i-1),-(i-1),m,0,1,0,0];
        end
    end
end




if gen(2)>0 && gen(3)<0 %sigma-lambda
    k=gen(1); n=gen(2); m=-gen(3);
    if mod(n,2)==0
        if 0<=k && k<=n-4 && mod(k,2)==0 && k~=n-2*m
            i=k/2;
            product{1}=1;
            Counter{1}=[n-2*i,i,-m,0,0,0,0];
        elseif k==n-2*m && k>=0 && m>=2 %Two generators
            i=k/2;
            product{1}=1;
            product{2}=4;
            Counter{1}=[n-2*i,i,-m,0,0,0,0];
            Counter{2}=[0,n/2,0,-m,0,0,0];
        elseif k==n-2*m && k<0 && m>=2
            product{1}=4;
            Counter{1}=[0,n/2,0,-m,0,0,0];
        elseif m==1 && k==n-2
            product{1}=2;
            Counter{1}=[0,n/2,0,-1,0,0,0];
        elseif n-2*m<k && k<n-3 && mod(k,2)==1
            i=(k-n+2*m+3)/2;
            product{1}=1;
            Counter{1}=[0,n/2,-(i-2),-(m-i),0,0,1];
        end
    else
        if k<=n-3 && mod(k,2)==0 &&k>=0
            i=k/2;
            product{1}=1;
            Counter{1}=[n-2*i,i,-m,0,0,0,0];
        elseif n-2*m<k && k<=n-4 && mod(k,2)==1
            i=(k-n+2*m+4)/2;
            product{1}=1;
            Counter{1}=[1,(n-1)/2,-(i-2),-(m-i),0,0,1];
        end
    end    
end