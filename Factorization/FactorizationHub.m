function [product,CounterNum,CounterDen,found,VisitedNum,VisitedDen]=FactorizationHub(gen,SwitchLimit,NumOrDen,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,product,CounterNum,CounterDen,VisitedNum,VisitedDen,k)
%[product,CounterNum,CounterDen,found,VisitedNum,VisitedDen]=FactorizationHub(gen,SwitchLimit,NumOrDen,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,product,CounterNum,CounterDen,VisitedNum,VisitedDen,k)
%
%INPUTS: array gen, ints SwitchLimit,NumOrDen,safety, cells
%Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles, int product,
%arrays of ints CounterNum, CounterDen, logical matrices VisitedNum,VisitedDen, int k
%
%OUTPUTS: int product, arrays of ints CounterNum, logical found, matrices
%VisitedNum, VisitedDen
%
%DESCRIPTION: FactorizationHub assigns the functions FactorizeNum and FactorizeDen to find
%factorizations for the numerator and denominators through a recursion 
%
%NumOrDen determines which of the two functions to use. SwitchLimit limits
%how many switches from FactorizeNum to FactorizeDen are possible in a
%single search (otherwise the recursion could be infinite)
%safety enables usage of VisitedNum,VisitedDen matrices that record which 
%generators we already checked through the recursion so we don't do it
%again (otherwise we can get infinite loops)
%found=1 if we found the factorization, 0 otherwise

do_the_other_as_well=1; % If Num/Den doesn't work try it with the other.

VisitedNumpre=VisitedNum;%Remember in case the search does not pan out, so we can revert to this
VisitedDenpre=VisitedDen;

while SwitchLimit>=0 && do_the_other_as_well>=0 %If we still have switches
    if NumOrDen==1 %Factorize numerator
        [gen,product,CounterNum,CounterDen,found,VisitedNum,VisitedDen]=FactorizeNum(gen,SwitchLimit-1,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,product,CounterNum,CounterDen,VisitedNum,VisitedDen,k);
        if found
            return %Got it!
        else %Do the other
            NumOrDen=-1;
            do_the_other_as_well=do_the_other_as_well-1;
        end
    else
        [gen,product,CounterNum,CounterDen,found,VisitedNum,VisitedDen]=FactorizeDen(gen,SwitchLimit-1,safety,Table,IndexedHomology,BasicIrreducibles,MoreIrreducibles,product,CounterNum,CounterDen,VisitedNum,VisitedDen,k);
        if found
            return
        else
            NumOrDen=1;
            do_the_other_as_well=do_the_other_as_well-1;
        end
    end
    VisitedNum=VisitedNumpre; %We don't want to influence the completely different iteration
    VisitedDen=VisitedDenpre;
end
%In case we run out of switches
found=0;
%The following values are requested but never used, so they don't matter
product=0; CounterNum=0; CounterDen=0;