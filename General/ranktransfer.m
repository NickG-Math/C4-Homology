function rank=ranktransfer(rank,level,order) 
%
%rank=RANKTRANSFER(rank,level,order) 
%
%INPUT: Array rank and ints level,order
%
%OUTPUT: Array rank
%
%DESCRIPTION: Transfers rank from bottom level to given level. 
%order is the order of the cyclic group.

rank(rank*level>order)=order/level;
end