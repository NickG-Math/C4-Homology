function rank=ranktransfer(rank,level,order) 
%
%INPUT: Array rank and ints level,order
%
%OUTPUT: Array transferred
%
%DESCRIPTION: Transfers rank from level 1 to given level. order is the order of the cyclic group

rank(rank*level>order)=order/level;
end