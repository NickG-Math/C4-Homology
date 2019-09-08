function transferred=ranktransfer(rank,originallevel) 
%Inputs: Array rank and int originallevel
%Outputs: Array transferred
%Description: Transfers rank on level higher. originallevel is 1 at bottom.
%Currently works only for C4.

transferlevel=2*originallevel;
halveit=(rank>4/transferlevel);
transferred=rank;
transferred(halveit)=transferred(halveit)/2;  %so a 4 in G/e transfers to a 2, and a 2 in G/C2 transfers to a 1
end