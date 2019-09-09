function transferred=ranktransfer(rank,originallevel,n) 
%Inputs: Array rank and ints originallevel, n
%Outputs: Array transferred
%Description: Transfers rank from originallevel to originallevel*2. originallevel is 1 at bottom. %If our group is G=C2^n then n=2 

halveit=(rank>n/(2*originallevel)); %Finds where rank>(n/2)*originallevel. These will need to be halved.  %eg for n=4, a 4 in C4/e transfers to a 2, and a 2 in C4/C2 transfers to a 1
transferred=rank;
transferred(halveit)=transferred(halveit)/2; %We halve them in a vectorized command
end