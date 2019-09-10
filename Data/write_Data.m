function write_Data(maxlengthA,maxlengthB)
%
%INPUT: ints maxlengthA,maxlengthB
%
%DESCRIPTION: Calls all "write" functions in the Data folder.

%write_C4_Standard_Chains(rangeN,rangeM); %Now deprecated
write_C4_Mackey_List;
write_Look_Up(2,maxlengthA,maxlengthB);

end