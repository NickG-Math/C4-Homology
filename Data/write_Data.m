function write_Data(rangeN,rangeM,maxlengthA,maxlengthB)
%Inputs: ints rangeN,rangeM,maxlengthA,maxlengthB
%Description: Calls all "write" functions in the Data folder.

write_C4_Standard_Chains(rangeN,rangeM);
write_C4_Mackey_List;
write_Look_Up(2,maxlengthA,maxlengthB);

end