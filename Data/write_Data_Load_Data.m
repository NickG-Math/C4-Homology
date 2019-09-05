function write_Data_Load_Data(rangeN,rangeM,maxlength)

%First we write
write_C4_Standard_Chains(rangeN,rangeM);
write_C4_Mackey_List;
write_C4_Change_Basis(maxlength);


%Then we load
load_Data;