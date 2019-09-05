clear

Data=struct();

C4_Mackey_List=load('Data\C4_Mackey_List.mat');
Data.MackeyList=C4_Mackey_List.MackeyList;

C4_Standard_Chains=load('Data\C4_Standard_Chains.mat');
Data.rankStandard=C4_Standard_Chains.rankStandard;
Data.ChainsStandard=C4_Standard_Chains.ChainsStandard;

C4_Change_Basis=load('Data\C4_Change_Basis.mat');

Data.longprimelist=C4_Change_Basis.longprimelist;
Data.linearToMatrix=C4_Change_Basis.linearToMatrix;
Data.primeToLinear=C4_Change_Basis.primeToLinear;
Data.ChangeBasis=C4_Change_Basis.ChangeBasis;

clear C4_Mackey_List C4_Standard_Chains C4_Mackey_List;