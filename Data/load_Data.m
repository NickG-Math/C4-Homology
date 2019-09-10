function Data=load_Data()
%
%OUTPUT: struct Data
%
%DESCRIPTION: Loads all the data variables into a single structure Data. This is passed as inpu to most of our functions 

Data=struct(); %The variable that has all the data

C4_Mackey_List=load('C4_Mackey_List.mat');
Data.MackeyList=C4_Mackey_List.MackeyList;

% C4_Standard_Chains=load('C4_Standard_Chains.mat'); %Now deprecated
% Data.rankStandard=C4_Standard_Chains.rankStandard;
% Data.ChainsStandard=C4_Standard_Chains.ChainsStandard;

C4_Change_Basis=load('C4_Change_Basis.mat');
Data.ChangeBasis=C4_Change_Basis.ChangeBasis;