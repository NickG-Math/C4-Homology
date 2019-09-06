function write_C4_Mackey_List
%action, Res, Transfer, Homology
MackeyList{1}={[0,0,0,0,0,0,0,0,0],"0"};
MackeyList{2}={[1,1,1,1,2,2,1,1,1],"Z"};
MackeyList{3}={[-1,-1,0,1,0,2,0,1,1],"Z_-"};
MackeyList{4}={[1,1,2,1,1,2,1,1,1],"p*L"};
MackeyList{5}={[-1,-1,0,1,1,2,2,1,1],"p*L_-"};
MackeyList{6}={[1,1,2,2,1,1,1,1,1],"L"};
MackeyList{7}={[-1,-1,0,2,1,1,2,1,1],"L_-"};
MackeyList{8}={[0,0,0,0,0,0,2,0,0],"Z/2"};
MackeyList{9}={[1,0,1,0,2,0,4,2,0],"Z/4"};
MackeyList{10}={[1,0,0,0,0,0,0,2,0],"V"};
MackeyList{11}={[1,0,0,0,1,0,2,2,0],"Q"};
MackeyList{12}={[1,0,0,0,0,0,2,2,0],"Z/2+V"};
MackeyList{13}={[-1,-1,0,1,0,2,2,1,1],"Z/2+Z_-"};
MackeyList{14}={[1,1,1,2,2,1,1,1,1],"L^o"};
MackeyList{15}={[-1,-1,0,2,0,1,0,1,1],"Z^o_-"};
MackeyList{16}={[1,0,1,0,0,0,2,2,0],"Q^o"};

MackeyList{17}={{[1,1,2,1,1,1], [2,1], [0,1], [0,2]} "L+Z/2"}; %action, Res(2), Tr(2), H(2),H(1) in the first matrix, H(4), Tr(4), Res(4) in the second
MackeyList{18}={{[1,1,2,1,1,1], [2,1], [1,1], [0,2]} "L+Z/2"}; %action, Res(2), Tr(2), H(2),H(1) in the first matrix, H(4), Tr(4), Res(4) in the second


save('C4_Mackey_List.mat');