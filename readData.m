clear
%load('data.mat');
data_fname = '../B1.dat' ;
file_id = fopen(data_fname, 'rb');
L=4*5e3;%每次读取的数据量
load('CBlist.mat');%读取CB1码表
prn_id=4;
CB=CBs(prn_id,:);

%    data_buffer=row_array(1:2:2*L)'+row_array(2:2:2*L)'*1i;


for m=1:20
    [row_array, ~] = fread(file_id, L*2, 'int8') ;
    data_buffer=row_array(1:2:2*L)'+row_array(2:2:2*L)'*1i;
    
    [freq_i,code_phase,rate]=catchB1([CB,zeros(1,3*2046)],data_buffer,0);

end

fclose(file_id);