clear
%load('data.mat');
data_fname = '../B1.dat' ;
file_id = fopen(data_fname, 'rb');
L=5e3;%ÿ�ζ�ȡ��������
load('CBlist.mat');%��ȡCB1���
prn_id=1;
CB=CBs(prn_id,:);

%    data_buffer=row_array(1:2:2*L)'+row_array(2:2:2*L)'*1i;
for m=1:10
    [row_array, ~] = fread(file_id, L*2, 'int8') ;
    data_buffer=row_array(1:2:2*L)'+row_array(2:2:2*L)'*1i;
    [freq_i,code_phase,rate]=catchB1(CB,data_buffer)
end
fclose(file_id);