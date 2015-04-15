clear ;
data_fname = '../B1.dat' ;
%load('CBlist.mat');
file_id = fopen(data_fname, 'rb');
L=4000*5e3;%每次读取的数据量
cdata = zeros(1,L) ; %存放读出的数据（复信号）
state=zeros(3,37);
for m=1:1
    [row_array, ele_count] = fread(file_id, L*2, 'int8') ;
    if ele_count < L*2
        break ;
    else
        cdata(:)=row_array(1:2:2*L)+row_array(2:2:2*L)*1i;
    end
    save('data.mat','cdata');
end

fclose(file_id);

clear file_id row_array L ele_count data_fname
