clc;
clear all;

data = importdata('recon_STR03713.mul');
data = data.data;

n_chan = size(data,2);
n_iter = round(size(data,1)/4096);

for i = 1:n_iter
    start = 1 + (i-1)*4096;
    endd =  i*4096;
    data_for_eph = data(start:endd,:);

    filename = strcat('cond', num2str(i),'.eph');

    row1 = [n_chan, 4096, 2048];

    dlmwrite(filename,row1,'delimiter','\t');
    dlmwrite(filename,data_for_eph','-append','delimiter','\t');
end

