clear all
close all
rng (0);

[X, map] = imread('barbara.png');
X=single(X);
X=X/255;

p=8;

row_lim=size(X(1:128,1:128),1);
col_lim=size(X(1:128,1:128),2);

row_lim=row_lim-p+1;
col_lim=col_lim-p+1;

data=zeros(p*p,row_lim*col_lim);

% 

count=1;
for i=1:row_lim,
    for j=1:col_lim,
       window = X(i:i+p-1,j:j+p-1);
       temp=window(:);
       data(:,count)= temp;
       count=count+1;
    end
end

my_ksvd (data,9,21);

