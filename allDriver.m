function [D_num,y_plot]=allDriver(param,minnumatoms,maxnumatoms,targetSparsity)
addpath('./common/export_fig/')
addpath('./common/')


x_plot=zeros(1,maxnumatoms-minnumatoms+1);


num_atom=[ 5, 12, 15, 20, 25, 30, 35, 40, 45];
for i=1:9
    [D,rmse_image]=ksvd_all_images(param,num_atom(i),5,5,8,targetSparsity);
%     x_plot(1,i)=5*i;
    if(i==minnumatoms)
        y_plot=rmse_image';
        D_num=size(D,2);
    else
        y_plot=[y_plot rmse_image'];
        D_num=[D_num size(D,2)];
    end
   
end

end

