function [y_plot]=alldriver(param,minnumatoms,maxnumatoms,targetSparsity)
addpath('./common/export_fig/')
addpath('./common/')


x_plot=zeros(1,maxnumatoms-minnumatoms+1);


for i=minnumatoms:maxnumatoms
    [~,rmse_image]=ksvd_all_images(param,5*i,5,5,8,targetSparsity);
%     x_plot(1,i)=5*i;
    if(i==minnumatoms)
        y_plot=rmse_image';
    else
        y_plot=[y_plot rmse_image'];
    end
   
end

end

