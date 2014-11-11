function [y_plot]=mydriver(param,minnumatoms,maxnumatoms)
addpath('../common/export_fig/')
addpath('../common/')
addpath('../')

y_plot=zeros(1,maxnumatoms-minnumatoms+1);

for i=minnumatoms:maxnumatoms
    [~,rmse_image]=nnsc_all_images(param,5*i,50,1,8,0);
   
    if(i==minnumatoms)
        y_plot=rmse_image';
    else
        y_plot=[y_plot rmse_image'];
    end
   
end





