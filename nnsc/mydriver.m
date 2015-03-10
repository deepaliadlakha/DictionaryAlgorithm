function [y_plot,y_plot2]=mydriver(param,minnumatoms,maxnumatoms)
addpath('../common/export_fig/')
addpath('../common/')
addpath('../')

y_plot=zeros(1,maxnumatoms-minnumatoms+1);

for i=minnumatoms:maxnumatoms
    [~,~,rmse_image]=nnsc_all_images_log(param,2,5*i,50,1,8,0.01,0.01);
   
    if(i==minnumatoms)
        y_plot=rmse_image';
    else
        y_plot=[y_plot rmse_image'];
    end
   
end

y_plot2=zeros(1,maxnumatoms-minnumatoms+1);

for i=minnumatoms:maxnumatoms
    [~,~,rmse_image]=nnsc_all_images_log(param,1,5*i,50,1,8,0.01,0.01);
   
    if(i==minnumatoms)
        y_plot2=rmse_image';
    else
        y_plot2=[y_plot2 rmse_image'];
    end
   
end





