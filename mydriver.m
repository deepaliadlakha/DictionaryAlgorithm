function mydriver(param,minnumatoms,maxnumatoms,sparseParam,targetSparsity)
addpath('./common/export_fig/')
addpath('./common/')


x_plot=zeros(1,maxnumatoms-minnumatoms+1);
y_plot=zeros(1,maxnumatoms-minnumatoms+1);

for i=minnumatoms:maxnumatoms
    [~,~,~,rmse_image]=ksvd_main(param,5*i,5,5,8,sparseParam,targetSparsity);
    x_plot(1,i)=5*i;
    y_plot(1,i)=rmse_image;
   
end


% scatter(x_plot,y_plot);
% pause;


