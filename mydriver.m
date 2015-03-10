function [x_plot,y_plot,y_plotk]=mydriver(param,minnumatoms,maxnumatoms)
addpath('./common/export_fig/')
addpath('./common/')


x_plot=zeros(1,maxnumatoms-minnumatoms+1);
y_plot=zeros(1,maxnumatoms-minnumatoms+1);
y_plotk=zeros(1,maxnumatoms-minnumatoms+1);

for i=minnumatoms:maxnumatoms
    [~,D,~,rmse_image]=ksvd_main(param,5*i,5,5,8);
    x_plot(1,i)=size(D,2);
    y_plot(1,i)=rmse_image;
    [~,~,~,rmse_image]=ksvd_main(3,size(D,2),5,5,8);
    y_plotk(1,i)=rmse_image;
end



