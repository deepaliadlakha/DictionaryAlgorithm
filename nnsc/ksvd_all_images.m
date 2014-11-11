function [D,rmse_image]=nnsc_all_images(param,dictsize,maxIter,numDisplay,p,mu,lambda)

    close all
    rng (0);

    numatoms=dictsize;
    
% %     Barbara 
    [X, ~] = imread('barbara.png');
    X_barbara=mat2gray(X);
    
    X = double(X)+randn(size(X,1),size(X,2))*20;
   
    
    
    X=mat2gray(X);
    save_image(X,'./output/noisy_barbara_', 0);
    

    row_lim=size(X,1);
    col_lim=size(X,2);

    row_lim=row_lim-p+1;
    col_lim=col_lim-p+1;

    data_barbara=zeros(p*p,row_lim*col_lim);

    count=1;
    for i=1:row_lim,
        for j=1:col_lim,
           window = X(i:i+p-1,j:j+p-1);
           temp=window(:);
           data_barbara(:,count)= temp;
           count=count+1;
        end
    end
    
% %     Boats 
    [X, ~] = imread('boats.png');
    X_boats=mat2gray(X);
    
    X = double(X)+randn(size(X,1),size(X,2))*15;
    
    
    X=mat2gray(X);
    save_image(X,'./output/noisy_boats_', 0);
    

    row_lim=size(X,1);
    col_lim=size(X,2);

    row_lim=row_lim-p+1;
    col_lim=col_lim-p+1;

    data_boats=zeros(p*p,row_lim*col_lim);

    count=1;
    for i=1:row_lim,
        for j=1:col_lim,
           window = X(i:i+p-1,j:j+p-1);
           temp=window(:);
           data_boats(:,count)= temp;
           count=count+1;
        end
    end
    
% %     House 
    [X, ~] = imread('house.png');
    X_house=mat2gray(X);
    
    X = double(X)+randn(size(X,1),size(X,2))*25;
    
    
    X=mat2gray(X);
    save_image(X,'./output/noisy_house_', 0);
    

    row_lim=size(X,1);
    col_lim=size(X,2);

    row_lim=row_lim-p+1;
    col_lim=col_lim-p+1;

    data_house=zeros(p*p,row_lim*col_lim);

    count=1;
    for i=1:row_lim,
        for j=1:col_lim,
           window = X(i:i+p-1,j:j+p-1);
           temp=window(:);
           data_house(:,count)= temp;
           count=count+1;
        end
    end
    
% %     Lena 
    [X, ~] = imread('lena.png');
    X_lena=mat2gray(X);
    
    X = double(X)+randn(size(X,1),size(X,2))*10;
    
    
    X=mat2gray(X);
    save_image(X,'./output/noisy_lena_', 0);
   

    row_lim=size(X,1);
    col_lim=size(X,2);

    row_lim=row_lim-p+1;
    col_lim=col_lim-p+1;

    data_lena=zeros(p*p,row_lim*col_lim);

    count=1;
    for i=1:row_lim,
        for j=1:col_lim,
           window = X(i:i+p-1,j:j+p-1);
           temp=window(:);
           data_lena(:,count)= temp;
           count=count+1;
        end
    end
    
    data=[data_barbara data_boats data_house data_lena];

    mean_data=mean(data);
    mean_data=repmat(mean_data,size(data,1),1);
    data_later=data;
    data=data-mean_data;
    var_data=var(data);
    var_data=var_data.*(var_data>0.01); %based on data

    dataSelected = data(:,logical(var_data));
    dataSelected = datasample(dataSelected,10000,2,'Replace',false);
    
 if(param==1)
     
    [A,S,~]=my_nnsc (dataSelected,dictsize,maxIter,p,numDisplay,mu,lambda);
%     [D]=my_ksvd (param,dataSelected,dictsize,maxIter,p,numDisplay,0,targetSparsity); %Threshold not clear take a clear patch
    filename=strcat('./output/D_',int2str(dictsize),'.mat');
    save(filename,'D');
    D_new=ones(size(D,1),size(D,2)+1);
    D_new=D_new/p;
    D_new(:,2:end)=D;
    
    
    
    coeff_barbara=ompCholesky(D_new,data_barbara,targetSparsity);
    coeff_boats=ompCholesky(D_new,data_boats,targetSparsity);
    coeff_house=ompCholesky(D_new,data_house,targetSparsity);
    coeff_lena=ompCholesky(D_new,data_lena,targetSparsity);
   
  
        
    
    filename=strcat('./output/coeff_barbara_',int2str(dictsize),'.mat');
    save(filename,'coeff_barbara');
    filename=strcat('./output/coeff_boats_',int2str(dictsize),'.mat');
    save(filename,'coeff_boats');
    filename=strcat('./output/coeff_house_',int2str(dictsize),'.mat');
    save(filename,'coeff_house');
    filename=strcat('./output/coeff_lena_',int2str(dictsize),'.mat');
    save(filename,'coeff_lena');
    'matrix saved'
    
 elseif(param==2)
    filename=strcat('./output/D_',int2str(dictsize),'.mat');
    load(filename,'D');
    D_new=ones(size(D,1),size(D,2)+1);
    D_new=D_new/p;
    D_new(:,2:end)=D;
    
    filename=strcat('./output/coeff_barbara_',int2str(dictsize),'.mat');
    load(filename,'coeff_barbara');
    filename=strcat('./output/coeff_boats_',int2str(dictsize),'.mat');
    load(filename,'coeff_boats');
    filename=strcat('./output/coeff_house_',int2str(dictsize),'.mat');
    load(filename,'coeff_house');
    filename=strcat('./output/coeff_lena_',int2str(dictsize),'.mat');
    load(filename,'coeff_lena');
    
    'matrix loaded'
    
%  else
%     [D]=my_ksvd (param,dataSelected,dictsize,maxIter,p,numDisplay,0,targetSparsity); %Threshold not clear take a clear patch
%     filename=strcat('./output/D_kmeans_',int2str(dictsize),'.mat');
%     save(filename,'D');
%     D_new=ones(size(D,1),size(D,2)+1);
%     D_new=D_new/p;
%     D_new(:,2:end)=D;
%     
%     coeff_barbara=ompCholesky(D_new,data_barbara,targetSparsity);
%     coeff_boats=ompCholesky(D_new,data_boats,targetSparsity);
%     coeff_house=ompCholesky(D_new,data_house,targetSparsity);
%     coeff_lena=ompCholesky(D_new,data_lena,targetSparsity);
%    
%   
%         
%     
%     filename=strcat('./output/coeff_barbara_kmeans_',int2str(dictsize),'.mat');
%     save(filename,'coeff_barbara');
%     filename=strcat('./output/coeff_boats_kmeans_',int2str(dictsize),'.mat');
%     save(filename,'coeff_boats');
%     filename=strcat('./output/coeff_house_kmeans_',int2str(dictsize),'.mat');
%     save(filename,'coeff_house');
%     filename=strcat('./output/coeff_lena_kmeans_',int2str(dictsize),'.mat');
%     save(filename,'coeff_lena');
%     'matrix kmeans saved'
     
 end
 
    %display_dictionary(D,p,numDisplay,size(D,2));

    Y_new=D_new*coeff_barbara;   
    Y=zeros(size(X,1),size(X,2));
    count=1;
    for i=1:row_lim,
        for j=1:col_lim,
            
           
           pCrosspMat=reshape(Y_new(:,count),p,p);
           Y(i+4,j+4)=pCrosspMat(4,4);
           count=count+1;
        end
    end
    
    rmse_count=[];
    rmse_image=(Y(:)-X_barbara(:)).^2;
    rmse_image=sqrt(sum(rmse_image)/(size(Y,1)*size(Y,2)));
    rmse_count=[rmse_count rmse_image];
    filename=strcat('./output/outpimage_barbara_',int2str(numatoms),'.png');
    save_image(Y,filename,0);
    
    Y_new=D_new*coeff_boats;   
    Y=zeros(size(X,1),size(X,2));
    count=1;
    for i=1:row_lim,
        for j=1:col_lim,
            
           
           pCrosspMat=reshape(Y_new(:,count),p,p);
           Y(i+4,j+4)=pCrosspMat(4,4);
           count=count+1;
        end
    end
    
    rmse_image=(Y(:)-X_boats(:)).^2;
    rmse_image=sqrt(sum(rmse_image)/(size(Y,1)*size(Y,2)));
    rmse_count=[rmse_count rmse_image];
    filename=strcat('./output/outpimage_boats_',int2str(numatoms),'.png');
    save_image(Y,filename,0);
    
    Y_new=D_new*coeff_house;   
    Y=zeros(size(X,1),size(X,2));
    count=1;
    for i=1:row_lim,
        for j=1:col_lim,
            
           
           pCrosspMat=reshape(Y_new(:,count),p,p);
           Y(i+4,j+4)=pCrosspMat(4,4);
           count=count+1;
        end
    end
    
    rmse_image=(Y(:)-X_house(:)).^2;
    rmse_image=sqrt(sum(rmse_image)/(size(Y,1)*size(Y,2)));
    rmse_count=[rmse_count rmse_image];
    filename=strcat('./output/outpimage_house_',int2str(numatoms),'.png');
    save_image(Y,filename,0);
    
    Y_new=D_new*coeff_lena;   
    Y=zeros(size(X,1),size(X,2));
    count=1;
    for i=1:row_lim,
        for j=1:col_lim,
            
           
           pCrosspMat=reshape(Y_new(:,count),p,p);
           Y(i+4,j+4)=pCrosspMat(4,4);
           count=count+1;
        end
    end
    
    rmse_image=(Y(:)-X_lena(:)).^2;
    rmse_image=sqrt(sum(rmse_image)/(size(Y,1)*size(Y,2)));
    rmse_count=[rmse_count rmse_image];
    filename=strcat('./output/outpimage_lena_',int2str(numatoms),'.png');
    save_image(Y,filename,0);
    
    
    
    
    
    
    %save_image(X, filename, 0);
    
    %display_dictionary(D,p,numDisplay,size(D,2));
    %imshow(Y);
    %pause
    
    
   
end