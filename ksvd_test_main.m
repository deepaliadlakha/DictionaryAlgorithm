function [Y,D,coeff,rmse_image]=ksvd_test_main(param,dictsize,maxIter,numDisplay,p)

    close all
    rng (0);
    
    
    features=zeros(9,10);
    band=ones(1,3);
    
    one_features=zeros(3,3);
    one_features(1,:)=band;
    features(:,1)=one_features(:);
    one_features=1-one_features;
    features(:,7)=one_features(:); 
    one_features=one_features';
    features(:,4)=one_features(:);
    one_features=1-one_features;
    features(:,8)=one_features(:);
    
    one_features=zeros(3,3);
    one_features(2,:)=band;
    features(:,2)=one_features(:);
    one_features=one_features';
    features(:,5)=one_features(:);
    
    one_features=zeros(3,3);
    one_features(3,:)=band;
    features(:,3)=one_features(:);
    one_features=1-one_features;
    features(:,9)=one_features(:);
    one_features=one_features';
    features(:,6)=one_features(:);
    one_features=1-one_features;
    features(:,10)=one_features(:);
    
    features=normc(features);
    
     n=10;
    
%     features = datasample(features,3,2,'Replace',false);
%     n=3;
    
    p=3;
    numDisplay=5;
    bins=linspace(1,n,n);
    combinations_three=nchoosek(bins,3);
    combinations_two=nchoosek(bins,2);
    combinations_one=nchoosek(bins,1);
    data=zeros(9,size(combinations_three,1)+size(combinations_two,1)+size(combinations_one,1));
%     data=zeros(9,size(combinations_three,1)+size(combinations_two,1));
    
    for i=1:size(combinations_three,1)
        three=features(:,combinations_three(i,:));
        data(:,i)=three(:,1)+three(:,2)+three(:,3);
    end
    
    for i=1:size(combinations_two,1)
        two=features(:,combinations_two(i,:));
        data(:,size(combinations_three,1)+i)=two(:,1)+two(:,2);
    end
    
    for i=1:size(combinations_one,1)
        one=features(:,combinations_one(i,:));
        data(:,size(combinations_three,1)+size(combinations_two,1)+i)=one(:,1);
    end
    
    dataSelected=normc(data);
    

%     numatoms=dictsize;
%     [X, ~] = imread('lena.png');
%     save_image(X, './output/original_lena.png', 0);
%     X = double(X)+randn(size(X,1),size(X,2))*5;
%     X=mat2gray(X);
%     save_image(X, './output/noisy_lena_0_1.png', 0);
%     
%     %X=X/255;
%     
% %     X=X(1:200,1:200);
% 
%     row_lim=size(X,1);
%     col_lim=size(X,2);
% 
%     row_lim=row_lim-p+1;
%     col_lim=col_lim-p+1;
% 
%     data=zeros(p*p,row_lim*col_lim);
% 
%     count=1;
%     for i=1:row_lim,
%         for j=1:col_lim,
%            window = X(i:i+p-1,j:j+p-1);
%            temp=window(:);
%            data(:,count)= temp;
%            count=count+1;
%         end
%     end
% 
%     mean_data=mean(data);
%     mean_data=repmat(mean_data,size(data,1),1);
%     data_later=data;
%     data=data-mean_data;
%     var_data=var(data);
%     var_data=var_data.*(var_data>0.01); %based on data
% 
%     
%     
%     dataSelected = data(:,logical(var_data));
    'size'
    size(dataSelected)
    numDisplay=3;
 if(param==1)
     
    
    [D]=my_ksvd (param,dataSelected,dictsize,maxIter,p,numDisplay); %Threshold not clear take a clear patch
    filename=strcat('./output/ksvd_test/D_',int2str(dictsize),'.mat');
    save(filename,'D');
%     D_new=ones(size(D,1),size(D,2)+1);
%     D_new=D_new/p;
%     D_new(:,2:end)=D;
    D_new=D;
    display_dictionary(D_new,p,numDisplay);
    colorbar;
    pause
    
    coeff=ompCholesky(D_new,dataSelected);
    Y_new=D_new*coeff;
    display_dictionary(Y_new(:,1:5),p,numDisplay);
    colorbar;
    norm(Y_new-dataSelected,'fro')
    pause
    
        
    
    filename=strcat('./output/coeff_',int2str(dictsize),'.mat');
    save(filename,'coeff');
    'matrix saved'
    
 elseif(param==2)
    filename=strcat('./output/D_',int2str(dictsize),'.mat');
    load(filename,'D');
    D_new=ones(size(D,1),size(D,2)+1);
    D_new=D_new/p;
    D_new(:,2:end)=D;
    filename=strcat('./output/coeff_',int2str(dictsize),'.mat');
    load(filename,'coeff');
    'matrix loaded'
    
 else
    [D]=my_ksvd (param,dataSelected,dictsize,maxIter,p,numDisplay); %Threshold not clear take a clear patch
    filename=strcat('./output/D_kmeans_',int2str(dictsize),'.mat');
    save(filename,'D');
    D_new=ones(size(D,1),size(D,2)+1);
    D_new=D_new/p;
    D_new(:,2:end)=D;
    
    
    coeff=ompCholesky(D_new,data_later);
    
    
    filename=strcat('./output/coeff_kmeans_',int2str(dictsize),'.mat');
    save(filename,'coeff');
    'matrix kmeans ------saved'
     
 end
 
    %display_dictionary(D,p,numDisplay,size(D,2));

    Y_new=D_new*coeff;   
    
%%%%%%%%%%%%%%% Averaging over the pixels %%%%%%%%%%%%%%%%%%%%%%%%%   
%     norm(dataSelected-D_new*coeff,'fro')
%     norm(dataSelected(:,1:5:end)-D_new*coeff(:,1:5:end),'fro')
%     imagesc (coeff); colorbar; pause, close    
%     count=1; 
%     Y=zeros(size(X,1),size(X,2));
%     count_Y=zeros(size(X,1),size(X,2));
%     all_ones=ones(p,p);
%     for i=1:row_lim,
%         for j=1:col_lim,
%            Y(i:i+p-1,j:j+p-1)=Y(i:i+p-1,j:j+p-1)+reshape(Y_new(:,count),p,p);
%            count_Y(i:i+p-1,j:j+p-1)=count_Y(i:i+p-1,j:j+p-1)+all_ones;
%            count=count+1;
%         end
%     end
    
%     
%     Y=zeros(size(X,1),size(X,2));
%     count=1;
%     for i=1:row_lim,
%         for j=1:col_lim,
%             
%            
%            pCrosspMat=reshape(Y_new(:,count),p,p);
%            Y(i+4,j+4)=pCrosspMat(4,4);
%            count=count+1;
%         end
%     end
%     
%     rmse_image=(Y(:)-X(:)).^2;
%     rmse_image=sqrt(sum(rmse_image)/(size(Y,1)*size(Y,2)));
%     
%     filename=strcat('./output/outpimage_',int2str(numatoms),'.png');
%     save_image(Y, filename, 0);
%     
%     display_dictionary(D,p,numDisplay);
%     colorbar;
%     pause
    
    
   
end