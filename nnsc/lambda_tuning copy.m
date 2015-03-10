function [error_val1]=lambda_tuning(dictsize,maxIter,numDisplay,p)

    addpath('../common/export_fig/')
    addpath('../common/')
    addpath('../')

    close all
     rng (0);
    
    
    noise_level=100;

% %     Barbara 
    [X, ~] = imread('boats.png');
    X_old=mat2gray(X);
    
%     X = double(X)+randn(size(X,1),size(X,2))*noise_level;
   
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
    
    data_barbara_old=zeros(p*p,row_lim*col_lim);
    count=1;
    for i=1:row_lim,
        for j=1:col_lim,
           window = X_old(i:i+p-1,j:j+p-1);
           temp=window(:);
           data_barbara_old(:,count)= temp;
           count=count+1;
        end
    end
    
%     
%     
% % %     Boats 
%     [X, ~] = imread('boats.png');
%     X_old=mat2gray(X);
%    
%     X = double(X)+randn(size(X,1),size(X,2))*noise_level;
%     
%     
%     X=mat2gray(X);
%     save_image(X,'./output/noisy_boats_', 0);
%     
% 
%     row_lim=size(X,1);
%     col_lim=size(X,2);
% 
%     row_lim=row_lim-p+1;
%     col_lim=col_lim-p+1;
% 
%     data_boats=zeros(p*p,row_lim*col_lim);
% 
%     count=1;
%     for i=1:row_lim,
%         for j=1:col_lim,
%            window = X(i:i+p-1,j:j+p-1);
%            temp=window(:);
%            data_boats(:,count)= temp;
%            count=count+1;
%         end
%     end
%     
%     data_boats_old=zeros(p*p,row_lim*col_lim);
%     count=1;
%     for i=1:row_lim,
%         for j=1:col_lim,
%            window = X_old(i:i+p-1,j:j+p-1);
%            temp=window(:);
%            data_boats_old(:,count)= temp;
%            count=count+1;
%         end
%     end
%     
% % %     House 
%     [X, ~] = imread('house.png');
%     X_old=mat2gray(X);
%    
%     X = double(X)+randn(size(X,1),size(X,2))*noise_level;
%     
%     
%     X=mat2gray(X);
%     save_image(X,'./output/noisy_house_', 0);
%     
% 
%     row_lim=size(X,1);
%     col_lim=size(X,2);
% 
%     row_lim=row_lim-p+1;
%     col_lim=col_lim-p+1;
% 
%     data_house=zeros(p*p,row_lim*col_lim);
% 
%     count=1;
%     for i=1:row_lim,
%         for j=1:col_lim,
%            window = X(i:i+p-1,j:j+p-1);
%            temp=window(:);
%            data_house(:,count)= temp;
%            count=count+1;
%         end
%     end
%     
%     data_house_old=zeros(p*p,row_lim*col_lim);
%     count=1;
%     for i=1:row_lim,
%         for j=1:col_lim,
%            window = X_old(i:i+p-1,j:j+p-1);
%            temp=window(:);
%            data_house_old(:,count)= temp;
%            count=count+1;
%         end
%     end
%     
% % %     Lena 
%     [X, ~] = imread('lena.png');
%     X_old=mat2gray(X);
%     
%     X = double(X)+randn(size(X,1),size(X,2))*noise_level;
%     
%     
%     X=mat2gray(X);
%     save_image(X,'./output/noisy_lena_', 0);
%    
% 
%     row_lim=size(X,1);
%     col_lim=size(X,2);
% 
%     row_lim=row_lim-p+1;
%     col_lim=col_lim-p+1;
% 
%     data_lena=zeros(p*p,row_lim*col_lim);
% 
%     count=1;
%     for i=1:row_lim,
%         for j=1:col_lim,
%            window = X(i:i+p-1,j:j+p-1);
%            temp=window(:);
%            data_lena(:,count)= temp;
%            count=count+1;
%         end
%     end
%     
%     data_lena_old=zeros(p*p,row_lim*col_lim);
%     count=1;
%     for i=1:row_lim,
%         for j=1:col_lim,
%            window = X_old(i:i+p-1,j:j+p-1);
%            temp=window(:);
%            data_lena_old(:,count)= temp;
%            count=count+1;
%         end
%     end
%     
%     data=[data_barbara data_boats data_house data_lena];
%     data_old=[data_barbara_old data_boats_old data_house_old data_lena_old];

    dataSelected=data_barbara;
    dataSelected = double(dataSelected)+randn(size(dataSelected,1),size(dataSelected,2))*0.5;
    dataSelected=dataSelected.*(dataSelected>0);
    data_old=data_barbara_old;
    
%     dataSelected=data;
%     minInt=ones(1,size(data,2));
%     for col=1:size(data,2)
%         minInt(1,col)=min(data(:,col));
%         dataSelected(:,col)=dataSelected(:,col)-minInt(1,col);
% 
%     end
%     data_old=mat2gray(dataSelected,prctile(dataSelected(:),[0,100]));
%     data_old
%     pause
%     data_old=normc(data_old);
% %     dataSelected = double(dataSelected)+randn(size(dataSelected,1),size(dataSelected,2))*0;
% %     dataSelected=dataSelected.*(dataSelected>0);
%     dataSelected=mat2gray(dataSelected,prctile(dataSelected(:),[7,93]));
%     dataSelected=normc(dataSelected);
    
    
    
    var_data=var(data_old);
    var_data=var_data.*(var_data>0.07); %based on data

    dataSelected = dataSelected(:,logical(var_data));
    
    data_oldSelected = data_old(:,logical(var_data));
    [dataSelected, order] = datasample(dataSelected,500,2,'Replace',false);
    data_oldSelected=data_oldSelected(:,order);
    m=size(dataSelected,2);
    
    rem=m;
    mask=zeros(5,m);
   
    mask(1,1:floor(rem/5))=ones(1,floor(rem/5));
    allot=floor(rem/5);
    rem=m-allot;
    
    mask(2,allot+1:allot+floor(rem/4))=ones(1,floor(rem/4));
    allot=allot+floor(rem/4);
    rem=m-allot;
    
    mask(3,allot+1:allot+floor(rem/3))=ones(1,floor(rem/3));
    allot=allot+floor(rem/3);
    rem=m-allot;
    
    mask(4,allot+1:allot+floor(rem/2))=ones(1,floor(rem/2));
    allot=allot+floor(rem/2);
    rem=m-allot;
    
    mask(5,allot+1:allot+rem)=ones(1,rem);
    
    lambda_val=[0,0.001,0.01,0.05,0.1,0.5,1,1.5,2,2.5,3,3.5,4];
    error_val=zeros(size(lambda_val,1),size(lambda_val,2));
    
    
    for i=1:size(lambda_val,2)
     for j=1:5
        train_set=dataSelected(:,~logical(mask(j,:)));
        val_set=dataSelected(:,logical(mask(j,:)));
        val_set_true=data_oldSelected(:,logical(mask(j,:)));
        
        [D,~,~]=my_nnsc (train_set,dictsize,maxIter,p,numDisplay,lambda_val(i));
%         D_new=ones(size(D,1),size(D,2)+1);
%         D_new=D_new/p;
%         D_new(:,2:end)=D;
%         D=D_new;
        
           coeff_val=rand(size(D,2),size(val_set,2));
%         coeff_val=ompCholesky(D,val_set,3);
% %         coeff_val
%        
         for dummy_c=1:100
             coeff_val=(coeff_val.*(D'*val_set))./((D'*D)*coeff_val+lambda_val(i));
         end
% %         coeff_val
% %         pause
         error_val(i)=error_val(i)+norm(val_set_true-D*coeff_val,'fro')/numel(val_set);
     end
     
      subplot(1,4,1),  display_dictionary(val_set_true(:,1:10),8,2);
      subplot(1,4,2),  display_dictionary(val_set(:,1:10),8,2);
      subplot(1,4,3), display_dictionary(D*coeff_val(:,1:10),8,2);
      coeff_val;
      error_val(i)
      subplot(1,4,4), display_dictionary(D,8,2);
      pause
     
%       display_dictionary(D,8,1);
    end
    
    error_val1=error_val;
%     error_val=zeros(size(lambda_val,1),size(lambda_val,2));
%     
%     dictsize=20;
%     for i=1:size(lambda_val,2)
%      for j=1:5
%         train_set=dataSelected(:,~logical(mask(j,:)));
%         val_set=data_oldSelected(:,logical(mask(j,:)));
%         
%         [D,~,~]=my_nnsc (train_set,dictsize,maxIter,p,numDisplay,lambda_val(i));
%         
%         coeff_val=rand(size(D,2),size(val_set,2));
%         
%        
%         for dummy_c=1:100
%             coeff_val=(coeff_val.*(D'*val_set))./((D'*D)*coeff_val+lambda_val(i));
%         end
%         error_val(i)=error_val(i)+norm(val_set-D*coeff_val,'fro')/numel(val_set);
%      end
% %      display_dictionary(D,8,1);
% %      pause
%     end
    
%     error_val2=error_val;
    
    
      

   
end