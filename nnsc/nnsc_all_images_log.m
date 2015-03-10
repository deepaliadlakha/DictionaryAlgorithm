function [D_new,coeff_barbara,rmse_count]=nnsc_all_images_log(param,second_param,dictsize,maxIter,numDisplay,p,lambda,noiseLevel)

    addpath('../common/export_fig/')
    addpath('../common/')
    addpath('../')

    close all
    rng (0);
    numatoms=dictsize;
    

% %     Barbara 
    [X, ~] = imread('lena.png');
    X_old=mat2gray(X);
    
    X = double(X_old)+randn(size(X_old,1),size(X_old,2))*noiseLevel;
    X=mat2gray(X,[0 1]); 
    
    if second_param==1
        save_image(X,'./output_log/noisy_barbara_', 0);
    else
        save_image(X,'./output/noisy_barbara_', 0);
    end
    

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
    
    

    dataSelected=data_barbara;
    data_old=data_barbara_old;
    
    minInt=ones(1,size(data_barbara,2));
    for col=1:size(data_barbara,2)
        minInt(1,col)=min(data_barbara(:,col));
%         dataSelected(:,col)=dataSelected(:,col)-minInt(1,col);
    end


%     dataSelected=data_barbara;
%     dataSelected = double(dataSelected)+randn(size(dataSelected,1),size(dataSelected,2))*0.5;
%     dataSelected=dataSelected.*(dataSelected>0);
%     data_old=data_barbara_old;
    
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
    var_data=var_data.*(var_data>0.05); %based on data

    dataSelected = dataSelected(:,logical(var_data));
    
    
%     data_oldSelected = data_old(:,logical(var_data));
    [dataSelected, ~] = datasample(dataSelected,1000,2,'Replace',false);
%     data_oldSelected=data_oldSelected(:,order);
    
 if(param==1)
     
     
    if(second_param==1)
        [D,~,~]=nnsc_and_log(dataSelected,dictsize,maxIter,p,numDisplay,lambda);
        filename=strcat('./output_log/D_',int2str(dictsize),'.mat');
        save(filename,'D');
    else
        [D,~,~]=my_nnsc(dataSelected,dictsize,maxIter,p,numDisplay,lambda);
        filename=strcat('./output/D_',int2str(dictsize),'.mat');
        save(filename,'D');
    end
        
%     [D]=my_ksvd (param,dataSelected,dictsize,maxIter,p,numDisplay,0,targetSparsity); %Threshold not clear take a clear patch
%     filename=strcat('./output/D_',int2str(dictsize),'.mat');
    save(filename,'D');
    D_new=ones(size(D,1),size(D,2)+1);
    D_new=D_new/p;
    D_new(:,2:end)=D;
    D_new=D;%%%%%%%%%%%%%%%%%%
    
    coeff_barbara=ones(size(D_new,2),size(data_barbara,2));
%     for i=1:100
%         coeff_barbara=(coeff_barbara.*(D_new'*data_barbara))./((D_new'*D_new)*coeff_barbara+lambda);
%     end

if(second_param==1)
    S=coeff_barbara;
    A=D_new;
    eps=0.00001;
   
    Sn=S;
    for cols=1:size(S,2)
        for ii=1:100
            W=diag(S(:,cols)+eps);
            An=A*W;
            S(:,cols)=(S(:,cols).*(An'*data_barbara(:,cols)))./((An'*An)*S(:,cols)+lambda);
        end
        Sn(:,cols)=W*S(:,cols);
    end
    
    S=Sn;
    coeff_barbara=S;
else
    for i=1:100
        coeff_barbara=(coeff_barbara.*(D_new'*data_barbara))./((D_new'*D_new)*coeff_barbara+lambda);
    end
end
    
%     coeff_boats=rand(size(D_new,2),size(data_boats,2));
%     for i=1:100
%         coeff_boats=(coeff_boats.*(D_new'*data_boats))./((D_new'*D_new)*coeff_boats+lambda);
%     end
%     
%     coeff_house=rand(size(D_new,2),size(data_house,2));
%     for i=1:100
%         coeff_house=(coeff_house.*(D_new'*data_house))./((D_new'*D_new)*coeff_house+lambda);
%     end
%     
%     coeff_lena=rand(size(D_new,2),size(data_lena,2));
%     for i=1:100
%         coeff_lena=(coeff_lena.*(D_new'*data_lena))./((D_new'*D_new)*coeff_lena+lambda);
%     end
%    
%     B=find(coeff_barbara<0);
%     assert(size(B,1)==0);
%     
%     B=find(coeff_boats<0);
%     assert(size(B,1)==0);
%     
%     B=find(coeff_house<0);
%     assert(size(B,1)==0);
%     
%     B=find(coeff_lena<0);
%     assert(size(B,1)==0);
      
    
    
    
%     coeff_barbara=ompCholesky(D_new,data_barbara,targetSparsity);
%     coeff_boats=ompCholesky(D_new,data_boats,targetSparsity);
%     coeff_house=ompCholesky(D_new,data_house,targetSparsity);
%     coeff_lena=ompCholesky(D_new,data_lena,targetSparsity);
   
  
        
    if(second_param==1)
        filename=strcat('./output_log/coeff_barbara_',int2str(dictsize),'.mat');
        save(filename,'coeff_barbara');
    else
        filename=strcat('./output/coeff_barbara_',int2str(dictsize),'.mat');
        save(filename,'coeff_barbara');
    end
    
%     filename=strcat('./output/coeff_boats_',int2str(dictsize),'.mat');
%     save(filename,'coeff_boats');
%     filename=strcat('./output/coeff_house_',int2str(dictsize),'.mat');
%     save(filename,'coeff_house');
%     filename=strcat('./output/coeff_lena_',int2str(dictsize),'.mat');
%     save(filename,'coeff_lena');
    'matrix saved'
    
 elseif(param==2) %change it
     
    if(second_param==1)
        filename=strcat('./output_log/D_',int2str(dictsize),'.mat');
        load(filename,'D');
    else
        filename=strcat('./output/D_',int2str(dictsize),'.mat');
        load(filename,'D');
    end
    
    D_new=ones(size(D,1),size(D,2)+1);
    D_new=D_new/p;
    D_new(:,2:end)=D;
    D_new=D;
    
    if(second_param==1)
        filename=strcat('./output_log/coeff_barbara_',int2str(dictsize),'.mat');
        load(filename,'coeff_barbara');
    else
        filename=strcat('./output/coeff_barbara_',int2str(dictsize),'.mat');
        load(filename,'coeff_barbara');
    end
    
    
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
%     for col=1:size(Y_new,2)
%         Y_new(:,col)=Y_new(:,col)+minInt(1,col);
%     end
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
    for ii=1:size(Y,1)
        for jj=1:size(Y,2)
            if Y(ii,jj)>1
                Y(ii,jj)=1;
            end
            if Y(ii,jj)<0
                Y(ii,jj)=0;
            end
        end
    end
%     Y=mat2gray(Y);
    rmse_image=(Y(:)-X_old(:)).^2;
    rmse_image=sqrt(sum(rmse_image)/(size(Y,1)*size(Y,2)));
    rmse_count=[rmse_count rmse_image];
    
    if(second_param==1)
        filename=strcat('./output_log/outpimage_barbara',int2str(numatoms),'.png');
        save_image(Y,filename,0);
    else
        filename=strcat('./output/outpimage_barbara',int2str(numatoms),'.png');
        save_image(Y,filename,0);
    end
    
%     Y_new=D_new*coeff_boats;   
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
%     Y=mat2gray(Y);
%     rmse_image=(Y(:)-X_boats(:)).^2;
%     rmse_image=sqrt(sum(rmse_image)/(size(Y,1)*size(Y,2)));
%     rmse_count=[rmse_count rmse_image];
%     filename=strcat('./output/outpimage_boats_',int2str(numatoms),'.png');
%     save_image(Y,filename,0);
%     
%     Y_new=D_new*coeff_house;   
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
%     Y=mat2gray(Y);
%     rmse_image=(Y(:)-X_house(:)).^2;
%     rmse_image=sqrt(sum(rmse_image)/(size(Y,1)*size(Y,2)));
%     rmse_count=[rmse_count rmse_image];
%     filename=strcat('./output/outpimage_house_',int2str(numatoms),'.png');
%     save_image(Y,filename,0);
%     
%     Y_new=D_new*coeff_lena;   
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
%     Y=mat2gray(Y);
%     rmse_image=(Y(:)-X_lena(:)).^2;
%     rmse_image=sqrt(sum(rmse_image)/(size(Y,1)*size(Y,2)));
%     rmse_count=[rmse_count rmse_image];
%     filename=strcat('./output/outpimage_lena_',int2str(numatoms),'.png');
%     save_image(Y,filename,0);
%     
    
    
    
    
    
    %save_image(X, filename, 0);
    
    %display_dictionary(D,p,numDisplay,size(D,2));
    %imshow(Y);
    %pause
    
    
   
end