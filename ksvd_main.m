function [Y,D,coeff,rmse_image]=ksvd_main(param,dictsize,maxIter,numDisplay,p,sparseParam,targetSparsity)

    close all
    rng (0);

    numatoms=dictsize;
    [X, ~] = imread('barbara.png');
    X=mat2gray(X);
    %X=X/255;
    
%     X=X(1:200,1:200);

    row_lim=size(X,1);
    col_lim=size(X,2);

    row_lim=row_lim-p+1;
    col_lim=col_lim-p+1;

    data=zeros(p*p,row_lim*col_lim);

    count=1;
    for i=1:row_lim,
        for j=1:col_lim,
           window = X(i:i+p-1,j:j+p-1);
           temp=window(:);
           data(:,count)= temp;
           count=count+1;
        end
    end

    mean_data=mean(data);
    mean_data=repmat(mean_data,size(data,1),1);
    data_later=data;
    data=data-mean_data;
    var_data=var(data);
    var_data=var_data.*(var_data>0.01); %based on data

    dataSelected = data(:,logical(var_data));
 if(param==1)
     
    
    [D]=my_ksvd (param,dataSelected,dictsize,maxIter,p,numDisplay,sparseParam,targetSparsity); %Threshold not clear take a clear patch
    filename=strcat('./output/D_',int2str(dictsize),'.mat');
    save(filename,'D');
    D_new=ones(size(D,1),size(D,2)+1);
    D_new=D_new/p;
    D_new(:,2:end)=D;
    
    if (sparseParam)
        coeff=sparsecodeNew(D_new,data_later,targetSparsity);
    else
        coeff=ompCholesky(D_new,data_later,targetSparsity);
    end
        
    
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
    [D]=my_ksvd (param,dataSelected,dictsize,maxIter,p,numDisplay,sparseParam,targetSparsity); %Threshold not clear take a clear patch
    filename=strcat('./output/D_kmeans_',int2str(dictsize),'.mat');
    save(filename,'D');
    D_new=ones(size(D,1),size(D,2)+1);
    D_new=D_new/p;
    D_new(:,2:end)=D;
    
    if (sparseParam)
        coeff=sparsecodeNew(D_new,data_later,targetSparsity);
    else
        coeff=ompCholesky(D_new,data_later,targetSparsity);
    end
    
    filename=strcat('./output/coeff_kmeans_',int2str(dictsize),'.mat');
    save(filename,'coeff');
    'matrix kmeans saved'
     
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
    
    
    Y=zeros(size(X,1),size(X,2));
    count=1;
    for i=1:row_lim,
        for j=1:col_lim,
            
           
           pCrosspMat=reshape(Y_new(:,count),p,p);
           Y(i+4,j+4)=pCrosspMat(4,4);
           count=count+1;
        end
    end
    
    rmse_image=(Y(:)-X(:)).^2;
    rmse_image=sqrt(sum(rmse_image)/(size(Y,1)*size(Y,2)));
    filename=strcat('./output/outpimage_',int2str(numatoms),'.png');
    
    imwrite(Y,filename);
    %save_image(X, filename, 0);
    
    %display_dictionary(D,p,numDisplay,size(D,2));
    %imshow(Y);
    %pause
    
    
   
end