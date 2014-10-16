function [Y,D,coeff,rmse_image]=nnsc_main(param,dictsize,maxIter,numDisplay,p,mu,lambda)

    addpath('../common/export_fig/')
    addpath('../common/')

    %nnsc_main(1,10,5,1,8,0.1,0.1)

    close all
    rng (0);

    numatoms=dictsize;
    [X, ~] = imread('../barbara.png');
    X=single(X);
    X=X/255;
    
    X=X(1:200,1:200);

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

    %mean_data=mean(data);
    %mean_data=repmat(mean_data,size(data,1),1);
    %data_later=data;
    %data=data-mean_data;
    %var_data=var(data);
    %var_data=var_data.*(var_data>0.01); %based on data

    %dataSelected = data(:,logical(var_data));
    
    
    % mask images for colour intensity regions
    
    dataSelected=data;
    minInt=ones(1,size(data,2));
    for col=1:size(data,2)
        minInt(1,col)=min(data(:,col));
        dataSelected(:,col)=dataSelected(:,col)-minInt(1,col);

    end
   
    
   
    
 if(param==1)
     
    
    [A,S,~]=my_nnsc (dataSelected,dictsize,maxIter,p,numDisplay,mu,lambda); %Threshold not clear take a clear patch
    filename=strcat('./output/A_',int2str(dictsize),'.mat');
    save(filename,'A');
    filename=strcat('./output/S_',int2str(dictsize),'.mat');
    save(filename,'S');
    
    
    
    'matrix saved'
    
 else
    filename=strcat('./output/A_',int2str(dictsize),'.mat');
    load(filename,'A');
    filename=strcat('./output/S_',int2str(dictsize),'.mat');
    load(filename,'S');
    'matrix loaded'
 end
 
    %display_dictionary(D,p,numDisplay,size(D,2));
    
    A_new=ones(size(A,1),size(A,2)+1);
    A_new=A_new/p;
    A_new(:,2:end)=A;
    A=A_new;
    
    S_new=ones(size(S,1)+1,size(S,2));
    
    
    S_new(1,:)=minInt(1,:);
    S_new(2:end,:)=S;
    S=S_new;
    
   
    
    for i=1:100
        S=(S.*(A'*data))./((A'*A)*S+lambda);
    end
   
    B=find(S<0);
    assert(size(B,1)==0);
    
    Y_new=A*S;   
    
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
    save_image(Y, filename, 0);
    
    display_dictionary(A,p,numDisplay,size(A,2));
    %pause
    
    
   
end