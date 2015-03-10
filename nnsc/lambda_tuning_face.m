function [lambda_val,error_val1,error_val2]=lambda_tuning_face(param,store_param,second_param,dictsize,maxIter,noiseLevel)

    addpath('../common/export_fig/')
    addpath('../common/')
    addpath('../')

    close all
    rng (0);
    
    thresh=0.001;

    if(param==1)
        sizeOfImage=92*112;
        noOfFaces=35;
        countEachFace=10;
%         countEachTestFace=5;
        height = 92;
        width = 112;
    else
        sizeOfImage=192*168;
        noOfFaces=39;
        countEachFace=2;
%         countEachTestFace=3;
        height = 192;
        width = 168;
    end
    
    
    if(param==2)
        noOfImages=(noOfFaces-1)*countEachFace;
%         noOfTestImages=(noOfFaces-1)*countEachTestFace;
    else
        noOfImages=noOfFaces*countEachFace;
%         noOfTestImages=noOfFaces*countEachTestFace;
    end
    
    
    X=zeros(sizeOfImage,noOfImages);
%     testImages=zeros(sizeOfImage,noOfTestImages);
    
   
    count=1;
%     testCount=1;
    
    for i=1:noOfFaces
        
        
        if(param==1)
            dir='../face_db/att_faces/s';
        else
            dir='../../../CroppedYale_Subset/CroppedYale_Subset/';
        end
        
        dir=strcat(dir,int2str(i),'/');
        
        for j=1:countEachFace
            pic=strcat(dir,int2str(j),'.pgm');
            img=imread(pic);
            img=double(img);
            img=img/255;
            X(:,count)=img(:);
            count=count+1;
        end
        
%         for j=1:countEachTestFace
%             pic=strcat(dir,int2str(j+countEachFace),'.pgm');
%             %pic=strcat(dir,int2str(j),'.pgm');
%             img=imread(pic);
%             img=double(img);
%             img=img/255;
%             %figure, imshow(mat2gray(log(1+ abs(fftshift(fft2(img))))))
%             %pause
%             testImages(:,testCount)=img(:);
%             testCount=testCount+1;
%         end 
    end
    
    data_oldSelected=X;
    dataSelected = double(X)+randn(size(X,1),size(X,2))*noiseLevel;
    
     
    
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
    
    lambda_val=[0,0.0001,0.001,0.01,0.1];
    error_val=zeros(size(lambda_val,1),size(lambda_val,2));
    error_val2=zeros(size(lambda_val,1),size(lambda_val,2));
    
    
    for i=1:size(lambda_val,2)
        
     for j=1:5
        train_set=dataSelected(:,~logical(mask(j,:)));
        val_set=dataSelected(:,logical(mask(j,:)));
        val_set_true=data_oldSelected(:,logical(mask(j,:)));
        
        
        if(second_param==1)
            [D,~,~]=nnsc_and_log(train_set,dictsize,maxIter,lambda_val(i));
        else
            [D,~,~]=my_nnsc(train_set,dictsize,maxIter,lambda_val(i));
        end
        
        D_new=ones(size(D,1),size(D,2)+1);
        D_new=D_new/sqrt(height*width);
        D_new(:,2:end)=D;


        coeff_val=ones(size(D_new,2),size(val_set,2));

        if(second_param==1)
            S=coeff_val;
            A=D_new;
            eps=0.001;

            Sn=S;
            for cols=1:size(S,2)
%                 for ii=1:10

                ii=1;
                W=diag(ones(size(S(:,cols))));
                err=norm(val_set(:,cols)-A*W*S(:,cols),'fro')+lambda_val(i)*sum(log(abs(S(:,cols))+eps));
                lasterr=0;
                while(abs(err-lasterr)>thresh*lasterr && ii<150)
                    lasterr=err;
                    W=diag(S(:,cols)+eps);
                    An=A*W;
                    S(:,cols)=(S(:,cols).*(An'*X(:,cols)))./((An'*An)*S(:,cols)+lambda_val(i));
                    err=norm(val_set(:,cols)-A*W*S(:,cols),'fro')+lambda_val(i)*sum(log(abs(S(:,cols))+eps));
%                     norm(X-An*S,'fro')
%                     S(:,cols)
                    ii=ii+1;
                end
                
%                norm(data(:,cols)-A*W*S(:,cols),'fro')
%                 pause
                Sn(:,cols)=W*S(:,cols);
            end
            S=Sn;
            coeff_val=S;
        else
            for ii=1:100
                coeff_val=(coeff_val.*(D_new'*val_set))./((D_new'*D_new)*coeff_val+lambda_val(i));
            end
        end
        
         error_val(i)=error_val(i)+norm(val_set_true-D_new*coeff_val,'fro')/sqrt(numel(val_set));
         error_val2(i)=error_val2(i)+norm(val_set-D_new*coeff_val,'fro')/sqrt(numel(val_set));
     end
     
%       subplot(1,4,1),  display_dictionary(val_set_true(:,1:10),8,2);
%       subplot(1,4,2),  display_dictionary(val_set(:,1:10),8,2);
%       subplot(1,4,3), display_dictionary(D*coeff_val(:,1:10),8,2);
%       coeff_val;
      '----------------'
      lambda_val(i)
      error_val(i)=error_val(i)/5;
      error_val(i)
%       error_val2(i)
%       subplot(1,4,4), display_dictionary(D,8,2);
       
     
%       display_dictionary(D,8,1);
    end
    
    
    close all;
    
    [~,minIndex]=min(error_val);
    corr_lambda=lambda_val(minIndex);
    corr_lambda
    pause
%      corr_lambda=0;


    if(second_param==1)
        [D,~,~]=nnsc_and_log(dataSelected,dictsize,maxIter,corr_lambda);
    else
        [D,~,~]=my_nnsc(dataSelected,dictsize,maxIter,corr_lambda);
    end
        
    display_dictionary(D,8,5);colorbar;
    pause
    coeff_val=rand(size(D,2),size(dataNoisy,2)); 
         for dummy_c=1:100
             coeff_val=(coeff_val.*(D'*dataNoisy))./((D'*D)*coeff_val+corr_lambda);
         end
%                    imagesc(coeff_val);colormap('default');colorbar;pause
         
        Y_new=D*coeff_val;   
        Y=zeros(size(X,1),size(X,2));
        count=1;
        for i=1:row_lim,
            for j=1:col_lim,
               pCrosspMat=reshape(Y_new(:,count),p,p);
               Y(i+4,j+4)=pCrosspMat(4,4);
               count=count+1;
            end
        end
        
        'klklk'
        imshow(X_old);colorbar;pause
        norm(Y(10:end-10,10:end-10)-X_old(10:end-10,10:end-10),'fro')/sqrt(numel(X_old(10:end-10,10:end-10)))
        imshow(Y);colorbar;pause;
        
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