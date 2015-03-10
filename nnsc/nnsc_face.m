function [D_new,coeff]=nnsc_face(param,store_param,second_param,dictsize,lambda,noiseLevel)
    
    addpath('../common/export_fig/')
    addpath('../common/')
    addpath('../')
    
    thresh=0.001;

    close all
    rng (0);
    
    if(param==1)
        sizeOfImage=92*112;
        noOfFaces=35;
        countEachFace=5;
        countEachTestFace=5;
        height = 92;
        width = 112;
    else
        sizeOfImage=192*168;
        noOfFaces=39;
        countEachFace=2;
        countEachTestFace=3;
        height = 192;
        width = 168;
    end
    
    
    if(param==2)
        noOfImages=(noOfFaces-1)*countEachFace;
        noOfTestImages=(noOfFaces-1)*countEachTestFace;
    else
        noOfImages=noOfFaces*countEachFace;
        noOfTestImages=noOfFaces*countEachTestFace;
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
 
    
%     mean_image=mean(X,1);
%     for i=1:size(X,1)
%         X(i,:)=X(i,:)-mean_image;
%     end
%     
%     std_dev=std(X,1);
%     for i=1:size(X,1)
%         X(i,:)=X(i,:)./std_dev;
%     end
%     
%     min_image=min(X);
%     min_image=min_image.*(min_image<0);
%     for i=1:size(X,1)
%         X(i,:)=X(i,:)-min_image;
%     end
    dataOld=X;
    data = double(X)+randn(size(X,1),size(X,2))*noiseLevel;
    
     
%     for i=1:size(data,2),
%        img=reshape(data(:,i),width,height);
%        filename=strcat('./output_log/noisy_face','_',int2str(i),'.png');
%        save_image(img,filename,0);
%     end
   
    
    if(store_param==1)
     
     
        if(second_param==1)
            [D,~,~]=nnsc_and_log(data,dictsize,50,lambda);
            filename=strcat('./output_log_face/D_',int2str(dictsize),'.mat');
            save(filename,'D');
        else
            [D,~,~]=my_nnsc(data,dictsize,50,lambda);
            filename=strcat('./output_face/D_',int2str(dictsize),'.mat');
            save(filename,'D');
        end

        D_new=ones(size(D,1),size(D,2)+1);
        D_new=D_new/sqrt(height*width);
        D_new(:,2:end)=D;


        coeff=ones(size(D_new,2),size(data,2));
    %     for i=1:100
    %         coeff_barbara=(coeff_barbara.*(D_new'*data_barbara))./((D_new'*D_new)*coeff_barbara+lambda);
    %     end

        if(second_param==1)
            S=coeff;
            A=D_new;
            eps=0.001;

            Sn=S;
            for cols=1:size(S,2)
                ii=1;
                W=diag(ones(size(S(:,cols))));
                err=norm(data(:,cols)-A*W*S(:,cols),'fro')+lambda*sum(log(abs(S(:,cols))+eps));
                lasterr=0;
                while(abs(err-lasterr)>thresh*lasterr && ii<150)
                    lasterr=err;
                    W=diag(S(:,cols)+eps);
                    An=A*W;
                    S(:,cols)=(S(:,cols).*(An'*data(:,cols)))./((An'*An)*S(:,cols)+lambda);
                    err=norm(data(:,cols)-A*W*S(:,cols),'fro')+lambda*sum(log(abs(S(:,cols))+eps))
                    ii=ii+1;
                end
%                 norm(data(:,cols)-A*W*S(:,cols),'fro')
%                 pause
                Sn(:,cols)=W*S(:,cols);
            end

            S=Sn;
            coeff=S;
        else
            for i=1:100
                coeff=(coeff.*(D_new'*data))./((D_new'*D_new)*coeff+lambda);
            end
        end

        if(second_param==1)
            filename=strcat('./output_log_face/coeff',int2str(dictsize),'.mat');
            save(filename,'coeff');
        else
            filename=strcat('./output_face/coeff',int2str(dictsize),'.mat');
            save(filename,'coeff');
        end
        
    elseif(store_param==2)

        if(second_param==1)
            filename=strcat('./output_log_face/D_',int2str(dictsize),'.mat');
            load(filename,'D');
        else
            filename=strcat('./output_face/D_',int2str(dictsize),'.mat');
            load(filename,'D');
        end

        D_new=ones(size(D,1),size(D,2)+1);
        D_new=D_new/sqrt(height*width);
        D_new(:,2:end)=D;

        if(second_param==1)
            filename=strcat('./output_log_face/coeff',int2str(dictsize),'.mat');
            load(filename,'coeff');
        else
            filename=strcat('./output_face/coeff',int2str(dictsize),'.mat');
            load(filename,'coeff');
        end

    end
    
    
    Y_new=D_new*coeff; 
    
    norm(dataOld-Y_new,'fro')
    dataOld(1:10,1:5)
    Y_new(1:10,1:5)
    pause
    
%     for i=1:size(Y_new,1)
%         Y_new(i,:)=Y_new(i,:)+min_image;
%     end
%     
%     for i=1:size(Y_new,1)
%         Y_new(i,:)=Y_new(i,:).*std_dev;
%     end
%     
%     for i=1:size(Y_new,1)
%         Y_new(i,:)=Y_new(i,:)+mean_image;
%     end
    
    
%     Y_new=dataOld;
%     for i=1:size(Y_new,2),
%        img=reshape(Y_new(:,i),width,height);
%        imshow(img);
%        filename=strcat('./output_log/face','_',int2str(i),'.png');
%        save_image(img,filename,0);
%     end
    
   
    
    
    
%     
%     for i=1:noOfTestImages
%         testImages(:,i)=testImages(:,i)-meanX;
%     end
%      
%     
%     
%     L=X'*X;
%     
%     
%     [v,D]=eig(L);  
%     [~, order] = sort(diag(D),'descend');
%     v = v(:,order);
%     %d
%     %pause;
%     eig_vec=X*v;
%     eig_vec=normc(eig_vec);
%     eig_vec=eig_vec(:,1:k);
%     %d=diag(d);
%     
%     
%     coeff=eig_vec'*X;
%     testCoeff=eig_vec'*testImages;
%     noOfHits=0;
%     
%     
%    
%     
%       
%     if(reconstruct==1)
%         %k =1;
%         noOfEigenFaces= 25;
%         Fourier = zeros(noOfEigenFaces,1);
%         %reconstructing X(:,k)
%         
%         'eigenfaces'
%         set(gca, 'LooseInset', get(gca,'TightInset'))
%         for i =1:noOfEigenFaces 
%             subplot(5,5,i);
%             prevSum = normc(eig_vec(:,i));%*sqrt(d(i));
%             Fourier(i) = log(1+ norm(fft(prevSum)));
%             pic=mat2gray(reshape(prevSum*255,width,height));
%             imshow(pic);
%             title(num2str(i));
%         end
%         save_image(pic,'../images/eigenfaces.png',0);
%         
%         
%         'fourier'
%         for i =1:noOfEigenFaces
%              
%             subplot(5,5,i);
%             prevSum = normc(eig_vec(:,i));
%             pic=mat2gray(log(1+ abs(fftshift(fft2(reshape(prevSum*255,width,height))))));
%             imshow(pic);
%         end
%         save_image(pic,'../images/fourier.png',0);
%         
%         'reconstruction'
%         recons_k = [2, 10, 20, 50, 75, 100, 125, 150, 175];
%         temp=zeros(size(coeff,1),1);
%         for i =1:9
%             subplot(3,3,i);
%             temp(1:recons_k(1,i),1)=coeff(1:recons_k(1,i),1)';
%             img_rec=meanX+eig_vec*temp;
%             pic=mat2gray(reshape(img_rec*255,width,height));
%             imshow(pic);
%         end
%         save_image(pic,'../images/reconstruction.png',0);
%         
%         
%     end
%    
%     for i=1:noOfTestImages
%         
%        temp=coeff;
%        for j=1:noOfImages
%           temp(:,j)=temp(:,j)-testCoeff(:,i);
%        end
%        %sum(temp(:,i))
%        %pause;
%        temp=temp.^2;
%        temp=sum(temp,1);
%        [~,ind]=min(temp);
%        
%        
%        if(floor((ind-1)/countEachFace)==floor((i-1)/countEachTestFace))
%            noOfHits=noOfHits+1;
%        end
%     end
%     recog_rate=noOfHits/noOfTestImages;
%    
%     
%    if(checkRecognition==1)
%        count=1;
%        for i=36:40,
%             new_X=zeros(sizeOfImage,50);
%             dir='../../../att_faces/s';
%             dir=strcat(dir,int2str(i),'/');
% 
%             for j=1:10
%                 pic=strcat(dir,int2str(j),'.pgm');
%                 img=imread(pic);
%                 img=img/255;
%                 new_X(:,count)=img(:);
%                 count=count+1;
%             end
%        end
%        
%        for i=1:50
%         new_X(:,i)=new_X(:,i)-meanX;
%        end
%     
%        
%        
%        
%        totalTestSet=zeros(sizeOfImage,225);
%        totalTestSet(:,1:175)=testImages;
%        totalTestSet(:,176:225)=new_X;
%        totaltestCoeff=eig_vec'*totalTestSet;
%        notRecognised=0;
%        noOfRecognised=0;
%             falsePositive =0;
%             falseNegative = 0;
%        
%        for i=1:225
%  
%             temp=coeff;
%             
%             dotProducts = zeros(size(coeff,2),1);
%             for j=1:size(coeff,2)
%                 dotProducts(j) = abs(sum(temp(:,j).*totaltestCoeff(:,i)))/(norm(temp(:,j))*norm(totaltestCoeff(:,i)));
%                  
%             end
%             [~,ind]=max(dotProducts);
%             
%             threshold = 0.7 ;
%             if(dotProducts(ind) > threshold)
%                 
%                 if(i<=175)
%                     if(floor((ind-1)/5)==floor((i-1)/5))
%                         noOfRecognised=noOfRecognised+1;
%                     else
%                         notRecognised = notRecognised+1;
%                     end
%                 else
%                     falsePositive =falsePositive+1;
%                     notRecognised = notRecognised+1;
%                 end
%                 
%             else
%                 if(i>175)
%                     noOfRecognised=noOfRecognised+1;
%                     
%                 else
%                     falseNegative =falseNegative+1;
%                     notRecognised = notRecognised+1;
%                 end
%             end
%        end
%         
%     'rate'
%     rate=noOfRecognised/(noOfRecognised+notRecognised);
%     rate
%     'False Positive'
%     falsePositive
%     'False Negative'
%     falseNegative
%     pause
%     
% 
%        

end
 
