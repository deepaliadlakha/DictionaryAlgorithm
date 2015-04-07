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
    
    [D_new,coeff]=modified_nnsc(store_param,second_param,data,dataOld,dictsize,lambda,width,height);


end
 
