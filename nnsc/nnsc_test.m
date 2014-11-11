function [error_val1]=nnsc_test(maxIter,numDisplay)
%     rng(0);

    rng(0);
    addpath('./common/export_fig/')
    addpath('./common/')
    addpath('../')

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
    
    %display_dictionary(features,3,5,size(features,2));
    %pause
    %filename=strcat('./output_nnsc/original_pic.png');
    %save_image(dictionaryPic, filename, 0);
    
    n=10;
    
%     features = datasample(features,3,2,'Replace',false);
%     n=3;
    
    p=3;
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
    
    data_later=data;
    data=normc(data);
    dataSelected=data;
%    dataSelected = datasample(data,size(data,2),2,'Replace',false);
%     display_dictionary(dataSelected,3,1);
%     pause
    
    
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
    lambda_val=[0,1,10,100,1000,10000,10000];
    error_val=zeros(size(lambda_val,1),size(lambda_val,2));
    
    dictsize=45;
    for i=1:size(lambda_val,2)
     for j=1:1
%         train_set=dataSelected(:,~logical(mask(j,:)));
%         val_set=dataSelected(:,logical(mask(j,:)));
        
        train_set=dataSelected(:,size(combinations_three,1)+1:size(combinations_two,1)+size(combinations_three,1));
%         val_set=data_later(:,size(combinations_two,1)+size(combinations_three,1)+1:end);
        val_set=data_later(:,1:size(combinations_three,1));
        
        
        [D,~,~]=my_nnsc (train_set,dictsize,maxIter,p,numDisplay,lambda_val(i));
%         coeff_val=ompCholesky(D,train_set,3);
%         coeff_val
        'kklll'
        
       
        coeff_val=ompCholesky(val_set,D,3);
        coeff_val
        
%         rng(0);
%          coeff_val=rand(size(D,2),size(val_set,2));
        coeff_val=ompCholesky(D,val_set,3);
% %         coeff_val
%        
%         for dummy_c=1:100
%             coeff_val=(coeff_val.*(D'*val_set))./((D'*D)*coeff_val+lambda_val(i));
%         end
% %         coeff_val
% %         pause
         error_val(i)=error_val(i)+norm(val_set-D*coeff_val,'fro')/numel(val_set);
%         
       
     end
%       display_dictionary(D,3,1,size(D,2));
%       pause

      subplot(1,3,1),  display_dictionary(val_set,3,5);
      subplot(1,3,2), display_dictionary(D*coeff_val,3,5);
      coeff_val;
      error_val(i)
      subplot(1,3,3), display_dictionary(D,3,5);
      pause

%       subplot(1,2,1),  display_dictionary(features,3,1);
%       error_val(i)
%       subplot(1,2,2), display_dictionary(D,3,1);
%       pause
    
    
    end
    
    error_val1=error_val;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
%     error_val=zeros(size(lambda_val,1),size(lambda_val,2));
    
%     dictsize=5;
%     for i=1:size(lambda_val,2)
%      for j=1:5
%         train_set=dataSelected(:,~logical(mask(j,:)));
%         val_set=dataSelected(:,logical(mask(j,:)));
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
% %      display_dictionary(D,8,1,size(D,2));
% %      pause
%     end
%     
%     error_val2=error_val;
   
    
    
    

    
    
    