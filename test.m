function test(param,mu,lambda,max_iter,dictsize,numdisplay)
    rng(0);
    
    addpath('./nnsc/')
    addpath('./common/')

    
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
    bins=linspace(1,n,n);
    combinations_three=nchoosek(bins,3);
    combinations_two=nchoosek(bins,2);
    combinations_one=nchoosek(bins,1);
    data=zeros(9,size(combinations_three,1)+size(combinations_two,1)+size(combinations_one,1));
    
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
    
    
    data=normc(data);
    
    if(param==1)
        [A,~,~,A_kmeans]=my_nnsc (data,dictsize,max_iter,3,5,mu,lambda); %Threshold not clear take a clear patch
    else
        [A,~,~,A_kmeans]=my_ksvd(data,dictsize,max_iter,3,numdisplay);
    end
    
    subplot(1,3,1), display_dictionary(features,3,5,size(features,2));
    subplot(1,3,2), display_dictionary(A,3,numdisplay,size(A,2));
    subplot(1,3,3), display_dictionary(A_kmeans,3,numdisplay,size(A_kmeans,2));
    
    pause
    %filename=strcat('./output_nnsc/dict_',int2str(100*mu),'_',int2str(100*lambda),'.png');
    %save_image(dictionaryPic, filename, 0);
   
    
    
    

    
    
    