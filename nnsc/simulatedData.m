function [D_new1,D_new2,coeff1,coeff2]=simulatedData(noiseLevel,lambda)

    rng(0);
    addpath('../common/export_fig/')
    addpath('../common/')
    addpath('../')

    features=zeros(100,20);
   
    
    for i=1:10
       one_features=zeros(10,10);
       one_features(i,:)=ones(1,10);
       features(:,i)=one_features(:);
    end
    
    for i=1:10
       one_features=zeros(10,10);
       one_features(:,i)=ones(10,1);
       features(:,i+10)=one_features(:);
    end
    
    features=normc(features);
    
    n=10;
    
    
    bins=linspace(1,n,n);
    combinations_two=nchoosek(bins,2);
    combinations_one=nchoosek(bins,1);
    all_features=zeros(100,size(combinations_two,1)+size(combinations_one,1));
    all_features1=zeros(100,size(combinations_two,1)+size(combinations_one,1));
    
    for i=1:size(combinations_two,1)
        two=features(:,combinations_two(i,:));
        all_features(:,i)=two(:,1)+two(:,2);
        
        two=features(:,10+combinations_two(i,:));
        all_features1(:,i)=two(:,1)+two(:,2);
        
    end
    
    for i=1:size(combinations_one,1)
        one=features(:,combinations_one(i,:));
        all_features(:,size(combinations_two,1)+i)=one(:,1);
        
        one=features(:,10+combinations_one(i,:));
        all_features1(:,size(combinations_two,1)+i)=one(:,1);
    end
    
    all_features=[all_features all_features1];
    all_features=normc(all_features);
    
    
    n=110;
    bins=linspace(1,n,n);
    combinations_two=nchoosek(bins,2);
    combinations_one=nchoosek(bins,1);
    data=zeros(100,size(combinations_two,1)+size(combinations_one,1));
    
    for i=1:size(combinations_two,1)
        two=all_features(:,combinations_two(i,:));
        data(:,i)=two(:,1)+two(:,2);
    end
    
    for i=1:size(combinations_one,1)
        one=all_features(:,combinations_one(i,:));
        data(:,size(combinations_two,1)+i)=one(:,1);
    end
    
    data=normc(data);
    data=datasample(data,500,2,'Replace',false);
    
    data_later=data;
    data = double(data)+randn(size(data,1),size(data,2))*noiseLevel;
    
    
    data=double(data)*255;
    data_later=double(data)*255;
    
    [D_new1,coeff1]=modified_nnsc(1,1,data,data_later,100,lambda,10,10);
    'of'
    
%     display_dictionary(D_new1,10,10);
%     pause
    [D_new2,coeff2]=modified_nnsc(1,2,data,data_later,100,lambda,10,10);
    'of'
    
    

    
    
    
    
%     dataSelected=data;
%     dataSelected = double(dataSelected)+randn(size(dataSelected,1),size(dataSelected,2))*1;
%     dataSelected=dataSelected.*(dataSelected>0);
%     dataSelected=normc(dataSelected);
%     
%     display_dictionary(dataSelected(:,1:5),3,5);
%     pause
    
%    dataSelected = datasample(data,size(data,2),2,'Replace',false);
%     display_dictionary(dataSelected,3,1);
%     pause
    
    
end

    
    
    