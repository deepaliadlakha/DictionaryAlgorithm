function [D,X,errors]=ksvd_test
    clear all
    close all
    rng (0);

    sigma=0.1;
    identity=eye(10);
    temp=repmat(linspace(1,10,10),100,1);
    data=temp(:);
    result = cellfun(@(x) imnoise(identity(:,x), 'gaussian', 0, sigma*sigma), num2cell(data, 1), 'UniformOutput', false);
    data=cell2mat(result);

    %data
    %imagesc (data); colorbar; pause, close
    
    
    mean_data=mean(data);
    mean_data=repmat(mean_data,size(data,1),1);
    data=data-mean_data;
    imagesc (data); colorbar; pause, close
    
    var_data=var(data);
    threshold = prctile (var_data, 10);
    %imagesc (var_data); colorbar; pause, close
    var_data=var_data.*(var_data>threshold);
    imagesc (var_data); colorbar; pause, close
    
    

    [D,X,errors]=my_ksvd (data(:,logical(var_data)),10,sigma*sqrt(numel(data)),100); %Threshold not clear
end

