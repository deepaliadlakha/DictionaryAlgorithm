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

    [D,X,errors]=my_ksvd (data,10,sigma*sqrt(numel(data)),100); %Threshold not clear
end

