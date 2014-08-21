function [D,X]=ksvd_test
    clear all
    close all
    rng (0);


    identity=eye(10);
    data=zeros(64,10*100);
    result = cellfun(@(x) imnoise(identity(:,randi(10)), 'gaussian', 0, 0.005), num2cell(data, 1), 'UniformOutput', false);
    data=cell2mat(result);

    %data

    [D,X]=my_ksvd (data,9,21);
end

