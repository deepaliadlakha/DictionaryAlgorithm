function [D,X,errors]=my_ksvd(data,dictsize,thresh_norm,maxiter)

% Initialising dictionary %

n=size (data,1);
m=size (data,2);
D_initial = zeros (n, dictsize);
D_initial(:,1)= data(:,randi(m));
errors=nan(1,maxiter);


dist=zeros(dictsize,m);

for i=1:dictsize-1,
    diff=data-D_initial(:,i*ones(1,m));
    dist(i,:)=sum(diff.^2);  
    [~, maxind]=max(min(dist(1:i,:)));
    D_initial(:,i+1)=data(:,maxind);
end


[~,D] = kmeans(data',dictsize, 'start', D_initial'); %There is a problem what if data has less than k clusters
%size (D)
%D = D + 0.2 * randn (size (D));
%[~,D] = kmeans(data',dictsize);
D=D';
%imagesc(D);colorbar;
%pause, close



%D = eye (size (data,1), dictsize);

% normalize the dictionary %
D = normc(D);


% main loop %
for it_count=1:maxiter
    
    '-----------------------------------'
    X = sparsecode(D,data);
    
    %new_data=num2cell(data,1);
    %X = cellfun(@(x) sparsecode_single(x,D), new_data, 'UniformOutput', false);
    
    j=1;
    threshold = prctile (abs (X(:)), 90)
    for j_count = 1:dictsize
    [D,X,j] = optimize_atom(data,D,j,X,threshold);
    %X(j,data_indices) = X_j;
    end
    curr_error=norm(data-D*X,'fro')
    errors(1,it_count)=curr_error;
    %imagesc (D); colorbar; pause, close
    if curr_error<=thresh_norm
        break;
    end

end

errors=errors(1:it_count);
%imagesc(D);colorbar;
X = sparsecode(D,data);
%pause   
end




%helper functions

function [D,X,j_new]=optimize_atom(Y,D,j,X,threshold)

data_indices = find(abs(X(j,:)) > threshold); % thresholds
if(size(data_indices,2)<0.05*size(X,2)) % threshold
    temp=linspace(1,size(D,2),size(D,2));
    temp=(temp~=j).*temp;
    'scarce data item'
    %X(j,:)
    size(D)
    size(X)
    pause
    D=D(:,logical(temp));
    X=X(logical(temp),:);
    j_new=j;
    return
    
end

j_new=j+1;
X_j = X(j,data_indices);
smallX = X(:,data_indices);
Dj = D(:,j);
[D(:,j),s,X_j] = svds(Y(:,data_indices) - D*smallX + Dj*X_j, 1);
X(j,data_indices) = s*X_j;

new_D=num2cell(D,1);
tempD = cellfun(@(x) corr(x,D(:,j)), new_D, 'UniformOutput', false);
tempDmat=cell2mat(tempD);
if(max(tempDmat(1:j-1))>0.9) %vague threshold
    temp=linspace(1,size(D,2),size(D,2));
    temp=(temp~=j).*temp;
    'correlated data item'
    size(D)
    size(X)
    pause
    D=D(:,logical(temp));
    X=X(logical(temp),:);
    j_new=j;
    return
end

end

