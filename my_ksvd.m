function [D,X]=my_ksvd(data,dictsize,thresh_norm)

% Initialising dictionary %

n=size (data,1);
m=size (data,2);
D_initial = zeros (n, dictsize);
D_initial(:,1)= data(:,randi(m));

dist=zeros(dictsize,m);

for i=1:dictsize-1,
    diff=data-D_initial(:,i*ones(1,m));
    dist(i,:)=sum(diff.^2);  
    [~, maxind]=max(min(dist(1:i,:)));
    D_initial(:,i+1)=data(:,maxind);
   
    
end

[~,D] = kmeans(data',dictsize, 'start', D_initial');
D=D';

%D = eye (size (data,1), dictsize);

% normalize the dictionary %
D = normc(D);


% main loop %
while true
    
    X = sparsecode(D,data);
    
    %new_data=num2cell(data,1);
    %X = cellfun(@(x) sparsecode_single(x,D), new_data, 'UniformOutput', false);
   
    for j = 1:dictsize
    [D(:,j),X_j,data_indices] = optimize_atom(data,D,j,X);
    X(j,data_indices) = X_j;
    end
    error=norm(data-D*X,'fro')
    if error<=thresh_norm
        break;
    end

end
    
X = sparsecode(D,data);
    
end




%helper functions

function [atom,X_j,data_indices]=optimize_atom(Y,D,j,X)

data_indices = find(X(j,:) > 1e-6 | X(j,:) < -1e-6);

X_j = X(j,data_indices);
smallX = X(:,data_indices);
Dj = D(:,j);
[atom,s,X_j] = svds(Y(:,data_indices) - D*smallX + Dj*X_j, 1);
X_j = s*X_j;

end

%{
function Y = colnorms_squared(X)

% compute in blocks to conserve memory
Y = zeros(1,size(X,2));
blocksize = 2000;
for i = 1:blocksize:size(X,2)
  blockids = i : min(i+blocksize-1,size(X,2));
  Y(blockids) = sum(X(:,blockids).^2);
end

end
%}
