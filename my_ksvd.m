function [D]=my_ksvd(param,data,dictsize,maxiter,p,numDisplay,sparseParam,targetSparsity)

    n=size (data,1);
    m=size (data,2);

    %%%%%%%%%%%%%%%%%%%%%%%% Farthest Point Clustering %%%%%%%%%%%%%%%%%%%%%%%
    D_initial = zeros (n, dictsize);
    D_initial(:,1)= data(:,randi(m));
    errors=nan(1,maxiter);

    dist=zeros(dictsize,m);

    for i=1:dictsize-1,
        diff=data-D_initial(:,i*ones(1,m));
        dist(i,:)=sum(diff.^2);  
        if i==1
            [~, maxind]=max(dist(1,:));
        else
            [~, maxind]=max(min(dist(1:i,:)));
        end
        D_initial(:,i+1)=data(:,maxind);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Kmeans %%%%%%%%%%%%%%%%
    [~,D] = kmeans(data',dictsize, 'start', D_initial'); 
    D=D';
    D = normc(D);
    D_kmeans=D;
    'kmeans done';
    
%     filename=strcat('./output/D_kmeans_',int2str(dictsize),'.mat');
%     save(filename,'D');
    'matrix kmeans saved'
    if(param==3)
        return;
    end
    
    %kmeansD=D;
    %kmeansdictsize=dictsize;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Displaying Initial Dictionary %%%%%%%%%%%%%%

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for it_count=1:maxiter

        if (sparseParam)
            X = sparsecodeNew(D,data,targetSparsity);
         else
            X = ompCholesky(D,data);
        end
        
% res=data(:,1)-D*X(:,1);
% res'
% res=data(:,1)-D*X1(:,1);
% res'
% 'residual'
% res=X(:,1);
% res'
% res=X1(:,1);
% res'
% pause

        
        residual = data - D * X;
        curr_error = sqrt (mean (residual(:).^2));
        errors(1,it_count)=curr_error;
        threshold = prctile (abs (X(:)), 85);
     
        for j_count = 1:dictsize
            [D,X,reset_flag] = optimize_atom(data,D,j_count,X,threshold);
            if(reset_flag==true)
                dictsize=dictsize-1;
                it_count=it_count-1;
                break;
            end
           
        end
     
    end
  
    errors=errors(1:it_count);
 
end




%helper functions

function [D,X,reset_flag]=optimize_atom(Y,D,j,X,threshold)

    data_indices = find(abs(X(j,:)) > threshold); % thresholds
    reset_flag=false;
   
    
    if(size(data_indices,2)<0.05*size(X,2)) % threshold
        temp=linspace(1,size(D,2),size(D,2));
        temp=(temp~=j).*temp;
        'scarce data item';
        D=D(:,logical(temp));
        reset_flag=true;
        return

    end

    X_j = X(j,data_indices);
    smallX = X(:,data_indices);
    Dj = D(:,j);
    [D(:,j),s,X_j] = svds(Y(:,data_indices) - D*smallX + Dj*X_j, 1);
    X(j,data_indices) = s*X_j;
    
    tempDmat=D(:,j)'*D;
    
    [maxValue,~]=max(abs(tempDmat(1:j-1)));
    
    
    if(j>1 && maxValue>0.9) %vague threshold
        temp=linspace(1,size(D,2),size(D,2));
        temp=(temp~=j).*temp;
        'correlated data item';
        D=D(:,logical(temp));
        reset_flag=true;
        return
    end


end



