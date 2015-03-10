
function [A,S,errors]=nnsc_and_log(data,dictsize,maxiter,lambda) %removed two things from the parameters
    
    thresh=0.005;
    n=size(data,1);
    m=size(data,2);
    %%%%%%%%%%%%%%%%%%%%%%%% Farthest Point Clustering %%%%%%%%%%%%%%%%%%%%%%%
    D_initial = zeros (n, dictsize);
    rng(0);
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
   
    'kmeans done';
    
    
    filename=strcat('./output_log/D_kmeans_',int2str(dictsize),'.mat');
    save(filename,'D');
    'matrix kmeans saved';
    
%     X = ompCholesky(D,data,3);
%     
%     for i=1:size(D,2),
%         tempX=X(i,:);
%         a=(tempX>0).*tempX;
%         if(find(a)<size(X,2)/2)
%             D(:,i)=-1*D(:,i);
%             b=(tempX<0).*tempX;
%             X(i,:)=-1*b;
%         else
%             X(i,:)=a;
%         end
%         
%     end

    X=rand(size(D,2),size(data,2));

    for dummy_c=1:100
        X=(X.*(D'*data))./((D'*D)*X+lambda);
    end
    
%     imagesc(X);colorbar;axis on;
%     curr_error=norm(data-D*X,'fro')
%     pause

        %X is data, A is D, S is coeff
    
    S=X; %coeff   
    X=data; %data
    A=D; %dictionary
    A_kmeans=A;
    mu=0.05;
    
    eps=0.001;
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for it_count=1:maxiter

        '-----------------------------------';

        fro_norm=norm(X-A*S,'fro');
        curr_fro_error=0.5*fro_norm*fro_norm;
        curr_error=curr_fro_error+lambda*sum(S(:));
        
        errors(1,it_count)=curr_error;

        Aprev=A;
        A;
        %Sprev=S;
        
        A=A-mu*(A*S-X)*S';
        A=(A>0).*A;
        A=normc(A);
        
        ai=1;
        while(0.5*norm(X-A*S,'fro')*norm(X-A*S,'fro')>curr_fro_error && ai<200)
            mu=mu/2;
            A=Aprev-mu*(Aprev*S-X)*S';
            A=(A>0).*A;
            A=normc(A);
            ai=ai+1;
        end
        
%         for ii=1:10 
%             Sr=(S.*(A'*X))./((A'*A)*S+lambda);
%         end
%         Sr(:,1:10)
%         pause
        
%         for ci=1:10

        it_count

            Sn=S;
            for cols=1:size(S,2)
%                 for ii=1:10

                ii=1;
                W=diag(ones(size(S(:,cols))));
                err=norm(data(:,cols)-A*W*S(:,cols),'fro')+lambda*sum(log(abs(S(:,cols))+eps));
                lasterr=0;
                while(abs(err-lasterr)>thresh*lasterr && ii<150)
                    lasterr=err;
                    W=diag(S(:,cols)+eps);
                    An=A*W;
                    S(:,cols)=(S(:,cols).*(An'*X(:,cols)))./((An'*An)*S(:,cols)+lambda);
                    err=norm(data(:,cols)-A*W*S(:,cols),'fro')+lambda*sum(log(abs(S(:,cols))+eps));
%                     norm(X-An*S,'fro')
%                     S(:,cols)
                    ii=ii+1;
                end
                
%                norm(data(:,cols)-A*W*S(:,cols),'fro')
%                 pause
                Sn(:,cols)=W*S(:,cols);
            end
            S=Sn;
            
%          end
%         S(:,1:10)
%         pause
        
        mu=1.2*mu;       
     
    end  
    %subplot(1,2,1), display_dictionary(initA,p,numDisplay,size(initA,2));
    %subplot(1,2,2), display_dictionary(A,p,numDisplay,size(A,2));
      
    errors=errors(1:it_count);
  
end





