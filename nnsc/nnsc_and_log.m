
function [A,S,errors]=nnsc_and_log(data,dictsize,maxiter,lambda) %removed two things from the parameters
    
    thresh=0.001;
    Thresh=0.01;
    
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
%     A_kmeans=A;
    mu=0.05;
    
    eps=0.001;
    
    

    it_count=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Err=norm(data-A*S,'fro')+lambda*sum(log(abs(S(:))+eps));
%     Lasterr=0;
%     while(abs(Err-Lasterr)>Thresh*abs(Lasterr) && it_count<50)

      while(it_count<20)
        
        Err
        it_count
        
%         Lasterr=Err;
        it_count=it_count+1;
        
        '-----------------------------------';

        fro_norm=norm(X-A*S,'fro');
        curr_fro_error=0.5*fro_norm*fro_norm;
        curr_error=curr_fro_error+lambda*sum(S(:));
        
        errors(1,it_count)=curr_error;

        Aprev=A;
        
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

        setGlobalD_newt(A)
        for cols=1:size(S,2)
            W=1./(S(:,cols)+eps);
            setGlobalWt(W)
            setGlobaly(data(:,cols))

            x0=ones(size(S(:,cols)))';
            x0=[x0 ones(size(S(:,cols)))'];
            
            options = optimoptions(@fmincon,'Algorithm','interior-point');
            [x,~]=fmincon(@f,x0,[],[],[],[],-Inf,+Inf,@c,options);
            'returned'
            pause
            S(:,cols)=x(1,1:size(S(:,cols),2))';
        end
       

       

%             Sn=S;
%             for cols=1:size(S,2)
% %                 for ii=1:10
% 
%                 ii=1;
%                 W=diag(ones(size(S(:,cols))));
%                 err=norm(data(:,cols)-A*W*S(:,cols),'fro')+lambda*sum(log(abs(S(:,cols))+eps));
%                 lasterr=0;
% %                 while(abs(err-lasterr)>thresh*lasterr && ii<150)
%                   while(ii<150)
%                       
%                     lasterr=err;
%                     W=diag(S(:,cols)+eps);
%                     An=A*W;
%                     
%                     for si=1:100
%                         S(:,cols)=(S(:,cols).*(An'*X(:,cols)))./((An'*An)*S(:,cols)+lambda);
%                     end
%                     
%                     err=norm(data(:,cols)-An*S(:,cols),'fro')+lambda*sum(log(abs(S(:,cols))+eps));
%                     
%                     'nn'
%                     err
%                     ii=ii+1;
%                 end
%                 
%                 'done'
%                 pause
%                 
%                 
% %                norm(data(:,cols)-A*W*S(:,cols),'fro')
% %                 pause
%                 Sn(:,cols)=W*S(:,cols);
%             end
%             S=Sn;
            
%          end
%         S(:,1:10)
%         pause
        
        mu=1.2*mu;  
        Err=norm(data-A*S,'fro')+lambda*sum(log(abs(S(:))+eps));
     
    end  
    %subplot(1,2,1), display_dictionary(initA,p,numDisplay,size(initA,2));
    %subplot(1,2,2), display_dictionary(A,p,numDisplay,size(A,2));
      
    errors=errors(1:it_count);
  
end


function fans = f(x)
    tempW=getGlobalWt();
    fans = tempW'*x(1,size(x,2)/2+1:end)';
end

function [cans,ceq] = c(x)

    tempD=getGlobalD_newt();
    tempy=getGlobaly;
    eps=16000;
    myX=x(1,1:size(x,2)/2)';
    myU=x(1,size(x,2)/2+1:end)';
    
%     (tempy-tempD*myX)'*(tempy-tempD*myX)
    
    
    cans=(tempy-tempD*myX)'*(tempy-tempD*myX)-eps;
    cans2=myX-myU;
    cans3=-myX-myU;
    cans4=-myX;
    
    
    ceq=0;
    cans=[cans cans2' cans3' cans4'];
%     cans(1,1)
    (tempy-tempD*myX)'*(tempy-tempD*myX)
    cans(1,1)
    '-------------'
    if cans(1,1)<=0
        'paused'
        cans
        pause
    end
    
     if cans<=zeros(size(cans))
        'second paused'
        pause
    end
    
end



function setGlobalD_newt(val)
    global D_newt
    D_newt = val;
end

function r = getGlobalD_newt
    global D_newt
    r = D_newt;
end

function setGlobaly(val)
    global y
    y = val;
end

function r = getGlobaly
    global y
    r = y;
end

function setGlobalWt(val)
    global Wt
    Wt = val;
end

function r = getGlobalWt
    global Wt
    r = Wt;
end





