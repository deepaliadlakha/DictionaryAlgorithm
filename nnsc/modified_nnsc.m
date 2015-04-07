function [D_new,coeff]=modified_nnsc(store_param,second_param,data,dataOld,dictsize,lambda,width,height)
    

    addpath('../common/export_fig/')
    addpath('../common/')
    addpath('../')
    
    thresh=0.001;

    close all
    rng (0);
    
    if(store_param==1)
        if(second_param==1)
            [D,~,~]=nnsc_and_log(data,dictsize,50,lambda);
            filename=strcat('./output_log_face/D_',int2str(dictsize),'.mat');
            save(filename,'D');
        else
            [D,~,~]=my_nnsc(data,dictsize,50,lambda);
            filename=strcat('./output_face/D_',int2str(dictsize),'.mat');
            save(filename,'D');
        end

        
        D_new=ones(size(D,1),size(D,2)+1);
        D_new=D_new/sqrt(height*width);
%         D_new(:,2:end)=D;
        D_new=D;
        setGlobalD_newt(D_new)
        

        coeff=ones(size(D_new,2),size(data,2));
  
        if(second_param==1)
            S=coeff;
%             A=D_new;
            eps=0.001;

%             Sn=S;
            for cols=1:size(S,2)
                global W;
                W=1./(S(:,cols)+eps);
                setGlobaly(data(:,cols))
                
                x0=[ones(size(S(:,cols)));ones(size(S(:,cols)))];
                options = optimoptions(@fmincon,'Algorithm','interior-point');
                [x,~]=fmincon(@f,x0,@c,options);
                S(:,cols)=x(1,:);
                
%                 err=norm(data(:,cols)-A*W*S(:,cols),'fro')+lambda*sum(log(abs(S(:,cols))+eps));
%                 lasterr=0;
%                 while(abs(err-lasterr)>thresh*abs(lasterr) && ii<150)
% %                     ii
% %                     abs(err-lasterr)/abs(lasterr)
%                     lasterr=err;
%                     W=diag(S(:,cols)+eps);
%                     An=A*W;
%                     
%                     for si=1:100
%                         S(:,cols)=(S(:,cols).*(An'*data(:,cols)))./((An'*An)*S(:,cols)+lambda);
%                     end
%                     
%                     err=norm(data(:,cols)-A*W*S(:,cols),'fro')+lambda*sum(log(abs(S(:,cols))+eps));
%                     ii=ii+1;
%                 end
% %                 norm(data(:,cols)-A*W*S(:,cols),'fro')
% %                 pause
%                 Sn(:,cols)=W*S(:,cols);
            end

%             S=Sn;
            coeff=S;
        else
            for i=1:100
                coeff=(coeff.*(D_new'*data))./((D_new'*D_new)*coeff+lambda);
            end
        end

        if(second_param==1)
            filename=strcat('./output_log_face/coeff',int2str(dictsize),'.mat');
            save(filename,'coeff');
        else
            filename=strcat('./output_face/coeff',int2str(dictsize),'.mat');
            save(filename,'coeff');
        end
        
    elseif(store_param==2)

        if(second_param==1)
            filename=strcat('./output_log_face/D_',int2str(dictsize),'.mat');
            load(filename,'D');
        else
            filename=strcat('./output_face/D_',int2str(dictsize),'.mat');
            load(filename,'D');
        end

        D_new=ones(size(D,1),size(D,2)+1);
        D_new=D_new/sqrt(height*width);
%         D_new(:,2:end)=D;

        D_new=D;

        if(second_param==1)
            filename=strcat('./output_log_face/coeff',int2str(dictsize),'.mat');
            load(filename,'coeff');
        else
            filename=strcat('./output_face/coeff',int2str(dictsize),'.mat');
            load(filename,'coeff');
        end

    end
    
    
    Y_new=D_new*coeff; 
    
    norm(dataOld-Y_new,'fro')
    dataOld(1:10,1:5)
    Y_new(1:10,1:5)
    pause

    
    
%     Y_new=dataOld;
%     for i=1:size(Y_new,2),
%        img=reshape(Y_new(:,i),width,height);
%        imshow(img);
%        filename=strcat('./output_log/face','_',int2str(i),'.png');
%        save_image(img,filename,0);
%     end
    
   
    


end

function fans = f(x)
    fans = W'*x(2,:);
end

function cans = c(x)
    tempD=getGlobalD_newt();
    tempy=getGlobaly;
    eps=0.001;
    cans=(tempy-tempD*x(1,:))'*(y-tempD*x(1,:))-eps;
    
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


 
