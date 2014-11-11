
function a=sparsecodeNew(D,data,targetSparsity)

%inverse thing changed to pinv
    
    m = size(data,2); %m is the no. of training examples
    n = size(D,1); %n is the dimension of the vector space
    k = size(D,2); %k is the no. of atoms in the dictionary
    a = zeros(k,m);
    identity=eye(n,'single');
    
    for i=1:m,
       X=data(:,i);
       R=data(:,i);
       phi = zeros(n,k);
       count = 1;
       flag = zeros(1,k,'uint8');
       gamma=zeros(k,1);
       I=zeros(1,k);
      
       %while(count<=min(k/3,15))
       while(count<=min(targetSparsity,k))
       
        %'-----------------------------------'
        %flagIndexZero = find (1-flag);
        
        temp=linspace(1,k,k);
        temp=temp(~logical(flag));
        Dunused = D (:,~logical(flag));
        %Dunused = D (:,flagIndexZero);
        [ ~, maxDotAtom ] = max (abs (Dunused' * R));
        maxDotAtom = temp(maxDotAtom);
      
        I(1,count)=maxDotAtom;
        flag(1,maxDotAtom)=1;
        %a(maxDotAtom,i)=D(:,maxDotAtom)' * R;
        
        phi(:,count)=D(:,maxDotAtom);
        
        phi_current=phi(:,1:count);
        %imagesc (phi_current); colorbar;
        intermMatrix=phi_current'*phi_current;
        %pause
        gamma(I(1,1:count),1)=(intermMatrix \ phi_current')*X;
        P=phi_current*gamma(I(1,1:count),1);
        %P=phi_current*pinv(intermMatrix)*phi_current';
        count=count+1;
        R=X-P;
       

       end
       a(:,i)=gamma;
    end
end
