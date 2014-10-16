


function a=ompCholesky(D,data,targetSparsity)

    m = size(data,2); %m is the no. of training examples
    n = size(D,1); %n is the dimension of the vector space
    k = size(D,2); %k is the no. of atoms in the dictionary
    a = zeros(k,m);
    
    for i=1:m,
       x=data(:,i);
       R=x;
       count = 1;
       I = [];
       flag=zeros(1,k);
       L=zeros(k,k);
       L(1:1)=1;
       gamma=zeros(k,1);
       alpha=D'*x;
       
       
       
      
       while(count<=targetSparsity)
       
        %'-----------------------------------'
        %flagIndexZero = find (1-flag);
        
        temp=linspace(1,k,k);
        temp=temp(~logical(flag));
        Dunused = D (:,~logical(flag));
        [ ~, maxDotAtom ] = max (abs (Dunused' * R));
        maxDotAtom = temp(maxDotAtom);
        flag(1,maxDotAtom)=1;
      
       epsilon = 1e-5;
        
        if(count>1)
            
            w=L(1:count-1,1:count-1)\(D(:,I)'*D(:,maxDotAtom));
            
            L(count,count) = sqrt (max (epsilon, 1 - w'*w));
            L(count,1:count-1)=w';
            
        end
        I=[I maxDotAtom];
        gamma(I,1)=L(1:count,1:count)'\(L(1:count,1:count)\alpha(I,1));
        R=double(x)-double(D(:,I)*gamma(I,1));
       
        
   
        count=count+1;

       end
       
       
       a(:,i)=gamma;
       
       
       
    end
end
