


function a=ompCholesky(D,data)

    m = size(data,2); %m is the no. of training examples
    n = size(D,1); %n is the dimension of the vector space
    k = size(D,2); %k is the no. of atoms in the dictionary
    a = zeros(k,m);
    optsLT.LT = true;
    optsUT.UT = true;
    
    for i=1:m,
       x=data(:,i);
       R=x;
       count = 1;
       I = zeros(1,k);
       flag=zeros(1,k);
       L=zeros(k,k);
       L(1:1)=1;
       gamma=zeros(k,1);
       
       gamma=single(gamma);
       x=single(x);
       D=single(D);
       
       alpha=D'*x;
       
       
       
      
%        while(count<=0.1*size(D,2))
        while(count<=min(1,size(D,2)))
       
        %'-----------------------------------'
        %flagIndexZero = find (1-flag);
        
%         temp=linspace(1,k,k);
%         temp=temp(~logical(flag));
        Dunused = D (:,~logical(flag));
        [ ~, maxDotAtom ] = max (abs (Dunused' * R));
        
        maxDotAtom=find((cumsum(1-flag)==maxDotAtom),1);
%         maxDotAtom = temp(maxDotAtom)
%         pause
        flag(1,maxDotAtom)=1;
      
        epsilon = 1e-5;
        
        if(count>1)
            
            w=linsolve(L(1:count-1,1:count-1),(D(:,I(1,1:count-1))'*D(:,maxDotAtom)),optsLT);
%             w=L(1:count-1,1:count-1)\(D(:,I(1,1:count-1))'*D(:,maxDotAtom));
            
            L(count,count) = sqrt (max (epsilon, 1 - w'*w));
            L(count,1:count-1)=w';
            
        end
        I(1,count)=maxDotAtom;
        %I=[I maxDotAtom];
%         gamma(I(1,1:count),1)=L(1:count,1:count)'\(L(1:count,1:count)\alpha(I(1,1:count),1));


    
        
        
        gamma(I(1,1:count),1)= linsolve(L(1:count,1:count)',linsolve(L(1:count,1:count),alpha(I(1,1:count),1),optsLT),optsUT);
        
        R=x-D(:,I(1,1:count))*gamma(I(1,1:count),1);
        count=count+1;

       end
       
       
       a(:,i)=gamma;
       
       
       
    end
end
