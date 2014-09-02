
function a=sparsecode(D,data)
   
    
    m = size(data,2); %m is the no. of training examples
    n = size(D,1); %n is the dimension of the vector space
    k = size(D,2); %k is the no. of atoms in the dictionary
    a = zeros(k,m);
    identity=eye(n,'single');
    
    for i=1:m,
       R=data(:,i);
       phi = zeros(n,k);
       count = 1;
       flag = zeros(1,k,'uint8');
      
       while(count<=k)
       
        %'-----------------------------------'
        %flagIndexZero = find (1-flag);
        
        temp=linspace(1,k,k);
        temp=temp(~logical(flag));
        Dunused = D (:,~logical(flag));
        %Dunused = D (:,flagIndexZero);
        [ ~, maxDotAtom ] = max (abs (Dunused' * R));
        maxDotAtom = temp(maxDotAtom);
      
        flag(1,maxDotAtom)=1;
        a(maxDotAtom,i)=D(:,maxDotAtom)' * R;
        
        phi(:,count)=D(:,maxDotAtom);
        
        phi_current=phi(:,1:count);
        %imagesc (phi_current); colorbar;
        intermMatrix=phi_current'*phi_current;
        %pause
        P=phi_current*(intermMatrix \ phi_current');
        count=count+1;
        R=(identity-P)*R;
       

       end
    end
end
