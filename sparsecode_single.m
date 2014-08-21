
function a=sparsecode_single(data_item,D)
   
    

    n = size(D,1); %n is the dimension of the vector space
    k = size(D,2); %k is the no. of atoms in the dictionary
    a = zeros(k,1);
    identity=eye(n,'single');


    R=data_item;
    phi = zeros(n,k);
    count = 1;
    flag = zeros(1,k,'uint8');

    while(count<=k)

        flagIndexZero = find (1-flag);
        %Dunused = D (:,~logical(flag));
        Dunused = D (:,flagIndexZero);
        [ ~, maxDotAtom ] = max (abs (Dunused' * R));
        maxDotAtom = flagIndexZero (maxDotAtom);

        flag(1,maxDotAtom)=1;
        a(maxDotAtom)=D(:,maxDotAtom)' * R;

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
