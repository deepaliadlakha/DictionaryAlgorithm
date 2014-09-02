for i=1:400000
    'hgjgjgjh'
    k=10;
    n=64;
    flag=zeros(1,k);
    flag(1,5)=1;
    flag(1,8)=1;
    phi = zeros(n,k);
     temp=linspace(1,k,k);
     temp=temp(~logical(flag));
    A=magic(10);
    A=A(~logical(flag));
    B=magic(10);
    
    phi=arrayfun(@(x) 3,phi);
    X=phi'*phi;
end