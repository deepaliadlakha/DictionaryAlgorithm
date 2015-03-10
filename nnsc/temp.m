% S=coeff_barbara;
rng(0);
A=D;
eps=0.00001;
lambda=0.1;
data=randi(255,64,1);
Sn=ones(size(A,2),1);

for ci=1:10
    
    W=diag(1/(Sn+eps));
    An=A*W;
    An=A;
    for ii=1:10 
        Sn=(Sn.*(An'*data))./((An'*An)*Sn+lambda);
    end

   
end
