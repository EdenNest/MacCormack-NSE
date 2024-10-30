function FirstDerivative=ddx_central(f,dx)
    [nx,ny]=size(f);
    
    A=diag(0*ones(1,nx)) +...
        diag(1*ones(1,nx-1),1) +...
        diag(-1*ones(1,nx-1),-1);
    A(1,1:3)=[-3,4,-1];
    A(end,nx-2:end)=[1,-4,3];
    % A(1,end)=-1;
    % A(end,1)=1;
    A=sparse(A);
    f=double(f);
    FirstDerivative=1/2/dx*A*f;


end