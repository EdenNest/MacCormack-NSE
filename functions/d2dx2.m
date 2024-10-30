function d2fdx2 = d2dx2(f,dx)

    [nx,ny]=size(f);

    A=diag(-2*ones(1,nx)) +...
        diag(1*ones(1,nx-1),1) +...
        diag(1*ones(1,nx-1),-1);
    A(1,1:4)=[2,-5,4,-1];
    A(end,nx-3:end)=[-1,4,-5,2];
    A=sparse(A);


    d2fdx2=1/dx^2*A*f;


end