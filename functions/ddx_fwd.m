function FirstDerivative=ddx_fwd(f,dx)

    % first order forward
    
    
    [nx,ny]=size(f);
    A=diag(-1*ones(1,nx)) +...
        diag(1*ones(1,nx-1),1) +...
        diag(0*ones(1,nx-2),2);
    
    % A(end,nx-1:end)=[-1,1];
    A(end,nx-2:end)=[1,-4,3];
    
    FirstDerivative=1/dx*A*f;
    
    
    % second order forward
    
        % A=diag(-3*ones(1,nx)) +...
        %     diag(4*ones(1,nx-1),1) +...
        %     diag(-1*ones(1,nx-2),2);
        % A(end,nx-2:end)=[1,-4,3];
        % A(end-1,nx-3:end)=[1,-4,3,0];
        % A=A/2;
        % A=sparse(A);


end
