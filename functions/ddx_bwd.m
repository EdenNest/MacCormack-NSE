
function FirstDerivative=ddx_bwd(f,dx)
[nx,ny]=size(f);

% first order backward
A=diag(1*ones(1,nx)) +...
    diag(-1*ones(1,nx-1),-1) +...
    diag(0*ones(1,nx-2),-2);
A(1,1:3)=[-3,4,-1];
A=sparse(A);
f=double(f);

FirstDerivative=1/dx*A*f;

% second order backward
% A=diag(3*ones(1,nx)) +...
%     diag(-4*ones(1,nx-1),-1) +...
%     diag(1*ones(1,nx-2),-2);
% A(1,1:3)=[-3,4,-1];
% A(2,1:4)=[0,-3,4,-1];
% A=A/2;



end
