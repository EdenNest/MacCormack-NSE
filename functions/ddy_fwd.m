function dfdy = ddy_fwd(f,dy)

% FIRST ORDER
[nx,ny]     = size(f);
  A=diag(-1*ones(1,ny)) +...
        diag(1*ones(1,ny-1),1) +...
        diag(0*ones(1,ny-2),2);
    A(end,ny-2:end)=[1,-4,3];


A=sparse(A);
dfdy=transpose((1/dy*A*f'));


%SECOND ORDER
    % A=diag(-3*ones(1,ny)) +...
    %     diag(4*ones(1,ny-1),1) +...
    %     diag(-1*ones(1,ny-2),2);
    % A(end,ny-2:end)=[1,-4,3];
    % A(end-1,ny-3:end)=[1,-4,3,0];

end