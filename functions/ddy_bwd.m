
function dfdy=ddy_bwd(f,dy)
[nx,ny]=size(f);
% first order backward
A=diag(1*ones(1,ny)) +...
    diag(-1*ones(1,ny-1),-1) +...
    diag(0*ones(1,ny-2),-2);
A(1,1:2)=[-1,1];
A=sparse(A);
f=double(f);

% dfdy=permute((1/dy*A*permute(f,[2,1])),[2,1]);
 dfdy=transpose((1/dy*A*f'));
 
% second order backward
% A=diag(3*ones(1,ny)) +...
%     diag(-4*ones(1,ny-1),-1) +...
%     diag(1*ones(1,ny-2),-2);
% A(1,1:3)=[-3,4,-1];
% A(2,1:4)=[0,-3,4,-1];
% A=A/2;



%% oliver
    % determine field size
    % [nx,ny]     = size(f);
    % 
    % % allocate return field
    % dfdy        = zeros(nx,ny);
    % 
    % % backward difference
    % for i=1:nx
    %     for j=2:ny
    %         dfdy(i,j) = (f(i,j)-f(i,j-1))/dy;
    %     end
    % end
    % 
    % % forward difference for first point
    % j = 1;
    % for i=1:nx
    %     dfdy(i,j) = (f(i,j+1)-f(i,j))/dy;
    % end
    % 



end




