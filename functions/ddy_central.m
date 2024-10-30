function dfdy = ddy_central(f,dy)

%% mine

    [nx,ny]=size(f);
    
    A=diag(0*ones(1,ny)) +...
        diag(1*ones(1,ny-1),1) +...
        diag(-1*ones(1,ny-1),-1);
    A(1,1:3)=[-3,4,-1];
    A(end,ny-2:end)=[1,-4,3];
    % A(1,end)=-1;
    % A(end,1)=1;
    A=sparse(A);
    A=A/2;

    % dfdy=permute((1/dy*A*permute(f,[2,1])),[2,1]);
    dfdy=transpose((1/dy*A*f'));
end