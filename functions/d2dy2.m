function d2dy2=d2dy2(f,dy)
  [nx,ny]=size(f);
    
    A=diag(-2*ones(1,ny)) +...
        diag(1*ones(1,ny-1),1) +...
        diag(1*ones(1,ny-1),-1);
    A(1,1:4)=[2,-5,4,-1];
    A(end,ny-3:end)=[-1,4,-5,2];
    A=sparse(A);

    d2dy2=permute((1/dy^2*A*permute(f,[2,1])),[2,1]);

  
end