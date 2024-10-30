function U = prim2cons(u,v,T,rho,cv)
    [nx,ny]=size(rho);
    U=zeros(4,nx,ny);
    e=cv*T;
    U(1,:,:)=rho;
    U(2,:,:)=rho.*u;
    U(3,:,:)=rho.*v;
    U(4,:,:)=rho.*(e+1/2*(u.^2+v.^2));
end