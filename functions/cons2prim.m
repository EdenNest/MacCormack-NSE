function [u,v,T,p,rho,e,Et] = cons2prim(U,R,cv)
  
    rho=squeeze(U(1,:,:));
    
    u=squeeze(U(2,:,:))./rho;
    v=squeeze(U(3,:,:))./rho;
    Et=squeeze(U(4,:,:));
    e= (Et./rho) - 1/2.*(u.^2+v.^2);
  
    T=e./cv;
    p=rho.*R.*T;

end