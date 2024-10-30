% Eden - Midterm - 290C

clc;clear;close all
addpath("functions")
global iteration Adiabatic SchlierenDensity
   

%% USER INPUT FOR VISUALIZATION



    % Plot every # snapshot :
         % Insert your desired gap between animation snapshots
         % 0 = plotting only the last snapshot
        Animation_plot_gap = 15;


    % Plot convergence animation (1 = active, 0 = inactive)
        Plot_convergence = 1; 

    % Schlieren plot for Density ( 1 = active, 0 = in-active )
        SchlierenDensity = 1;

    % Adiabatic Choice   ( 0 = constant wall temp | 1 = adiabatic wall)
        Adiabatic = 0;  

    % Activate CFL condition  (1 = active, 0 = in-active (constant dt))  
        CFLdt = 1;

    % Save results afterwards
        save_result = 0;
    
    % Iteration 
        iteration = 1500 ; 
    
    % Mach Number
        Mach = 4;
     

%%  CONSTANT PARAMETERS
    global dx dy 
    global p0 T0 u_inf 
    global R Pr cp cv S1 mu0 gamma


    M = Mach;
    nx = 75;
    ny = 80;
    L = 1e-5;
    H = 8e-6;
    dt_cnst = 2.35e-11; 
    
    % Adiabatic doesn't converge with CFL. Enforce constant dt
    if Adiabatic==1
          CFLdt = 0;
    end
  

   % standard sea level
    p0 = 101300;
    rho0 = 1.225; 
    R = 287;
    cp = 1005;
    cv = 718;
    gamma = 1.4;
    Pr = 0.71;
    mu0 = 1.735e-5; 
    S1 = 110.4;
    T0 = 288.15;
    u_inf = M*sqrt(gamma*R*T0);
  

    

%% DOMAIN

    x = linspace(0,L,nx);
    y = linspace(0,H,ny);
    [X,Y] = ndgrid(x,y);
    dx = x(2) - x(1);
    dy = y(2) - y(1);


%% INITIAL CONDITION
    
    % 3D matrices for Primitive u v p T With 3rd dimention as time

    u = ones(nx,ny,1) .* u_inf; v = zeros(nx,ny,1);
    T = ones(nx,ny,1) .* T0;    p = ones(nx,ny,1) .* p0;
    rho = ones(nx,ny,1).*rho0;  e = cv*T;
  

%% SOLVER STARTER

    % Starting the primitive variables inside the loop:
        u_n = squeeze(u); v_n = squeeze(v);
        T_n = squeeze(T); p_n = squeeze(p);
        rho_n = squeeze(rho); U_n = prim2cons(u_n,v_n,T_n,rho_n,cv);

    
    % plotting related parameters 
        time = 0;
        sum_u = sum(u_n,'all'); % convergence variable 1
        absdudt = sum(abs(u_n)/dt_cnst,'all');  % convergence variable 2
        dtarray = dt_cnst; % array of every dt calculated in each iteration
        





 %% SOLVER LOOP %%%%%%%%%%%%


    for n=1:iteration

      fprintf('computing iteration %d / %d \n' , n, iteration)



        % calculating time step     - with CFL or not 
            if CFLdt == 1; dt = CFL(u_n,v_n,T_n,rho_n);
            else;          dt = dt_cnst;
            end
            


       %%%%%% PREDICTOR %%%%%%%%
           
            % FD biasing 
            [E,F] = GenerateBiasedEF( 'fwd' , u_n , v_n , T_n ,p_n, rho_n);
        
            % main predictor time march
            U_bar = predictor(U_n,E,F,dt) ;

            % calculating primitives from U bar
            [u_bar,v_bar,T_bar,p_bar,rho_bar,~,~] = ...
                cons2prim(U_bar,R,cv);

            % Imposing Boundary Conditions on Primitives (also updates rho)
            [u_bar,v_bar,T_bar,p_bar,rho_bar] = ...
                imposeBC(u_bar,v_bar,T_bar,p_bar,rho_bar);
         
            % Update U bar
            U_bar = prim2cons(u_bar,v_bar,T_bar,rho_bar,cv);


            
        %%%%% CORRECTOR %%%%%%

            % FD biasing
                [E,F] = GenerateBiasedEF('bwd',u_bar,v_bar,T_bar,p_bar,rho_bar);

            % main corrector time march
                U_n = corrector(U_bar,U_n,E,F,dt);

            % Calculating primitives from U
                [u_n,v_n,T_n,p_n,rho_n,~,~] = cons2prim(U_n,R,cv);

            % Imposing BC
                [u_n,v_n,T_n,p_n,rho_n] = imposeBC(u_n,v_n,T_n,p_n,rho_n);
            
            % Update U
                U_n = prim2cons(u_n,v_n,T_n,rho_n,cv);
            
            % getting the final e
                Et_n=squeeze(U_n(4,:,:));
                e_n=(Et_n./rho_n)-1/2*(u_n.^2+v_n.^2);
            

            % Saving the calculated snapshots
                u = cat(3,u,u_n); v = cat(3,v,v_n);
                T = cat(3,T,T_n); p = cat(3,p,p_n);
                e = cat(3,e,e_n); rho = cat(3,rho,rho_n);
                
         
                
        %%%% Visualization %%%%%%
                if Plot_convergence ~= 0 %&& mod(n-1,10)==0 
                    figure(2)
                    ConvergencePlot(sum_u,absdudt,dtarray,n)
                   
                end
                if Animation_plot_gap ~= 0 && mod(n-2,Animation_plot_gap)==0 
                    
                    fig = figure(1);
                    if n == 2;set(fig ,'WindowState','maximized');end
                    visualize(n,X,Y,u,v,T,p,e,rho,time(end))
                   
                end
            
               % Updating convergence variables and time arrays for plotting
                 sum_u = cat(2,sum_u,sum(u_n,'all'));
                 absdudt = cat(2,absdudt,sum(abs(u(:,:,end)-u(:,:,end-1))/dt,'all'));
                 time = cat( 2 , time , time(end)+dt) ;
                 dtarray = cat( 2 , dtarray , dt) ;
      
    end


%%

    %%%%%%% End of calculation loop %%%%%%

  
% Plotting only the last snapshot if desired
    
    if Animation_plot_gap == 0
        figure('WindowState','maximized');
        visualize(iteration,X,Y,u,v,T,p,e,rho,time(end))
    end


% Saving the results for easier plotting later on.
    if save_result == 1
    save(sprintf('M = %0.1f Ad = %d.mat',M,Adiabatic), 'u','v','T',...
        'p','e','rho','time','dtarray','sum_u','absdudt','X','Y')
    end

%% %%%% FUNCTIONS %%%%%%%%%

function [txx,tyy,txy,tyx,qx_dot,qy_dot]=InnerDerivatives...
    (xdir,ydir,u,v,T)

    % Calculates stresses and fluxes using the given xdir and ydir as the
    % desired direction for x derivative and y derivative. (Doesn't do
    % biasing. it just simply recieves the directions)

    global dx dy mu0 cp Pr T0 S1  Adiabatic

    % precalculates derivatives according to the input direction
    switch xdir
        case 'bwd'
            ux = ddx_bwd(u,dx);
            vx = ddx_bwd(v,dx);
            Tx = ddx_bwd(T,dx);
        case 'fwd'
            ux = ddx_fwd(u,dx);
            vx = ddx_fwd(v,dx);
            Tx = ddx_fwd(T,dx);
        case 'ctr'
            ux = ddx_central(u,dx);
            vx = ddx_central(v,dx);
            Tx = ddx_central(T,dx);
    end
    % 
   
    % 
    switch ydir
        case 'bwd'
            uy = ddy_bwd(u,dy);
            vy = ddy_bwd(v,dy);
            Ty = ddy_bwd(T,dy);
        case 'fwd'
            uy = ddy_fwd(u,dy);
            vy = ddy_fwd(v,dy);
            Ty = ddy_fwd(T,dy);
        case 'ctr'
            uy = ddy_central(u,dy);
            vy = ddy_central(v,dy);
            Ty = ddy_central(T,dy);
    end

     if Adiabatic == 1
        Ty(:,1) = 0;
     end


    % Get Temperature dependent paramteres

    mu_now = mu0*((T/T0).^(3/2)).*((T0+S1)./(T+S1));
    k_now = cp/Pr.*mu_now;

    % normal and shear stresses 
    txx = 2*mu_now.*(ux - 1/3*(ux + vy));
    tyy = 2*mu_now.*(vy - 1/3*(ux + vy));
    txy = mu_now.*(uy + vx);
    tyx = txy ;
    
    % heat fluxes 
    qx_dot = - k_now .* Tx;
    qy_dot = - k_now .* Ty;

    
end






function [E,F] = GenerateBiasedEF(direction, u , v , T , p , rho)
    
    % Computes E and F using the InnerDerivatives function. It
    % automatically finds the correct biasing to insert into
    % InnerDerivatives function based on "direction" input. "direction"
    % should be the direction of the outer derivative. for example it
    % should be forward for predictor level and backward for the corrector
    % level.
    % example: ydirE = direction for y derivative inside of E.



    % Biasing Logic Block
        ydirE = 'ctr';    
        xdirF = 'ctr';
        switch direction
          case 'fwd'                 % means predictor
              xdirE = 'bwd';
              ydirF = 'bwd';
          case 'bwd'                 % means corrector
              xdirE = 'fwd';   
              ydirF = 'fwd';
        end


    % preallocate E and F
    E = zeros([4,size(u)]);
    F = zeros([4,size(u)]);

    % Get inner derivatives for F
        [~,tyy,txy,~,~,qy_dot]=...
            InnerDerivatives(xdirF,ydirF,u,v,T);

        
    global cv 
    e = cv*T;
    Et = rho.*(e+1/2*(u.^2+v.^2));

    % Compute F
        F(1,:,:) = rho.*v; 
        F(2,:,:) = rho.*u.*v - txy; 
        F(3,:,:) = rho.*v.^2 + p - tyy; 
        F(4,:,:) = (Et+p).*v - u.*txy -v.*tyy + qy_dot;

    % Get inner derivatives for E
        [txx,~,txy,~,qx_dot,~]=...
            InnerDerivatives(xdirE,ydirE,u,v,T);

    % Compute E
        E(1,:,:) = rho.*u; 
        E(2,:,:) = rho.*u.^2 + p - txx;
        E(3,:,:) = rho.*u.*v - txy; 
        E(4,:,:) = (Et+p).*u - u.*txx -v.*txy + qx_dot;
       
end


function [u,v,T,p,rho]=imposeBC(u,v,T,p,rho)

    % Imposes boundary Conditions and updates rho

    global T0 p0 u_inf R Adiabatic
  
    % In-flow
    u(1,:) = u_inf;
    v(1,:) = 0;
    T(1,:) =  T0;
    p(1,:) = p0;
    % Farfield
    u(:,end) = u_inf;
    v(:,end) = 0;
    T(:,end) = T0;
    p(:,end) = p0;
    % Out-flow (right)
    u(end,:) = (4*u(end-1,:) - u(end-2,:))/3;
    v(end,:) = (4*v(end-1,:) - v(end-2,:))/3;
    p(end,:) = (4*p(end-1,:) - p(end-2,:))/3;
    T(end,:) = (4*T(end-1,:) - T(end-2,:))/3;
    % Wall
    u(:,1) = 0;
    v(:,1) = 0;
    p(:,1) = (4*p(:,2) - p(:,3))/3;
    T(:,1) = T0; 

    % enforce ADIABATIC 
            if Adiabatic == 1
                T(:,1) = (T(:,3)-4*(T(:,2))) / (-3);
              
            end
            
    % Leading Edge
    u(1,1) = 0;
    v(1,1) = 0;
    p(1,1) = p0;
    T(1,1) = T0; 


    % Update rho (efficiently)
    rho(1,:) = p(1,:) ./ T(1,:) ./ R;
    rho(:,1) = p(:,1) ./ T(:,1) ./ R;
    rho(:,end) = p(:,end) ./ T(:,end) ./ R;
    rho(end,:) = p(end,:) ./ T(end,:) ./ R;
        

end




function U_bar = predictor(U,E,F,dt) 
    U_bar=zeros(size(U));
    global dx dy
    for i=1:4
        Ei=squeeze(E(i,:,:));
        Fi=squeeze(F(i,:,:));
        Ui = squeeze(U(i,:,:));
        U_bar(i,:,:) = Ui - dt .* ( ddx_fwd(Ei,dx) + ddy_fwd(Fi,dy));
    end

end

function newU = corrector(U_bar,U,E_bar,F_bar,dt)
    global dx dy
    newU = zeros(size(U));
    for i=1:4
        Ei = squeeze(E_bar(i,:,:));
        Fi = squeeze(F_bar(i,:,:));
        Ui = squeeze(U(i,:,:));
        Ubari = squeeze(U_bar(i,:,:));

        newU(i,:,:) = 1/2 .* (Ui + Ubari ...
        - dt.*( ddx_bwd(Ei,dx) + ddy_bwd(Fi,dy)));
    end
end



 function dt = CFL(u,v,T,rho)
   global gamma Pr R T0 S1 mu0 dx dy
   mu = mu0*((T/T0).^(3/2)).*((T0+S1)./(T+S1));
   a = @(T) sqrt(gamma*R*T);

   v_prime = max(4/3.*mu.*(gamma*mu/Pr)./rho);
   dtCFL = [abs(u)./dx+abs(v)./dy+a(T).*sqrt(1/dx^2+1/dy^2)...
       +2*ones(size(u)).*v_prime.*(1/dx^2+1/dy^2)].^(-1);

   dt = min(0.5*dtCFL,[],"all");

 end



 %% Plotting Functions


function visualize(n,X,Y,u,v,T,p,e,rho,timenow)
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultTextInterpreter','latex');


   
   global iteration SchlierenDensity
    tcl = tiledlayout(2,3);
   
    titlepos = [0.5e-5,8.7e-6,0];
    nexttile
    pcolor(X,Y,u(:,:,n));shading interp;axis equal tight;
    xlabel('x');ylabel('y');title('x-velocity ($u$)','Position',titlepos);
    colorbarEden('u $\left[\frac{m}{s}\right]$');

    nexttile
    pcolor(X,Y,v(:,:,n));shading interp;axis equal tight;
    xlabel('x');ylabel('y');title('y-velocity ($v$)','Position',titlepos);
    colorbarEden('v $\left[\frac{m}{s}\right]$');

    nexttile
    pcolor(X,Y,p(:,:,n));shading interp;axis equal tight;
    xlabel('x');ylabel('y');title('Pressure','Position',titlepos);
    colorbarEden('p $\left[Pa\right]$');

    nexttile
    pcolor(X,Y,T(:,:,n));shading interp;axis equal tight;
    xlabel('x');ylabel('y');title('Temperature','Position',titlepos);
    colorbarEden('T $\left[K\right]$');
   
    a = nexttile();
    if SchlierenDensity == 1
        pcolor(X,Y,Schlieren(squeeze(rho(:,:,n))));
        colorbarEden('$\rho\ \left[\frac{kg}{m^3}\right]$',0,1);
        colormap(a,gray)
    else
    pcolor(X,Y,rho(:,:,n));shading interp;axis equal tight;
    colorbarEden('$\rho\ \left[\frac{kg}{m^3}\right]$');
    end

    shading interp;axis equal tight;
    xlabel('x');ylabel('y');
    title('Density\ ($\rho$) ','Position',titlepos);
    
    nexttile
    pcolor(X,Y,e(:,:,n));shading interp;axis equal tight;
    xlabel('x');ylabel('y');title('\ \ \ \  Internal\ Energy\ ($e$)','Position',titlepos);
    colorbarEden('$e\ \left[ \frac{J}{kg} \right]$');
    title(tcl,{sprintf(' t = %0.3g  ( %d / %d )', timenow , n , iteration ), '  '},...
        'FontSize',19,'Interpreter','latex');
    drawnow

end

function ConvergencePlot(sumu,dudtabs,time,n)
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultTextInterpreter','latex');

    % convergence block 
    converge=gcf;
    subplot(3,1,1)
    plot(3:n,sumu(3:end),'LineWidth',4);title('$\sum{u}_{ij}$');
    ylabel('$\frac{m}{s}$')

    subplot(3,1,2)
    plot(3:n,dudtabs(3:end),'LineWidth',4);
    title('$\sum{ \frac{|\Delta u|}{\Delta t}  }$');
    ylabel('$\frac{m}{s}$')
  
    subplot(3,1,3)
    plot(3:n,time(3:end),'LineWidth',4);title('$\Delta$t CFL');
    xlabel('iteration');ylabel('s')

    sgtitle('converging variables and $\Delta$t plots ')
    set(findall(converge,'-property','FontSize'),'FontSize',15)
    drawnow
     
end

function S = Schlieren(rho)
    nx = 75;
    ny = 80;
    L = 1e-5;
    H = 8e-6;
     dx = L / nx;
     dy = H / ny;
     beta = 0.8 ; kappa = 15;
     dummy = sqrt((ddx_fwd(rho,dx)).^2+(ddy_fwd(rho,dy)).^2);
     S = beta.* exp(-kappa.*dummy./max(dummy,[],"all")); 
end
