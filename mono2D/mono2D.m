%2D Time-Dependent RT solver
%Based on Jiang14 algorithm
%Monotonic Upwind Interpolation
%Mixed frame to O(v/c)

%RENAME RHO_A RHO_S TO SIGMAS

%------------------------ PARAMETERS ------------------ %
clear all;
close all;
nx = 32;
ny = 32;
ntheta = 4; %quadarture is only defined up to 12 in each direction, must be even
%order of the quadrature, N, refers to the number of mu-levels in the interval [-1, 1].

lx = 1.0;
ly = 1.0;
c = 1.0;
dx = lx/nx;
dy = ly/ny;
dt = 0.002;
nt = 30;
out_count = 0;
%Upwind monotonic interpolation scheme
method = 'van Leer'; 
time_centering = 0; %explicit = 0, implicit = 1, CN=1/2
delta_c = 1e-10; %error tolerance for NR convergence

%------------------------ BUILD MESH ------------------ %
%Fixed Background Fluid Velocity 
v = zeros(nx,ny,2); 
%Fixed fluid temperature
temp = ones(nx,ny);  %do we need to use the stefan boltzmann constant?
%Gas density
density = zeros(nx,ny);

%Radiation Angular Discretization, S_n discrete ordinates. 
[na, mu, pw, lw] = angular_quad2D(ntheta);

%Radiation Spatial Discretization
%Radiation points are centered on the fluid cells
xx=linspace(0,lx,nx)';
yy=linspace(0,ly,ny)';

%Monochromatic specific intensity, boundary conditions at 2,nz-1 
intensity = zeros(nx,ny,na); 
mean_intensity = zeros(nx,ny); 

%--------------- JIANG14 variables, dimensionless constants -------------- %
a_r =  5.670373e-8; %stefan boltzman constant
adiabatic = 5/3; %gamma=adiabatic index
R = 1; %8.311; %ideal gas constant

%characteristic values (arbitrary?)
a_0 = 0.1; %characteristic velocity
T_0 = 1.0; %characteristic temperature
P_0 = a_r;%1.0; %characteristic gas pressure

C = c/a_0; %dimensionlesss speed of light
P = a_r*T_0^4/P_0;  %measure of relative importance of rad and gas pressure

v(:,:,1) = 0.3*C; 
%------------------------ PROBLEM SETUP ------------------ %
%For 2D and verification of Jiang14, we no longer use the source function
%or extinction probability 
%Absorption opacities
rho_a = ones(nx,ny);
%Scattering opacities
rho_s = zeros(nx,ny);
%Fluid density
density = ones(nx,ny);

%------------------------ PRE-TIMESTEPPING SETUP ------------------ %
%Calculate Radiation CFL numbers
cfl_mu = C*dt*abs(mu)*[1/dx 1/dy 0]';
%set dt to have max cfl of 0.4
dt = dt*0.4/max(cfl_mu);
%Recalculate radiation CFL
cfl_mu = C*dt*abs(mu)*[1/dx 1/dy 0]';

%------------------------ VELOCITY  ------------------------------ %
%Calculate dot product of local fluid velocity and radiation rays
nv = zeros(nx,ny,na);
%Calculate loval velocity magnitude
vCsquare = zeros(nx,ny);
vsquare = zeros(nx,ny);
absV = zeros(nx,ny);
for i=1:nx
    for j=1:ny
        vsquare(i,j) = (v(i,j,1)^2 + v(i,j,2)^2);         
        vCsquare(i,j) = (v(i,j,1)^2 + v(i,j,2)^2)/C^2; 
        absV(i,j) = sqrt(v(i,j,1)^2 + v(i,j,2)^2);
        for k=1:na
            nv(i,j,k) = v(i,j,1)*mu(k,1) + v(i,j,2)*mu(k,2);
        end
    end
end


%Explicit-Implicit operator splitting scheme
%-------------------------------------------
%assert(min(abs(cfl_mu) <= ones(ntheta,1))); 
intensity_new = intensity;

%initial condition for equilibirium test
%intensity(:,:,:) = 1/(4*pi);
% %this is a projection!
% hold on;
% for i=1:na
%     if mu(i,3) > 0.5
%         %mu_z = 0.88
%         quiver(0,0,intensity(2,2,i).*mu(i,1), intensity(2,2,i).*mu(i,2),'-b');
%     else
%         quiver(0,0,intensity(2,2,i).*mu(i,1), intensity(2,2,i).*mu(i,2),'-k');
%     end
% 
% end
% 
% return;
for i=1:nt
    %Reapply boundary conditions
    intensity(nx,:,:) = 0.0;
    intensity(:,ny,:) = 0.0;
    intensity(:,1,:) = 0.0;
    intensity(1,:,:) = 0.0;

    %intensity(1,2:ny-1,:) = 0.0; %all first row elements (physically, left bndry)

    %Setup for 2D crossing beams in Jiang14 and Davis12
    %Periodic boundary conditions on y-boundaries, only for appropriate
    %rays
    %intensity(:,1,1:na/2) = intensity(:,ny-1,1:na/2); 
    %intensity(:,ny,na/2:na) = intensity(:,2,na/2:na); 
    %Inject beam at center of LHS boundary with mu_x=0.333, -0.3333
    intensity(1,ny/2,1) = 1.0;
    %intensity(1,ny/2,7) = 1.0;
    
    %Box periodic boundary conditions
    %intensity(:,1,1:na/2) = intensity(:,ny-1,1:na/2); 
    %intensity(:,ny,na/2:na) = intensity(:,2,na/2:na); 
    %mu_x arent easily divided into positive and negative...
    %intensity(1,:,1:3) = intensity(nx-1,:,1:3); 
    %intensity(1,:,7:9) = intensity(nx-1,:,7:9); 
    %intensity(nx,:,4:6) = intensity(2,:,4:6); 
    %intensity(nx,:,10:12) = intensity(2,:,10:12); 
    
    %Calculate current mean intensity
    mean_intensity = zeros(nx,ny);
    for k=1:nx
        for l=1:ny
            for m=1:na
                mean_intensity(k,l) = mean_intensity(k,l) + intensity(k,l,m)*pw(m);
            end
        end
    end
    
    for j=1:na
        %Substep #1: Explicitly advance transport term
        %Split the transport term into diffusion and advection terms
        %Photon diffusion at nearly the speed of light
        tau = (10*(dx/mu(j,1) + dy/mu(j,2))*(rho_a+rho_s)).^2; %Eq 14 in Jiang14
        %Correct for division by zero 
        for k=1:nx
            for l=1:ny
                if tau(k,l) > 0 
                    alpha(k,l) = sqrt((1-exp(-tau(k,l)))./(tau(k,l)));
                else
                    alpha(k,l) = 1.0;
                end
            end
        end
        %coefficients vs u_i for van leer interpolation?
%          i_flux = upwind_interpolate2D(intensity(:,:,j),method,mu(j,:),dt,dx,dy,alpha*C);
         i_flux = upwind_interpolate2D(intensity(:,:,j) - 3*nv(:,:,j).*mean_intensity,method,mu(j,:),dt,dx,dy,alpha*C);

         A = circshift(i_flux(:,:,1),[-1 0]);
         B = circshift(i_flux(:,:,2),[0 -1]);
         intensity(2:nx-1,2:ny-1,j) = intensity(2:nx-1,2:ny-1,j) + ...
                         (mu(j,1)*dt*C/dx*(i_flux(2:nx-1,2:ny-1,1) - A(2:nx-1,2:ny-1)) + ...
             mu(j,2)*dt*C/dy*(i_flux(2:nx-1,2:ny-1,2) - B(2:nx-1,2:ny-1))); 
             %alpha(2:nx-1,2:ny-1).*(mu(j,1)*dt*C/dx*(i_flux(2:nx-1,2:ny-1,1) - A(2:nx-1,2:ny-1)) + ...
  
 
         %symmetry test
         %fliplr(intensity(2:nx-1,2:ny-1,j)) - intensity(2:nx-1,2:ny-1,j) 
 
         %Advection with fluid at fluid velocity
          i_flux = upwind_interpolate2D(3*nv(:,:,j).*mean_intensity,method,mu(j,:),dt,dx,dy,absV); 
          intensity(:,:,j) = intensity(:,:,j) + ...
              (mu(j,1)*dt/dx*(i_flux(:,:,1) - circshift(i_flux(:,:,1),[-1, 0, 0])) + ...
              mu(j,2)*dt/dy*(i_flux(:,:,2) - circshift(i_flux(:,:,2),[0 -1 0 ])));  
          end
        %Substep #2: Implicitly advance the absorption source terms at each
        %cell, all angles in a cell
%         for k=2:nx-1
%             for l=2:ny-1
%             %Use Newton-Raphson. System is nonlinear in T^n+1. Use previous
%             %values as initial guesses. Need to do implicitly because
%             %thermalization is fast
%             I_old = squeeze(intensity(k,l,:));
%             temp_old = temp(k,l);
% 
%             %Tensor-like quadrature term in the equations
%             vvnn = zeros(na,1);
%             for r=1:na 
%                 vvnn(r) = v(k,l,1)^2*mu(r,1)^2 + v(k,l,2)^2*mu(r,2)^2 + ...
%                     2*v(k,l,1)*mu(r,1)*v(k,l,2)*mu(r,2);
%             end;     
%             
%             %REDUCE TO 1D NR solve
%             
%             A = zeros(na,1);
%             B = zeros(na,1);
%             D = zeros(na,1);
%             A(na) = (1 +dt*rho_a(k,l)*C - dt*rho_a(k,l)*nv(k,l,na));
%             for r=1:na-1
%                 A(r) = A(na)/(1 +dt*rho_a(k,l)*C - dt*rho_a(k,l)*nv(k,l,r)); 
%                 B(r) = 3.0*dt*rho_a(k,l)*(nv(k,l,r) - nv(k,l,na))/(4*pi*(1+dt*rho_a(k,l)*C*(1-nv(k,l,r)/C)));
%                 D(r) = (intensity(k,l,r) - intensity(k,l,na))/(1+dt*rho_a(k,l)*C*(1-nv(k,l,r)/C));
%             end
%             A(na) = 1.0;
%             B(na) = 0.0;
%             D(na) = 0.0;
%   
%             %My equilibrium function, and derivative
% %             E = zeros(na,1);
% %             for r=1:na 
% %                 E(r) = (I_old(r) + C*rho_a(k,l)*temp_old^4/(4*pi))/(1 +dt*rho_a(k,l)*C);
% %             end;
% %             
% %             fn_handle = @(t) deal(t - temp_old + dt*rho_a(k,l)*P*C/(density(k,l)*R)*(t^4 - ...
% %                 4*pi*E'*pw), 1 + dt*rho_a(k,l)*P*C/(density(k,l)*R)*(1- 1/(1+dt*rho_a(k,l)*C))*4*t^3);      
%             
%             %Match coefficients in hydro_to_rad.c
%             Jfactor = P*(2+dt*rho_a(k,l)*C*(1+vCsquare(k,l)));
%             J1 = - density(k,l)*R/(4*pi*(adiabatic-1)*Jfactor); 
%             J2 = dt*rho_a(k,l)*P*C*(1+vCsquare(k,l))/(4*pi*Jfactor);
%             J3 = 2*P*mean_intensity(k,l)/Jfactor + density(k,l)*R*temp_old/(4*pi*(adiabatic-1)*Jfactor); 
%             
%             coef1 = 1 + dt*rho_a(k,l)*C*(1-nv(k,l,na)/C) + dt*rho_a(k,l)*((vsquare(k,l) + vvnn(r)).*A)'*pw/C;
%             coef2 = dt*rho_a(k,l)*C*(1+3*nv(k,l,na)/C)/(4*pi) - dt*rho_a(k,l)*((vsquare(k,l) + vvnn(r)).*B)'*pw/C;
%             coef3 = intensity(k,l,na) - dt*rho_a(k,l)*((vsquare(k,l) + vvnn(r)).*D)'*pw/C;
%             
%             Tcoef = coef1*J1;
%             T4coef = coef1*J2;
%             Tconst = coef1*J3;
% 
%             T4coef = T4coef - (coef1*B + coef2*A)'*pw;
%             Tconst = Tconst - (coef1*D + coef3*A)'*pw;
% 
%             Tmin = max(-0.25*Tcoef/T4coef,0.0).^(1/3); 
%             Tmax = (0.5*density(k,l)*vsquare(k,l) +4*pi*mean_intensity(k,l))*(adiabatic -1.0)/(density(k,l)*R) + temp_old; 
%             %For now, assume that my function (which has no need to pass in
%             %the coeff) is correct. use modified rtsafe
%             %temp(k,l) = rtsafe(fn_handle,Tmin,Tmax,delta_c);
%             %ATHENA equilibrium function, and derivative
%             fn_handle = @(t,coef1,coef2,coef3,coef4) deal(coef1*t^4 + coef2*t + coef3,4.0*coef1*t^3 + coef2); 
%             temp(k,l) = rtsafe_orig(fn_handle,Tmin,Tmax,delta_c, T4coef, Tcoef, Tconst, 0.0);
%                    
%             %Explicitly update all other rays 
%             intensity(k,l,na) = (coef2*temp(k,l)^4 + coef3)/coef1;
%             intensity(k,l,1:na-1) = A(1:na-1)*intensity(k,l,na) + B(1:na-1) + D(1:na-1); 
%             
%             %Fake update gas quantities as it absorbs internal energy
%             %density?
%             %change in gas internal energy
%             dEt = (density(k,l)*R/adiabatic)*(temp(k,l) - temp_old);
%             dGasMomentum = zeros(2,1);
%             for r=1:2
%                 for n=1:na
%                     dGasMomentum(r) = dGasMomentum(r) - P/C*(intensity(k,l,n) - I_old(n))*mu(n,r)*pw(n);                 
%                 end
%             end
%             %no change in momentum due to isotropy
%             dGasKE =((density(k,l)*squeeze(v(k,l,:)) + dGasMomentum)'*(density(k,l)*squeeze(v(k,l,:)) + dGasMomentum) - ...
%                 (density(k,l)*squeeze(v(k,l,:)))'*(density(k,l)*squeeze(v(k,l,:))))/(2*density(k,l));
%             
%             end
%         end
                    mean_intensity = zeros(nx,ny);
    for k=1:nx
        for l=1:ny
            for m=1:na
                mean_intensity(k,l) = mean_intensity(k,l) + intensity(k,l,m)*pw(m);
            end
        end
    end
        time = dt*i; %#ok<NOPTS>
        
        %Substep #3: Implicitly advance the scattering source terms
    
%------------------------ OUTPUT ------------------ %
    if ~mod(i,5)
%      out_count = out_count +1; 
%     time_plot(out_count) = time; 
%     y_out(out_count) = 4*pi*mean_intensity(nx/2,ny/2);
%     y_out2(out_count) = temp(nx/2,ny/2)^4;
% % %     figure(1);
%       plot(time_plot,y_out,'-o',time_plot,y_out2,'-o');
%       legend('E_r','T^4');

%      xlabel('time (s)');
%      ylabel('E_r,T^4');
%      title('\sigma_a = 1');
   
%     figure(2);
 %      pcolor(xx(2:nx-1,1),yy(2:ny-1,1),mean_intensity(2:nx-1,2:ny-1)')
   %colorbar
%         hi = subplot(2,3,1); 
%         h = pcolor(xx,yy,intensity(:,:,1)');
%         set(h, 'EdgeColor', 'none');
%         x_label = sprintf('mu =(%f, %f)',mu(1,1),mu(1,2));
%         xlabel(x_label);
%         colorbar
             %hi = subplot(1,2,1); 

 %            h = pcolor(yy(2:ny-1),xx(2:nx-1),mean_intensity(2:nx-1,2:ny-1));
             %set(h, 'EdgeColor', 'none');
             %x_label = sprintf('mu =(%f, %f)',mu(1,1),mu(1,2));
%             time_title = sprintf('t = %f (s)',time);
             %xlabel(x_label);
 %            title(time_title);
  %           colorbar
%                           hi = subplot(1,2,2); 
%              h = pcolor(yy,xx,intensity(:,:,7));
%              set(h, 'EdgeColor', 'none');
%              x_label = sprintf('mu =(%f, %f)',mu(7,1),mu(7,2));
%              time_title = sprintf('t = %f (s)',time);
%              xlabel(x_label);
%              title(time_title);
%              colorbar
            %Plot each angular intensity (recall, ntheta must be even)
              for j=1:na/2
                  hi = subplot(2,3,j); 
                  h = pcolor(xx,yy,intensity(:,:,j)');
                  set(h, 'EdgeColor', 'none');
                  x_label = sprintf('mu =(%f, %f)',mu(j,1),mu(j,2));
                  time_title = sprintf('t = %f (s)',time);
                  xlabel(x_label);
                  title(time_title);
                  colorbar
              end

%Section 5.2 tests
%this is a projection!
% for i=1:na
%     if mu(i,3) > 0.5
%         %mu_z = 0.88
%         quiver(0,0,intensity(2,2,i).*mu(i,1), intensity(2,2,i).*mu(i,2),'-b');
%     else
%         quiver(0,0,intensity(2,2,i).*mu(i,1), intensity(2,2,i).*mu(i,2),'-k');
%     end
% hold on;
% 
% end
% hold off; 
            pause(2.0);
    end
end

