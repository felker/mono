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
nt = 6000;

%Upwind monotonic interpolation scheme
method = 'van Leer'; 
time_centering = 0; %explicit = 0, implicit = 1, CN=1/2
MAX_ITER = 1000; % maximum number of iterations for Newton-Raphson
delta_c = 1e-4; %error tolerance for NR convergence

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
%------------------------ PROBLEM SETUP ------------------ %

%For 2D and verification of Jiang14, we no longer use the source function
%or extinction probability 
%Absorption opacities
rho_a = ones(nx,ny);
%Scattering opacities
rho_s = zeros(nx,ny);
%Fluid density
density = ones(nx,ny);

%Open questions: how does one convert from the thermal source function to
%rho_a, rho_s for LTE
% How can we ensure that the cfl number per angle is the same without
% changing the timestep?
% Why is stability CFL=0.5?
v(:,:,2) = 0.0; 

%------------------------ PRE-TIMESTEPPING SETUP ------------------ %

%Calculate Radiation CFL numbers
cfl_mu = c*dt*abs(mu)*[1/dx 1/dy 0]';
%set dt to have max cfl of 0.4
dt = dt*0.4/max(cfl_mu);
cfl_mu = c*dt*abs(mu)*[1/dx 1/dy 0]';

%--------------- JIANG14 variables, dimensionless constants -------------- %
a_r =  5.670373e-8; %stefan boltzman constant
adiabatic = 5/3; %gamma=adiabatic index
R = 8.311; %ideal gas constant

%characteristic values (arbitrary?)
a_0 = 0.1; %characteristic velocity
T_0 = 1.0; %characteristic temperature
P_0 = 1/a_r;%1.0; %characteristic gas pressure

C = c/a_0; %dimensionlesss speed of light
P = a_r*T_0^4/P_0;  %measure of relative importance of rad and gas pressure

%Calculate dot product of local fluid velocity and radiation rays
nv = zeros(nx,ny,na);
%Calculate loval velocity magnitude
vC = zeros(nx,ny);
for i=1:nx
    for j=1:ny
        vC(i,j) = (v(i,j,1)^2 + v(i,j,2)^2)/C^2; 
        for k=1:na
            nv(i,j,k) = v(i,j,1)*mu(k,1) + v(i,j,2)*mu(k,2);
        end
    end
end

%Explicit-Implicit operator splitting scheme
%-------------------------------------------
%assert(min(abs(cfl_mu) <= ones(ntheta,1))); 
intensity_new = intensity;
for i=1:nt
    %Reapply boundary conditions
    %intensity(nx,:,:) = 0.0;
    %intensity(:,ny,:) = 0.0;
    %intensity(:,1,:) = 0.0;
    %intensity(1,2:ny-1,:) = 0.0; %all first row elements (physically, left bndry)

    %Setup for 2D crossing beams in Jiang14 and Davis12
    %Periodic boundary conditions on y-boundaries, only for appropriate
    %rays
    %intensity(:,1,1:na/2) = intensity(:,ny-1,1:na/2); 
    %intensity(:,ny,na/2:na) = intensity(:,2,na/2:na); 
    %Inject beam at center of LHS boundary with mu_x=0.333, -0.3333
    %intensity(1,ny/2,1) = 1.0;
    %intensity(1,ny/2,7) = 1.0;
    
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
    
        i_flux = upwind_interpolate2D(intensity(:,:,j),method,mu(j,:),dt,dx,dy,c);
        A = circshift(i_flux(:,:,1),[-1 0]);
        B = circshift(i_flux(:,:,2),[0 -1]);
        intensity(2:nx-1,2:ny-1,j) = intensity(2:nx-1,2:ny-1,j) + ...
            alpha(2:nx-1,2:ny-1).*(mu(j,1)*dt*C/dx*(i_flux(2:nx-1,2:ny-1,1) - A(2:nx-1,2:ny-1)) + ...
            mu(j,2)*dt*c/dy*(i_flux(2:nx-1,2:ny-1,2) - B(2:nx-1,2:ny-1)));  

        %fliplr(intensity(2:nx-1,2:ny-1,j)) - intensity(2:nx-1,2:ny-1,j) 

        %Advection with fluid at fluid velocity
        i_flux = upwind_interpolate2D(3*nv(:,:,j).*mean_intensity,method,mu(j,:),dt,dx,dy,c); 
        intensity(:,:,j) = intensity(:,:,j) + ...
            sqrt(v(:,:,1).^2 + v(:,:,2).^2).*(mu(j,1)*dt/dx*(i_flux(:,:,1) - circshift(i_flux(:,:,1),[-1, 0, 0])) + ...
            mu(j,2)*dt/dy*(i_flux(:,:,2) - circshift(i_flux(:,:,2),[0 -1 0 ])));  
        end
        %Substep #2: Implicitly advance the absorption source terms at each
        %cell
        for k=2:nx-1
            for l=2:ny-1
        %Use Newton-Raphson. System is nonlinear in T^n+1. Use previous
        %values as initial guesses. Need to do implicitly because
        %thermalization is fast
            temp_step = temp(k,l);
            IN_step = intensity(k,l,na);
            I_old = squeeze(intensity(k,l,:));
            temp_old = temp_step;
            for m =1:MAX_ITER
                %Build RHS vector representing the original nonlinear
                %function, restricted to only the final two components
                F = zeros(2,1); 
                J = zeros(2,2);
                %Tensor like quadrature term in the equations
                vvnn = zeros(na,1);
                A = zeros(na,1);
                for r=1:na 
                    vvnn(r) = v(k,l,1)^2*mu(r,1)^2 + v(k,l,2)^2*mu(r,2)^2 + ...
                        2*v(k,l,1)*mu(r,1)*v(k,l,2)*mu(r,2);
                    %Corresponds to Am[] in ATHENA
                    A(r) = (1 +dt*rho_a(k,l)*C - dt*rho_a(k,l)*nv(k,l,r)); 
                end;
                J(1,1) = A(na) + rho_a(k,l)*dt*vC(k,l)*(((A(na)*ones(na-1,1))./A(1:na-1))'*pw(1:na-1) + ...
                    pw(na)) + rho_a(k,l)*dt/C*(((vvnn(1:na-1).*pw(1:na-1))./A(1:na-1))'*A(na)*ones(na-1,1) + ...
                    pw(na)*vvnn(na));
                J(1,2) = -4*dt*rho_a(k,l)*C*(1+nv(k,l,na)/C)*temp_step^3/pi; 
                J(2,1) = -dt*(adiabatic -1)/(R*density(k,l))*(-8*pi*P*rho_a(k,l)*(((A(na)*squeeze(nv(k,l,1:na-1))./A(1:na-1))'*pw(1:na-1)...
                    + nv(k,l,na)*pw(na)) - c*vC(k,l)*(A(na)*ones(1,na-1)*((pw(1:na-1)./A(1:na-1))) + vvnn(na)*pw(na))) + ...
                    P*C*(1- vC(k,l))*rho_a(k,l)*4*pi*(pw(1:na-1)'*(A(na)*ones(na-1,1)./A(1:na-1)) + pw(na)));
                J(2,2) = 1 +4*(adiabatic -1)/(R*density(k,l))*dt*P*C*(1-vC(k,l))*rho_a(k,l)*temp_step^3; 
                
                %this item is from subtracting na-th line from lth line
                D = ((A(na)*IN_step -I_old(na))*ones(na-1,1) + I_old(1:na-1))./A(1:na-1); 
                
                F(1) = IN_step - I_old(na) - dt*(C*rho_a(k,l)*(temp_step^4/(4*pi) - IN_step) + ...
                nv(k,l,na)*rho_a(k,l)*(3*temp_step^4/(4*pi) + IN_step) - rho_a(k,l)*vC(k,l)*(pw(1:na-1)'*D + ...
                pw(na)*IN_step) - rho_a(k,l)/C*((pw(1:na-1).*vvnn(1:na-1))'*D + pw(na)*vvnn(na)*IN_step));
            
                F(2) = temp_step - temp_old - (adiabatic -1)*dt/(R*density(k,l))*(-8*pi*P*rho_a(k,l)*( ...
                    (squeeze(nv(k,l,1:na-1)).*pw(1:na-1))'*D + nv(k,l,na)*pw(na)*IN_step - C*vC(k,l)*(D'*pw(1:na-1) + pw(na)*IN_step) - ...
                    1/C*((vvnn(1:na-1).*pw(1:na-1))'*D + vvnn(na)*pw(na)*IN_step)) - P*C*(1-vC(k,l))*rho_a(k,l)*( ...
                    temp_step^4 - 4*pi*(pw(1:na-1)'*D + pw(na)*IN_step)));

                %Solve 2x2 system for I^n+1_na and T^n+1 steps
                delta_step = J\-F; 
                IN_step = IN_step + delta_step(1);
                temp_step = temp_step + delta_step(2);
                %check convergence criteria. Maybe change this to delta
                %intensity?
                if max(abs(delta_step./[IN_step temp_step]')) || (delta_step(1) == 0 && delta_step(2) == 0 && IN_step == 0 && temp_step ==0) < delta_c
                    break;
                end  
            end %done with NR solve at cell
               %Explicitly update all other rays 
                %for n=1:na-1
                    intensity(k,l,1:na-1) = ((A(na)*IN_step -I_old(na))*ones(na-1,1) + I_old(1:na-1))./A(1:na-1); 

                    %intensity(k,l,n) = (1-dt*(1+nv(k,l,n)/c)*rho_a(k,l))^-1*((IN_step - I_old(na))/dt + ...
                    %    rho_a(k,l)*IN_step + 3*temp_step^4/(4*pi*c)*rho_a(k,l)*(nv(k,l,n) - nv(k,l,na)) ...
                    %    - nv(k,l,na)/c*rho_a(k,l)*IN_step);           
                %end
                intensity(k,l,na) = IN_step; 
                temp(k,l) = temp_step;
                                return;
            end
        end
        
        
        %Substep #3: Implicitly advance the scattering source terms
    
    time = dt*i; %#ok<NOPTS>
%------------------------ OUTPUT ------------------ %
    if ~mod(i,100)
     %  pcolor(xx(2:nx-1,1),yy(2:ny-1,1),mean_intensity(2:nx-1,2:ny-1)')
   %colorbar
%         hi = subplot(2,3,1); 
%         h = pcolor(xx,yy,intensity(:,:,1)');
%         set(h, 'EdgeColor', 'none');
%         x_label = sprintf('mu =(%f, %f)',mu(1,1),mu(1,2));
%         xlabel(x_label);
%         colorbar
             %hi = subplot(1,2,1); 
             h = pcolor(yy,xx,mean_intensity);
             set(h, 'EdgeColor', 'none');
             %x_label = sprintf('mu =(%f, %f)',mu(1,1),mu(1,2));
             time_title = sprintf('t = %f (s)',time);
             %xlabel(x_label);
             title(time_title);
             colorbar
%                           hi = subplot(1,2,2); 
%              h = pcolor(yy,xx,intensity(:,:,7));
%              set(h, 'EdgeColor', 'none');
%              x_label = sprintf('mu =(%f, %f)',mu(7,1),mu(7,2));
%              time_title = sprintf('t = %f (s)',time);
%              xlabel(x_label);
%              title(time_title);
%              colorbar
            %Plot each angular intensity (recall, ntheta must be even)
%              for j=1:na
%              hi = subplot(2,6,j); 
%              h = pcolor(xx,yy,intensity(:,:,j)');
%              set(h, 'EdgeColor', 'none');
%              x_label = sprintf('mu =(%f, %f)',mu(j,1),mu(j,2));
%              time_title = sprintf('t = %f (s)',time);
%              %legend(method);
%              xlabel(x_label);
%              title(time_title);
%              colorbar
%              end
            pause(1.0);
    end
end

