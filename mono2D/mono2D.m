%2D Time-Dependent RT solver
%Based on Jiang14 algorithm
%Monotonic Upwind Interpolation
%Mixed frame to O(v/c)

%------------------------ PARAMETERS ------------------ %
clear all;
close all;
nx = 40;
ny = 40;
ntheta = 4; %quadarture is only defined up to 12 in each direction, must be even
%order of the quadrature, N, refers to the number of mu-levels in the interval [-1, 1].

lx = 1.0;
ly = 1.0;
c = 1.0;
dx = lx/nx;
dy = ly/ny;
dt = 0.002;
nt = 1000;

%Upwind monotonic interpolation scheme
method = 'van Leer'; 
time_centering = 0; %explicit = 0, implicit = 1, CN=1/2

%------------------------ BUILD MESH ------------------ %
%Fixed Background Fluid Velocity 
v = zeros(nx,ny,2); 

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
rho_a = zeros(nx,ny);
%Scattering opacities
rho_s = zeros(nx,ny);

%Open questions: how does one convert from the thermal source function to
%rho_a, rho_s for LTE
% How can we ensure that the cfl number per angle is the same without
% changing the timestep?
% Why is stability CFL=0.5?
v(:,:,1) = 0.0; 

%------------------------ PRE-TIMESTEPPING SETUP ------------------ %

%Calculate Radiation CFL numbers
cfl_mu = c*dt*abs(mu)*[1/dx 1/dy 0]';
%set dt to have max cfl of 0.4
dt = dt*0.4/max(cfl_mu);
cfl_mu = c*dt*abs(mu)*[1/dx 1/dy 0]';

%Calculate dot product of local fluid velocity and radiation rays
nv = zeros(nx,ny,na);
for i=1:nx
    for j=1:ny
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
    intensity(1,2:ny-1,:) = 1.0; %all first row elements (physically, left bndry)

    %Setup for 2D crossing beams in Jiang14 and Davis12
    %Periodic boundary conditions on x-boundaries, flip directions 
    %intensity(nx,:,:) = intensity(2,:,:); 
    %intensity(1,:,:) = intensity(nx-1,:,:); 
    %Inject beam at center of bottom boundary with mu_x=0.333, -0.3333
    %intensity(:,1,1) = 1.0;
    %intensity(:,1,4) = 1.0;
    
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
            alpha(2:nx-1,2:ny-1).*(mu(j,1)*dt*c/dx*(i_flux(2:nx-1,2:ny-1,1) - A(2:nx-1,2:ny-1)) + ...
            mu(j,2)*dt*c/dy*(i_flux(2:nx-1,2:ny-1,2) - B(2:nx-1,2:ny-1)));  

        %fliplr(intensity(2:nx-1,2:ny-1,j)) - intensity(2:nx-1,2:ny-1,j) 

        %Advection with fluid at fluid velocity
        i_flux = upwind_interpolate2D(3*nv(:,:,j).*mean_intensity,method,mu(j,:),dt,dx,dy,c); 
        %might need to remove the mu(j) here
        intensity(:,:,j) = intensity(:,:,j) + ...
            sqrt(v(:,:,1).^2 + v(:,:,2).^2).*(mu(j,1)*dt/dx*(i_flux(:,:,1) - circshift(i_flux(:,:,1),[-1, 0, 0])) + ...
            mu(j,2)*dt/dy*(i_flux(:,:,2) - circshift(i_flux(:,:,2),[0 -1 0 ])));   
        %Substep #2: Implicitly advance the absorption source terms

        %Substep #3: Implicitly advance the scattering source terms
 
    end
    time = dt*i %#ok<NOPTS>
    if ~mod(i,1000)
     %  pcolor(xx(2:nx-1,1),yy(2:ny-1,1),mean_intensity(2:nx-1,2:ny-1)')
   %colorbar
%         hi = subplot(2,3,1); 
%         h = pcolor(xx,yy,intensity(:,:,1)');
%         set(h, 'EdgeColor', 'none');
%         x_label = sprintf('mu =(%f, %f)',mu(1,1),mu(1,2));
%         xlabel(x_label);
%         colorbar

            %Plot each angular intensity (recall, ntheta must be even)
             for j=1:na
             hi = subplot(2,6,j); 
             h = pcolor(xx,yy,intensity(:,:,j)');
             set(h, 'EdgeColor', 'none');
             x_label = sprintf('mu =(%f, %f)',mu(j,1),mu(j,2));
             time_title = sprintf('t = %f (s)',time);
             %legend(method);
             xlabel(x_label);
             title(time_title);
             colorbar
             end
            pause(1.0);
    end
end

