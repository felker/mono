%2D Time-Dependent RT solver
%Based on Jiang14 algorithm
%Monotonic Upwind Interpolation
%Mixed frame to O(v/c)

%------------------------ PARAMETERS ------------------ %
clear all;
close all;
nx = 100;
ny = 100;
ntheta = 4; %quadarture is only defined up to 12 in each direction, must be even
%order of the quadrature, N, refers to the number of mu-levels in the interval [-1, 1].

lx = 1.0;
ly = 1.0;
c = 1.0;
dx = lx/nx;
dy = ly/ny;
dt = 0.001;
nt = 100;

%Upwind monotonic interpolation scheme
method = 'donor'; 
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

%------------------------ PRE-TIMESTEPPING SETUP ------------------ %

%Calculate Radiation CFL numbers
cfl_mu = c*dt*abs(mu)*[1/dx 1/dy 0]';

%Calculate dot product of local fluid velocity and radiation rays
nv = zeros(nx,ny,na);
%cat(3,v,zeros(nx,ny,1)).*mu'; non-looping way?
for i=1:nx
    for j=1:ny
        for k=1:na
            nv(i,j,k) = v(i,j,1)*mu(na,1) + v(i,j,2)*mu(na,2);
        end
    end
end
        
%Explicit-Implicit operator splitting scheme
%-------------------------------------------
%assert(min(abs(cfl_mu) <= ones(ntheta,1))); 
intensity_new = intensity;
for i=1:nt
    intensity = intensity_new; 
    %Calculate the mean intensity for the mixed frame equations
    
    %Reapply boundary conditions
    intensity(1,:,1) = 1.0;
    %intensity(nz,ntheta/2+1:ntheta) = 0.0;
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
        intensity(:,:,j) = intensity(:,:,j) + ...
            alpha.*(mu(j,1)*dt*c/dx*(i_flux(:,:,1) - circshift(i_flux(:,:,1),[-1, 0, 0])) + ...
            mu(j,2)*dt*c/dy*(i_flux(:,:,2) - circshift(i_flux(:,:,2),[0 -1 0 ])));  
        
        %Advection with fluid at fluid velocity
        %Calculate current mean intensity
        mean_intensity = zeros(nx,ny);
        for k=1:nx
            for l=1:ny
                for m=1:na
                    mean_intensity(k,l) = mean_intensity(k,l) + intensity(k,l,m)*pw(m);
                end
            end
        end
        i_flux = upwind_interpolate2D(3*nv(:,:,j).*mean_intensity,method,mu(j,:),dt,dx,dy,c); 
        intensity(:,:,j) = intensity(:,:,j) + ...
            alpha.*(mu(j,1)*dt*c/dx*(i_flux(:,:,1) - circshift(i_flux(:,:,1),[-1, 0, 0])) + ...
            mu(j,2)*dt*c/dy*(i_flux(:,:,2) - circshift(i_flux(:,:,2),[0 -1 0 ])));   
    return;
        %Substep #2: Implicit absorption source terms

        %Substep #3: Implicit scattering source terms
    end
    time = dt*i
    if ~mod(i,10)
            %Plot each angular intensity (recall, ntheta must be even)
            for j=1:ntheta
            hi = subplot(2,1,j); 
            plot(hi,zz,intensity(:,j,1),'-o');
            x_label = sprintf('mu = %f',mu(j));
            time_title = sprintf('t = %f (s)',time);
            legend(method);
            xlabel(x_label);
            title(time_title); 
            end
            pause(1.0);
    end
end

