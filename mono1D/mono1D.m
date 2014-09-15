%1D Time-Dependent RT solver
%Based on StoneMihalas92 algorithm
%Monotonic Upwind Interpolation

%Parameters
clear all;
close all;
nz = 100;
dt = 0.008; %cfl=1
%dt = 0.15; %cfl =15
%dt = 0.005; %cfl = 0.5

%nt = 160;
nt = 100; 
ntheta = 2; %quadarture is only defined up to 12 in each direction, must be even
%order of the quadrature, N, refers to the number of mu-levels in the interval [-1, 1].

lz = 1.0;
c = 1.0;
dz = lz/nz;

%Angular Discretization, Sn discrete ordinates. 
[mu, w] = angular_quad1D(ntheta);

%Spatial Discretization
%Radiation gridpoints are centered on the fluid cells
zz=linspace(0,lz,nz)';

%Monochromatic specific intensity, boundary conditions at 2,nz-1 
intensity = zeros(nz,ntheta); 

%Total opacity 
%X_tot = 10*ones(nz,1); 
%X_tot = 2*ones(nz,1); 
X_tot = zeros(nz,1);
%Photon destruction probability eps = X_abs/X_tot, LTE
destruction_probability = zeros(nz,1); 

%Thermal blackbody source function
T = 1000; % in K
temperatures = T*ones(nz,1);
B = 0; %1e-3;
thermal_source = B*temperatures;

%Calculate source function
source_function = destruction_probability.*thermal_source; 

cfl_mu = c*dt/dz.*mu 
intensity_new = intensity; 
%Upwind monotonic interpolation scheme
method = 'van Leer'; 
time_centering = 1; %explicit = 0, implicit = 1, CN=1/2

%Explicit-Implicit operator splitting scheme
%-------------------------------------------
assert(min(abs(cfl_mu) <= ones(ntheta,1))); 
intensity_substep = intensity;
for i=1:nt
    intensity = intensity_new; 
    %Reapply boundary conditions
    intensity(1,1:ntheta/2) = 1.0;
    intensity(nz,ntheta/2+1:ntheta) = 0.0;

    %Substep #1: Explicit advection terms
    %calculate the fluxes at the cell interfaces
    i_flux = upwind_interpolate(intensity,method,mu,dt,dz,c); 
    %i-th index corresponds to left bndry of ith cell. i+1th is right bndry
    intensity_substep = intensity + (ones(nz,1)*cfl_mu').*(i_flux - circshift(i_flux,-1)); %must be more efficient
    %Substep #2: Implicit material terms
    for j=1:ntheta
    intensity_new(:,j) = 1./(1 + c*dt.*X_tot*time_centering).*(intensity_substep(:,j) + ...
        c*dt.*X_tot.*(source_function - (1-time_centering).*intensity_substep(:,j)));
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


%Fully implicit scheme, donor cell interpolation
%---------------------
%At each timestep, solve Ax=b
% b = I^n, x=I^n+1
% the linear operator A depends on the interpolation scheme, the time
% centering, and the source terms

% for i=1:nt
%     intensity = intensity_new; 
%     %Reapply boundary conditions
%     intensity(1,1:ntheta/2) = 1.0;
%     intensity(nz,ntheta/2+1:ntheta) = 1.0;
%     for j=1:ntheta
%         b = intensity(:,j) + c*dt.*X_tot.*source_function;
%         %Construct A for donor cell method
%         %Main diagonal, I^n+1_i
%         %Assume time_centering=1 of fluxes
%         main_diag = (1 + abs(cfl_mu(j))).*ones(nz,1) - c*dt.*X_tot*time_centering;
%         right_diag = zeros(nz-1,1); 
%         left_diag = right_diag;
%         if mu(j) > 0
%             left_diag = -cfl_mu(j)*ones(nz-1,1);
%         else
%             right_diag = cfl_mu(j)*ones(nz-1,1);
%         end
%         A = diag(main_diag,0) + diag(left_diag,-1) + diag(right_diag,+1);
%         %Boundary conditions
%         if mu(j) > 0
%             A(1,1) = 1;
%         else 
%             A(nz,nz) = 1;
%         end
%         intensity_new(:,j) = A\b; 
%     end
%     time = dt*i
%         if ~mod(i,1)
%         %Plot each angular intensity (recall, ntheta must be even)
%         for j=1:ntheta
%         hi = subplot(2,1,j); 
%         plot(hi,zz,intensity(:,j,1),'-o');
%         x_label = sprintf('mu = %f',mu(j));
%         time_title = sprintf('t = %f (s)',time);
%         legend(method);
%         xlabel(x_label);
%         title(time_title); 
%         end
%         pause(1.0);
%     end
% end

%Fully implicit method, van Leer interpolation
% Nonfunctional-- gave up since Jacobian is not analytic
%----------------------------------
% MAX_ITER = 10000;
% delta_nr = 1e-5; %residual tolerance per Newton Raphson step
% intensity_iter = zeros(nz,1); %next iteration guess (per ray)
% for i=1:nt
%     intensity = intensity_new; 
%     %Reapply boundary conditions
%     intensity(1,1:ntheta/2) = 1.0;
%     intensity(nz,ntheta/2+1:ntheta) = 1.0;
%     for j=1:ntheta
%         b = intensity(:,j) + c*dt.*X_tot.*source_function;
%         %van Leer interpolation scheme is nonlinear
%         %must use iterative scheme to solve system
%         for k=1:MAX_ITER %Use previous time step as guess
%             %Evaluate function at current iteration
%             i_flux = upwind_interpolate(intensity(:,j),'van Leer',mu(j),dt,dz,c); 
%             function_iter = -b + abs(cfl_mu(j)).*(i_flux - circshift(i_flux,-1));
%             %Recheck the limiter condition at each point
%             vLslopes_condition = ((circshift(intensity(:,j),-1) - intensity(:,j)).* ...
%                 (intensity(:,j) - circshift(intensity(:,j),+1)) > 0);
%             %Construct the derivatives of the vL slopes
%             for l=2:nz-1
%                 left_deriv(l) = 2*(intensity(l,j) - intensity(l+1,j))/(intensity(l+1,j) - ...
%                     intensity(l-1,j)) - (intensity(l+
%              
%             %slope centered at i, derivative at i-1
%             
%             %slope centered at i, derivative at i
%             
%             %slope centered at i, derivative at i+1
%             %Construct Jacobian for van Leer interpolation
%             if mu(j) > 0
%                 main_diag = ones(nz,1) + abs(cfl_mu(j))*((1 - abs(cfl_mu(j)))/2.* ...
%                     circshift(vLslopes_condition,+1).*DERIVATIVE) + abs(cfl_mu(j)).*
%                 
%                 .*ones(nz,1) - c*dt.*X_tot*time_centering;
%                 
%                 J = diag(main_diag,0); 
%             end
%             step = -J\function_iter; %solve for NR step
%             intensity(:,j) = intesnity(:,j) + step; 
%             i_flux = upwind_interpolate(intensity(:,j),'van Leer',mu(j),dt,dz,c); 
%             residual = -b + abs(cfl_mu(j)).*(i_flux - circshift(i_flux,-1));
%             if max(residual) <= delta_nr
%                 break;
%             end
%         end
%     end
%     time = dt*i
%         if ~mod(i,1)
%         %Plot each angular intensity (recall, ntheta must be even)
%         for j=1:ntheta
%         hi = subplot(2,1,j); 
%         plot(hi,zz,intensity(:,j,1),'-o');
%         x_label = sprintf('mu = %f',mu(j));
%         time_title = sprintf('t = %f (s)',time);
%         legend(method);
%         xlabel(x_label);
%         title(time_title); 
%         end
%         pause(1.0);
%     end
% end
% should implement integration over domain of dependence

