%2D Time-Dependent RT solver
%Based on Jiang14 algorithm
%Monotonic Upwind Interpolation
%Mixed frame to O(v/c)

% TO-DO:
%RENAME RHO_A RHO_S TO SIGMAS
%NEED BETTER WAY OF PASSING READ-ONLY INFO TO OUTPUT FUNCTIONS

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
dt = 0.0025;
nt = 940;
%Upwind monotonic interpolation scheme
method = 'van Leer'; 
time_centering = 0; %explicit = 0, implicit = 1, CN=1/2
delta_c = 1e-14; %error tolerance for NR convergence

%------------------------ BOUNDARY CONDITIONS  ------------------ %
%%% van Leer interpolation requires two bc 
num_ghost = 2; 
is = 1+num_ghost; ie = nx+num_ghost;
js= 1+num_ghost; je= ny+num_ghost;
%number of physical cells
nx_r= nx;
ny_r = ny;
nx = nx+2*num_ghost;
ny = ny+2*num_ghost;
%------------------------ BUILD MESH ------------------ %
%Fixed Background Fluid Velocity 
v = zeros(nx,ny,2); 
%Fixed fluid temperature
temp = zeros(nx,ny); 
%Gas density
density = zeros(nx,ny);

%Radiation Angular Discretization, S_n discrete ordinates. 
[na, mu, pw, lw] = angular_quad2D(ntheta);

%Radiation Spatial Discretization
%Radiation points are centered on the fluid cells
xx=linspace(0,lx,nx_r)';
yy=linspace(0,ly,ny_r)';

%Monochromatic specific intensity, boundary conditions at 2,nz-1 
intensity = zeros(nx,ny,na); 

%--------------- JIANG14 variables, dimensionless constants --------------%
a_r =  5.670373e-8; %stefan boltzman constant
adiabatic = 5/3; %gamma=adiabatic index
R = 1; %8.311; %ideal gas constant

%characteristic values (arbitrary?)
a_0 = 0.1; %characteristic velocity
T_0 = 1.0; %characteristic temperature
P_0 = a_r;%1.0; %characteristic gas pressure

C = c/a_0; %dimensionlesss speed of light
P = a_r*T_0^4/P_0;  %measure of relative importance of rad and gas pressure

%------------------------ PROBLEM SETUP ------------------ %
%For 2D and verification of Jiang14, we no longer use the source function
%or extinction probability 
%Absorption opacities
rho_a = ones(nx,ny);
%Scattering opacities
rho_s = zeros(nx,ny);
%Fluid density, temperature
density = ones(nx,ny);
temp = ones(nx,ny); 

v(:,:,1) = 0.3*C; 
intensity(:,:,:) = 1.00000/(4*pi); 

%------------------------ PRE-TIMESTEPPING SETUP ------------------ %
%Calculate Radiation CFL numbers
cfl_mu = C*dt*abs(mu)*[1/dx 1/dy]';
%set dt to have max cfl of 0.4
dt = dt*1.0/max(cfl_mu);
%Recalculate radiation CFL
cfl_mu = C*dt*abs(mu)*[1/dx 1/dy]';
assert(min(abs(cfl_mu) <= ones(na,1))); 

dt = 0.0025;

%------------------------ VELOCITY TERMS ------------------------------ %
[nv, vvnn, vCsquare, vsquare, absV] = update_velocity_terms(v,mu,C);
GasMomentum = zeros(nx,ny,2); 
GasKE = zeros(nx,ny);
GasTE = R/(adiabatic-1)*density.*temp; 
for i=1:2
    GasMomentum(:,:,i) = (density.*v(:,:,i));
    GasKE(:,:) = GasKE(:,:) + 0.5*v(:,:,i).^2;
end
%Lorentz factor
gamma = 1./(sqrt(1-vCsquare)); 

%------------------------ OUTPUT VARIABLES------------------------------ %
output_interval = 20; 
num_output = 8; %number of data to output
num_pts = nt/output_interval; 
time_out = dt*linspace(0,nt+output_interval,num_pts+1); %extra pt for final step
y_out = zeros(num_pts+1,num_output);
%Explicit-Implicit operator splitting scheme
%-------------------------------------------
for i=0:nt
    time = dt*i %#ok<NOPTS> 
    %Substep 0: Time series output, boundary condition, moments update
    %Update moments
    [J,H,K,rad_energy,rad_flux,rad_pressure] = update_moments(intensity,mu,pw,c);
    photon_momentum = P*rad_flux/C;
    if ~mod(i,output_interval)
        temp(1,1)
        [y_out] = time_series_output(y_out,time_out,rad_energy,rad_flux,rad_pressure,v,nx,ny,C,GasMomentum,GasKE,GasTE,photon_momentum);
    end

    %Box periodic boundary conditions
    for j=1:num_ghost
        intensity(:,js-j,1:na/2) = intensity(:,je+1-j,1:na/2); 
        intensity(:,je+j,na/2:na) = intensity(:,js-1+j,na/2:na); 
        %mu_x arent easily divided into positive and negative
        intensity(is-j,:,1:3) = intensity(ie+1-j,:,1:3); 
        intensity(is-j,:,7:9) = intensity(ie+1-j,:,7:9); 
        intensity(ie+j,:,4:6) = intensity(is-1+j,:,4:6); 
        intensity(ie+j,:,10:12) = intensity(is-1+j,:,10:12); 
    end
    
     %Substep #1: Explicitly advance transport term
    net_flux = zeros(nx,ny,na);
    for j=1:na %do all nx, ny at once
        %Split the transport term into diffusion and advection terms
        %Photon diffusion at nearly the speed of light
        tau = (10*(dx)*(rho_a+rho_s)).^2; %Eq 14 in Jiang14
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
       
        %Advection with the fluid only if there are material terms
        advection_term = (3*nv(:,:,j).*J).*(abs(rho_a+rho_s)>0);
        i_flux = upwind_interpolate2D_decoupledV2(mu(j,1)*(intensity(:,:,j) - advection_term/C),method,dt,dx,alpha*C*sign(mu(j,1)),is,ie+1,js,je+1,1);
        net_flux(is:ie,js:je,j) = dt*C/dx*(i_flux(is+1:ie+1,js:je) - i_flux(is:ie,js:je));
        %should also turn off advection velocity if the projection along
        %axis is zero in rare cases
        %Also, should take into account velocity gradient between cells
        advection_term = (mu(j,1)^2*3*J).*(abs(rho_a+rho_s)>0);
        velocity_term = v(:,:,1) + mu(j,2)*v(:,:,2)/mu(j,1);
        i_flux = upwind_interpolate2D_decoupledV2(advection_term, method, dt, dx, velocity_term,is,ie+1,js,je+1,1);
        net_flux(is:ie,js:je,j) = net_flux(is:ie,js:je,j) + ...
            dt/dx*0.5*((velocity_term(is+1:ie+1,js:je)+velocity_term(is:ie,js:je))*i_flux(is+1:ie+1,js:je) -...
            (velocity_term(is-1:ie-1,js:je)+velocity_term(is:ie,js:je))*i_flux(is:ie,js:je));
    end %end of ray trace loop, x-direction
    
    for j=1:na %do all nx, ny at once
        %Split the transport term into diffusion and advection terms
        %Photon diffusion at nearly the speed of light
        tau = (10*(dy)*(rho_a+rho_s)).^2; %Eq 14 in Jiang14
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
         %Advection with the fluid only if there are material terms
        advection_term = (3*nv(:,:,j).*J).*(abs(rho_a+rho_s)>0);
        i_flux = upwind_interpolate2D_decoupledV2(mu(j,2)*(intensity(:,:,j) - advection_term/C),method,dt,dy,alpha*C*sign(mu(j,2)),is,ie+1,js,je+1,2);
        net_flux(is:ie,js:je,j) = dt*C/dy*(i_flux(is:ie,js+1:je+1) - i_flux(is:ie,js:je));
        %should also turn off advection velocity if the projection along
        %axis is zero in rare cases
        %Also, should take into account velocity gradient between cells
        advection_term = (mu(j,2)^2*3*J).*(abs(rho_a+rho_s)>0);
        velocity_term = v(:,:,2) + mu(j,1)*v(:,:,1)/mu(j,2);
        i_flux = upwind_interpolate2D_decoupledV2(advection_term, method, dt, dy, velocity_term,is,ie+1,js,je+1,2);
        net_flux(is:ie,js:je,j) = net_flux(is:ie,js:je,j) + ...
            dt/dy*0.5*((velocity_term(is:ie,js+1:je+1)+velocity_term(is:ie,js:je))*i_flux(is:ie,js+1:je+1) -...
            (velocity_term(is:ie,js-1:je-1)+velocity_term(is:ie,js:je))*i_flux(is:ie,js:je));
    end    %end of ray trace loop, y-direction
    intensity = intensity - net_flux; 
        [J,H,K,rad_energy,rad_flux,rad_pressure] = update_moments(intensity,mu,pw,c);
        % is the source term update the only meaningful part of this test?
    assert(max(max(max(net_flux))) == 0);
    
    %Substep 1.2: Estimate flow velocity at t= n+1/2
    new_velocity = zeros(nx,ny,2);
    for k=1:nx %rays must be solved together
        for l=1:ny
            for m=1:2
            new_velocity(k,l,m) = 1/(density(k,l) + 0.5*dt*P*(rho_a(k,l) + rho_s(k,l))*C/P*density(k,l) + 0.5*dt*P*(rho_a(k,l) + rho_s(k,l))*(rad_energy(k,l) + rad_pressure(k,l,m,m))/C)*...
                (density(k,l)*v(k,l,m) + 0.5*dt*P*(rho_a(k,l) + rho_s(k,l))*(C/P*density(k,l)*v(k,l,m) + rad_flux(k,l,m)));
            end
        end
    end
    old_velocity = v;
    v = new_velocity;
    [nv, vvnn, vCsquare, vsquare, absV] = update_velocity_terms(v,mu,C);
    gamma = 1./(sqrt(1-vCsquare)); 

    %Substep #2: Implicitly advance the absorption source terms at each
    %cell, for all rays in a cell 
    for k=1:nx 
        for l=1:ny
            I_old = squeeze(intensity(k,l,:)); 
            temp_old = temp(k,l);
            %Temperature is fixed in this method
%             for m=1:na %rays no longer need to be solved together
%                 intensity(k,l,m) = 1/(1+C*dt*rho_a(k,l)*gamma(k,l)*(1-nv(k,l,m)/C))*(I_old(m) + ...
%                     C*dt*rho_a(k,l)*gamma(k,l)^-3*(1-nv(k,l,m)/C)^-3*temp(k,l)^4/(4*pi));        
%             end
            %Temperature is implicitly updated in the following method
            
            %I^n+1_i only depends on T^n+1, I^n_i, and n_i*v
            
            A = zeros(na,1);
            for r=1:na
                A(r) = 1 + C*dt*rho_a(k,l)*gamma(k,l)*(1-nv(k,l,r)/C); 
            end
            
            %prefactor on RHS
            J1 = -(adiabatic-1)*C*dt*rho_a(k,l)*P/(density(k,l)*R);
            coef1=0;
            coef2=0;
            %coupling to rays
            for r = 1:na
                coef1 = coef1 +gamma(k,l)^-3*(1-nv(k,l,r)/C)^-2*pw(r)- gamma(k,l)^-2*(1-nv(k,l,r)/C)^-3*pw(r)/A(r)*C*dt*rho_a(k,l)*(...
                    1-2/C*nv(k,l,r) + (C^-2)*vvnn(k,l,r));
                %does the vvnn term make sense?
                coef2 = coef2 - 4*pi*gamma(k,l)*pw(r)/A(r)*I_old(r)*(1-2/C*nv(k,l,r) + (C^-2)*vvnn(k,l,r));
            end           
            Tcoef = -1; 
            T4coef = J1*coef1; 
            Tconst = temp_old + J1*coef2;            
            Tmax = (0.5*density(k,l)*vsquare(k,l) + P*rad_energy(k,l))*(adiabatic -1.0)/(density(k,l)*R) + temp_old; 

            if abs(T4coef) < 1e-16
                temp(k,l) = -Tconst/Tcoef;
            else
                Tmin = max(-0.25*Tcoef/T4coef,0.0).^(1./3); 
                fmin = T4coef*Tmin^4.0 + Tcoef*Tmin + Tconst;
                fmax = T4coef*Tmax^4.0 + Tcoef*Tmax + Tconst;
                if (T4coef * fmin > 0.0)
                    error('No solution in this case');
                else
                    %deal function and derivative
                    fn_handle = @(t,coef1,coef2,coef3,coef4) deal(coef1*t^4 + coef2*t + coef3,4.0*coef1*t^3 + coef2); 
                    temp(k,l) = rtsafe_orig(fn_handle,Tmin,Tmax,delta_c, T4coef, Tcoef, Tconst, 0.0);
                end
            end
             
            %Explicitly update all rays 
            intensity(k,l,:) = 1./A(:).*(I_old + C*dt*rho_a(k,l)*temp(k,l)^4./(4*pi*gamma(k,l)^3*(1-squeeze(nv(k,l,:)/C)).^3));
            
            %Revert to velocity at the beginning of the timestep
            v(k,l,:) = old_velocity(k,l,:);
            
            %Update gas quantities 
            dGasMomentum = zeros(2,1);
            for r=1:2 %4pi factor?? why isnt that in the paper??
                dGasMomentum(r) = -P/C*4*pi*((mu(:,r).*pw(:))'*squeeze(intensity(k,l,:)) - (mu(:,r).*pw(:))'*I_old);
                GasMomentum(k,l,r) = GasMomentum(k,l,r) + dGasMomentum(r);
            end
            dGasKE =((density(k,l)*squeeze(v(k,l,:)) + dGasMomentum)'*(density(k,l)*squeeze(v(k,l,:)) + dGasMomentum) - ...
                (density(k,l)*squeeze(v(k,l,:)))'*(density(k,l)*squeeze(v(k,l,:))))/(2*density(k,l));
            GasKE(k,l) = GasKE(k,l) + dGasKE;  
            %Explicit change of temperature when using fixed temperature in
%update
%             temp(k,l) = temp_old - (density(k,l)*R/(adiabatic-1))^-1*(P*4*pi*(...
%                 pw(:)'*squeeze(intensity(k,l,:)) - pw(:)'*I_old) + dGasKE); 
            %change in gas internal energy
            dGasTE = (density(k,l)*R/(adiabatic-1)*(temp(k,l) - temp_old));
            GasTE(k,l) = GasTE(k,l) + dGasTE;  
         end
    end
    temp(1,1)
    %Update fluid velocity
    v(:,:,1) = GasMomentum(:,:,1)./density(:,:);
    v(:,:,2) = GasMomentum(:,:,2)./density(:,:);
    %Update periodic boundary condtions, velocity
    for j=1:num_ghost
        v(:,js-j,:) = v(:,je+1-j,:);
        v(:,je+j,:) = v(:,js-1+j,:);
        v(is-j,:,:) = v(ie+1-j,:,:);
        v(ie+j,:,:) = v(is-1+j,:,:);        
    end
    [nv, vvnn, vCsquare, vsquare, absV] = update_velocity_terms(v,mu,C);
  %      Substep #3: Implicitly advance the scattering source terms
    
%------------------------ NON-TIME SERIES OUTPUT ------------------ %
 %   if ~mod(i,output_interval)
  %      static_output(); 
  %  end
end


