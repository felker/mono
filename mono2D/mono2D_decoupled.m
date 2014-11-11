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
dt = 0.002;
nt = 100;
%Upwind monotonic interpolation scheme
method = 'van Leer'; 
time_centering = 0; %explicit = 0, implicit = 1, CN=1/2
delta_c = 1e-10; %error tolerance for NR convergence

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
rho_a = 100*ones(nx,ny);
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
dt = dt*0.9/max(cfl_mu);
%Recalculate radiation CFL
cfl_mu = C*dt*abs(mu)*[1/dx 1/dy]';
assert(min(abs(cfl_mu) <= ones(na,1))); 

%------------------------ VELOCITY TERMS ------------------------------ %
[nv, vvnn, vCsquare, vsquare, absV] = update_velocity_terms(v,mu,C);
GasMomentum = zeros(nx,ny,2); 
GasKE = zeros(nx,ny);
%adiabatic -1? doesnt make sense with dGasTE formula
GasTE = R/(adiabatic-1)*density.*temp; 
for i=1:2
    GasMomentum(:,:,i) = (density.*v(:,:,i));
    GasKE(:,:) = GasKE(:,:) + 0.5*v(:,:,i).^2;
end

%------------------------ OUTPUT VARIABLES------------------------------ %
output_interval = 2; 
num_output = 6; %number of data to output
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
    %Substep 1.2: Estimate flow velocity at t= n+1/2
    new_velocity = zeros(nx,ny,2);
    for k=1:nx %rays must be solved together
        for l=1:ny
            for m=1:2
            vel_update = @(x) x - v(k,l,m) - 0.5*dt*P*(rho_a(k,l) + rho_s(k,l))*(C/P*v(k,l,m) + ...
                rad_flux(k,l,m)./density(k,l) - C/P*x - ...
                x./(C*density(k,l))*(rad_energy(k,l) + rad_pressure(k,l,m,m)));
            new_velocity(k,l,m) = fzero(vel_update, v(k,l,m));
            end
        end
    end
    old_velocity = v; 
    v = new_velocity;
    [nv, vvnn, vCsquare, vsquare, absV] = update_velocity_terms(v,mu,C);
    %Substep #2: Implicitly advance the absorption source terms at each
    %cell, for all rays in a cell 
    for k=1:nx %rays must be solved together
        for l=1:ny
            %Use Newton-Raphson. System is nonlinear in T^n+1. Use previous
            %values as initial guesses. Need to do implicitly because
            %thermalization is fast
            I_old = squeeze(intensity(k,l,:));
            temp_old = temp(k,l);
       
            %REDUCE TO 1D NR solve            
            A = zeros(na,1);
            B = zeros(na,1);
            D = zeros(na,1);
            A(na) = (1 +dt*rho_a(k,l)*C - dt*rho_a(k,l)*nv(k,l,na));
            for r=1:na-1
                A(r) = A(na)/(1 +dt*rho_a(k,l)*C - dt*rho_a(k,l)*nv(k,l,r)); 
                B(r) = 3.0*dt*rho_a(k,l)*(nv(k,l,r) - nv(k,l,na))/(4*pi*(1+dt*rho_a(k,l)*C*(1-nv(k,l,r)/C)));
                D(r) = (intensity(k,l,r) - intensity(k,l,na))/(1+dt*rho_a(k,l)*C*(1-nv(k,l,r)/C));
            end
            A(na) = 1.0;
            B(na) = 0.0;
            D(na) = 0.0; 
            
            %Match coefficients in hydro_to_rad.c
            Jfactor = P*(2+dt*rho_a(k,l)*C*(1+vCsquare(k,l)));
            J1 = - density(k,l)*R/(4*pi*(adiabatic-1)*Jfactor); 
            J2 = dt*rho_a(k,l)*P*C*(1+vCsquare(k,l))/(4*pi*Jfactor);
            J3 = 2*P*J(k,l)/Jfactor + density(k,l)*R*temp_old/(4*pi*(adiabatic-1)*Jfactor); 
            
            coef1 = 1 + dt*rho_a(k,l)*C*(1-nv(k,l,na)/C) + dt*rho_a(k,l)*((vsquare(k,l) + vvnn(k,l,r)).*A)'*pw/C;
            coef2 = dt*rho_a(k,l)*C*(1+3*nv(k,l,na)/C)/(4*pi) - dt*rho_a(k,l)*((vsquare(k,l) + vvnn(k,l,r)).*B)'*pw/C;
            coef3 = intensity(k,l,na) - dt*rho_a(k,l)*((vsquare(k,l) + vvnn(k,l,r)).*D)'*pw/C;
            
            Tcoef = coef1*J1;
            T4coef = coef1*J2;
            Tconst = coef1*J3;

            T4coef = T4coef - (coef1*B + coef2*A)'*pw;
            Tconst = Tconst - (coef1*D + coef3*A)'*pw;

            Tmax = (0.5*density(k,l)*vsquare(k,l) +rad_energy(k,l))*(adiabatic -1.0)/(density(k,l)*R) + temp_old; 

            %ATHENA equilibrium function, and derivative
            if abs(T4coef) < 1e-16
                temp(k,l) = -Tconst/Tcoef;
            else
                Tmin = max(-0.25*Tcoef/T4coef,0.0).^(1./3); 
                fmin = T4coef*Tmin^4.0 + Tcoef*Tmin + Tconst;
                fmax = T4coef*Tmax^4.0 + Tcoef*Tmax + Tconst;
                if (T4coef * fmin > 0.0)
                    error('No solution in this case');
                else
                    fn_handle = @(t,coef1,coef2,coef3,coef4) deal(coef1*t^4 + coef2*t + coef3,4.0*coef1*t^3 + coef2); 
                    temp(k,l) = rtsafe_orig(fn_handle,Tmin,Tmax,delta_c, T4coef, Tcoef, Tconst, 0.0);
                end
            end
            %Explicitly update all other rays 
            intensity(k,l,na) = (coef2*temp(k,l)^4 + coef3)/coef1 + net_flux(k,l,na);
            intensity(k,l,1:na-1) = A(1:na-1)*intensity(k,l,na) + B(1:na-1) + D(1:na-1) + squeeze(net_flux(k,l,1:na-1)); 
            v(k,l,:) = old_velocity(k,l,:); 
            %Update gas quantities as it absorbs internal energy
            %density?
            %change in gas internal energy
            dGasTE = (density(k,l)*R/adiabatic)*(temp(k,l) - temp_old);
            dGasMomentum = zeros(2,1);
            for r=1:2
                for n=1:na
                    dGasMomentum(r) = dGasMomentum(r) - P*(intensity(k,l,n) - I_old(n))*mu(n,r)*pw(n);                 
                end
                GasMomentum(k,l,r) = GasMomentum(k,l,r) + dGasMomentum(r);
            end
            %no change in momentum due to isotropy
            dGasKE =((density(k,l)*squeeze(v(k,l,:)) + dGasMomentum)'*(density(k,l)*squeeze(v(k,l,:)) + dGasMomentum) - ...
                (density(k,l)*squeeze(v(k,l,:)))'*(density(k,l)*squeeze(v(k,l,:))))/(2*density(k,l));
            GasKE(k,l) = GasKE(k,l) +dGasKE;  
            GasTE(k,l) = GasTE(k,l) +dGasTE;  
            
         end
    end
    %Change fluid velocity?
    v(:,:,1) = GasMomentum(:,:,1)./density(:,:);
    v(:,:,2) = GasMomentum(:,:,2)./density(:,:);
    %periodic boundary condtions, velocity
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


