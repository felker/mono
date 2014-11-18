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
xx=linspace(0,lx,nx)';
yy=linspace(0,ly,ny)';

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
dt = dt*0.8/max(cfl_mu);
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

    %Reapply boundary conditions
%      intensity(nx,:,:) = 0.0; 
%      intensity(:,ny,:) = 0.0;
%      intensity(:,1,:) = 0.0;
%      intensity(1,:,:) = 0.0;

    %intensity(1,2:ny-1,:) = 0.0; %all first row elements (physically, left bndry)

    %Setup for 2D crossing beams in Jiang14 and Davis12
    %Periodic boundary conditions on y-boundaries, only for appropriate
    %rays
    %intensity(:,1,1:na/2) = intensity(:,ny-1,1:na/2); 
    %intensity(:,ny,na/2:na) = intensity(:,2,na/2:na); 
    %Inject beam at center of LHS boundary with mu_x=0.333, -0.3333
    %intensity(1,ny/2,1) = 1.0;
    %intensity(1,ny/2,7) = 1.0;
    
    %Box periodic boundary conditions
    intensity(:,1,1:na/2) = intensity(:,ny-1,1:na/2); 
    intensity(:,ny,na/2:na) = intensity(:,2,na/2:na); 
    %mu_x arent easily divided into positive and negative
    intensity(1,:,1:3) = intensity(nx-1,:,1:3); 
    intensity(1,:,7:9) = intensity(nx-1,:,7:9); 
    intensity(nx,:,4:6) = intensity(2,:,4:6); 
    intensity(nx,:,10:12) = intensity(2,:,10:12); 
    
    %Substep #1: Explicitly advance transport term
    net_flux = zeros(nx,ny,na);
    for j=1:na %do all nx, ny at once
        %Split the transport term into diffusion and advection terms
        %Photon diffusion at nearly the speed of light
        %tau = (10*(dx/mu(j,1) + dy/mu(j,2))*(rho_a+rho_s)).^2; %Eq 14 in Jiang14
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
        %NEW open questions: 
        %- ATHENA and I place our mu terms on different sides of the
        %interpolation function: Does this matter in computaton of the van
        %Leer coefficients?
        
        % Same for the velocity projections, nv. On the advection term or the
        % advection velocity?
        
        % Can we do both directions at the same time?
        
        %where do we put the net flux? on the absorption term update?
                
        %KYLE FIX 11/9: divide nv*J by C !! this was in the paper, but I
        %missed it
        
        %do I need to change the - IV/C term in the first substep to only
        %happen if rho_a, rho_s are nonzero? This seems to make physical
        %sense. YES
        advection_term = (3*nv(:,:,j).*J).*(abs(rho_a+rho_s)>0);

        i_flux = upwind_interpolate2D(intensity(:,:,j) - advection_term/C,method,mu(j,:),dt,dx,dy,alpha*C);
%         i_flux2 = upwind_interpolate2D_decoupled(mu(j,1)*(intensity(:,:,j) - advection_term/C),method,dt,dx,alpha*C*sign(mu(j,1)),nx,ny,1);
%         intensity(:,:,j) - i_flux(:,:,1)

        A = circshift(i_flux(:,:,1),[-1 0]);
        B = circshift(i_flux(:,:,2),[0 -1]);
        net_flux(2:nx-1,2:ny-1,j) = (mu(j,1)*dt*C/dx*(i_flux(2:nx-1,2:ny-1,1) - A(2:nx-1,2:ny-1)) + ...
            mu(j,2)*dt*C/dy*(i_flux(2:nx-1,2:ny-1,2) - B(2:nx-1,2:ny-1)));
    
        %Advection with fluid at fluid velocity
        %Turn the advection term to zero in cells with zero material
        %interaction
        
        %should also turn off advection velocity if the projection along
        %axis is zero in rare cases
        %Also, should take into account velocity gradient between cells
        i_flux = upwind_interpolate2D(advection_term,method,mu(j,:),dt,dx,dy,absV);
        A = circshift(i_flux(:,:,1),[-1 0]);
        B = circshift(i_flux(:,:,2),[0 -1]);
        net_flux(2:nx-1,2:ny-1,j) = net_flux(2:nx-1,2:ny-1,j) + ...
            mu(j,1)*dt/dx*(i_flux(2:nx-1,2:ny-1,1) - A(2:nx-1,2:ny-1)) + ...
            mu(j,2)*dt/dy*(i_flux(2:nx-1,2:ny-1,2) - B(2:nx-1,2:ny-1)); 
    end %end of ray loop
    %Substep 1.2: Estimate flow velocity at t= n+1/2
    new_velocity = zeros(nx,ny,2);
    for k=2:nx-1 %rays must be solved together
        for l=2:ny-1
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
    for k=2:nx-1 %rays must be solved together
        for l=2:ny-1
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
            [nv, vvnn, vCsquare, vsquare, absV] = update_velocity_terms(v,mu,C);
        %Substep #3: Implicitly advance the scattering source terms
    
%------------------------ NON-TIME SERIES OUTPUT ------------------ %
    if ~mod(i,output_interval)
        %static_output(); 
    end
end


