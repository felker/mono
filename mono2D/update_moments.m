function [J,H,K,rad_energy,rad_flux,rad_pressure] = update_moments(intensity,mu,pw,c)
%UPDATE_MOMENTS Given the current intensity, calculate the angular moments
% of the radiation field 
[nx,ny,na] = size(intensity);
%Eddington moments
J = zeros(nx,ny);
H = zeros(nx,ny,2);
K = zeros(nx,ny,2,2);
for i=1:nx
    for j=1:ny
        for k=1:na
            J(i,j) = J(i,j) + intensity(i,j,k)*pw(k);
            for l=1:2
                H(i,j,l) = H(i,j,l) + intensity(i,j,k)*pw(k)*mu(k,l); 
                for m=1:2
                    K(i,j,l,m) = K(i,j,l,m) + intensity(i,j,k)*pw(k)*mu(k,l)*mu(k,m); 
                end
            end
        end
    end
end
rad_energy = 4*pi*J;
rad_flux = 4*pi*c*H; %c or C?
rad_pressure = 4*pi*K;
end

