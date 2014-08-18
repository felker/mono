function i_flux = upwind_interpolate(intensity, method, mu, dt, dz, c)
%UPWIND_INTERPOLATE For finite volume, pick one of several methods for
%advection step into cell
%   Detailed explanation goes here
%   Recall, intensity(nz,ntheta)
[nz,ntheta] = size(intensity); 
%Each index corresponds to LHS cell boundary
i_flux = zeros(nz,ntheta);

switch method
    case 'donor'
        %Donor Cell method, first order accurate
        for i=1:ntheta
            if mu(i) > 0
                i_flux(:,i) = circshift(intensity(:,i),[+1]); 
            elseif mu(i) < 0
                i_flux(:,i) = intensity(:,i);            
            end    
        end
    case 'van Leer'
        %Monotonized van Leer slopes
        vLslopes = zeros(nz,1);
        for i=1:ntheta
            for j=2:nz-1
                if ((intensity(j+1,i) - intensity(j,i))*(intensity(j,i) - ...
                        intensity(j-1,i)) > 0)
                    vLslopes(j) = 2*(intensity(j+1,i) - intensity(j,i))*(intensity(j,i) - ...
                        intensity(j-1,i))/ (intensity(j+1,i) - intensity(j-1,i));
                else
                    vLslopes(j) = 0; 
                end
            end
            if mu(i) > 0
                i_flux(:,i) = circshift(intensity(:,i),[+1]) + ...
                     (1-c*mu(i)*dt/dz)*circshift(vLslopes,1); 
            elseif mu(i) < 0
                i_flux(:,i) = intensity(:,i) - (1+c*mu(i)*dt/dz)*vLslopes/2 ;            
            end    
        end
    case 'CW' %colella and woodward
    otherwise
        warning('unexpected interpolation scheme')
end

