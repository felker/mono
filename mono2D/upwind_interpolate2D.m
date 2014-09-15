function i_flux = upwind_interpolate2D(intensity, method, mu, dt, dx, dy, c)
%UPWIND_INTERPOLATE2D For finite volume, pick one of several methods for
%advection step into cell. Does both directions at once
% Does only one ray at a time
[nx,ny] = size(intensity); 
%First index is LHS boundary
%Second index is top boundary 
i_flux = zeros(nx,ny,2);

switch method
    case 'donor'
        %Donor Cell method, first order accurate
        if mu(1) > 0
            i_flux(:,:,1) = circshift(intensity(:,:),[1,0]); 
        elseif mu(1) < 0
            i_flux(:,:,1) = intensity(:,:);            
        end
        if mu(2) > 0
            i_flux(:,:,2) = circshift(intensity(:,:),[0,1]); 
        elseif mu(2) < 0
            i_flux(:,:,2) = intensity(:,:);            
        end
    case 'van Leer'
        %Monotonized van Leer slopes
        vLslopes = zeros(nx,ny);
        for i=2:ny-1
            for j=2:nx-1
                if ((intensity(j+1,i) - intensity(j,i))*(intensity(j,i) - ...
                        intensity(j-1,i)) > 0)
                    vLslopes(j,i) = 2*(intensity(j+1,i) - intensity(j,i))*(intensity(j,i) - ...
                        intensity(j-1,i))/ (intensity(j+1,i) - intensity(j-1,i));
                else
                    vLslopes(j,i) = 0; 
                end
            end
        end
        if mu(1) > 0
            i_flux(:,:,1) = circshift(intensity(:,:),[+1, 0]) + ...
                     (1-c*mu(1)*dt/dx)*circshift(vLslopes,[1,0])/2; 
        elseif mu(1) < 0
            i_flux(:,:,1) = intensity(:,:) - (1+c*mu(1)*dt/dx)*vLslopes/2 ;            
        end
        
        vLslopes = zeros(nx,ny);
        for i=2:nx-1
            for j=2:ny-1
                if ((intensity(i,j+1) - intensity(i,j))*(intensity(i,j) - ...
                        intensity(i,j-1)) > 0)
                    vLslopes(i,j) = 2*(intensity(i,j+1) - intensity(i,j))*(intensity(i,j) - ...
                        intensity(i,j-1))/ (intensity(i,j+1) - intensity(i,j-1));
                else
                    vLslopes(i,j) = 0; 
                end
            end
        end
        if mu(2) > 0
            i_flux(:,:,2) = circshift(intensity(:,:),[0,+1]) + ...
                     (1-c*mu(2)*dt/dy)*circshift(vLslopes,[0,1])/2; 
        elseif mu(2) < 0
            i_flux(:,:,2) = intensity(:,:) - (1+c*mu(2)*dt/dy)*vLslopes/2 ;            
        end    
    case 'CW' %colella and woodward
    otherwise
        warning('unexpected interpolation scheme')
end

