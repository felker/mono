function i_flux = upwind_interpolate2D_decoupled(intensity, method, dt, ds, v, nx,ny,dir)
%UPWIND_INTERPOLATE1D_decoupled For finite volume, pick one of several methods for
%advection step into cell. Does both directions at once
% Does only one ray at a time
if size(intensity) ~= [nx,ny]
    error('incorrect passing to interpolation function');
end
    
%First index is LHS boundary
%Second index is top boundary 
i_flux = zeros(nx,ny);

switch method
    case 'donor'
    case 'van Leer'
        %Monotonized van Leer slopes
        vLslopes = zeros(nx,ny);
        if dir==1    
            for j=2:ny-1
                for i=2:nx-1 
                    if ((intensity(i+1,j) - intensity(i,j))*(intensity(i,j) - ...
                        intensity(i-1,j)) > 0)
                        vLslopes(i,j) = 2*(intensity(i+1,j) - intensity(i,j))*(intensity(i,j) - ...
                            intensity(i-1,j))/ (intensity(i+1,j) - intensity(i-1,j));
                    else
                        vLslopes(i,j) = 0; 
                    end
                end 
                vL_circ = circshift(vLslopes,[+1,0]);  
                i_circ = circshift(intensity,[+1,0]);  
                for i=1:nx
                    vel = v(i,j);
                    %vel = 0.5*(lhs_v(i-1)+ v(i,j));
                    if vel > 0
                        i_flux(i,j) = i_circ(i,j) + (1-vel*dt/ds)*vL_circ(i,j)/2;
                    else %interface velocity is negative
                        i_flux(i,j) = intensity(i,j) - (1+vel*dt/ds)*vLslopes(i,j)/2 ;  
                    end
                end
            end 
        elseif dir==2
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
                vL_circ = circshift(vLslopes,[0,+1]);  
                i_circ = circshift(intensity,[0,+1]); 
                for j=1:ny
                    vel = v(i,j);
                    %vel = 0.5*(lhs_v(i-1)+ v(i,j));
                    if vel > 0
                        i_flux(i,j) = i_circ(i,j) + (1-vel*dt/ds)*vL_circ(i,j)/2;
                    else %interface velocity is negative
                        i_flux(i,j) = intensity(i,j) - (1+vel*dt/ds)*vLslopes(i,j)/2 ;  
                    end
                end
            end
        end
    case 'CW' %colella and woodward
    otherwise
        warning('unexpected interpolation scheme')
end

