function [nv, vvnn, vCsquare, vsquare, absV] = update_velocity_terms(v,mu,C )
%UPDATE_VELOCITY_TERMS
[nx,ny,ndim] = size(v); 
[na, ndim] = size(mu);

%Calculate dot product of local fluid velocity and radiation rays
nv = zeros(nx,ny,na);
%Tensor component in the O(v/c) absorption quadrature terms
vvnn = zeros(nx,ny,na);
%Calculate loval velocity magnitude
vCsquare = zeros(nx,ny);
vsquare = zeros(nx,ny);
absV = zeros(nx,ny);
for i=1:nx
    for j=1:ny
        vsquare(i,j) = (v(i,j,1)^2 + v(i,j,2)^2);         
        vCsquare(i,j) = (vsquare(i,j))/C^2; 
        absV(i,j) = sqrt(vsquare(i,j));
        for k=1:na
            nv(i,j,k) = v(i,j,1)*mu(k,1) + v(i,j,2)*mu(k,2);
            vvnn(i,j,k) = v(i,j,1)^2*mu(k,1)^2 + v(i,j,2)^2*mu(k,2)^2 + ...
                    2*v(i,j,1)*mu(k,1)*v(i,j,2)*mu(k,2);
        end
    end
end

end

