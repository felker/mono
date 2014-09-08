function [num_rays,direction_cosines,point_weights,level_weights] = angular_quad3D(N)
%Input:
% N: Order of the quadrature = # of polar angles [0,pi] 
% only defined for N <= 12
%Output:
% num_rays
% direction_cosines 
% point_weights
% level_weights
%References: B.G. Carlson 1963, Bruls 1999
%For 2D, we assume the physical domain is homogeneous in the z-axis.
%Therefore, we only do the positive hemisphere about the midplane, where
%the z-axis is the polar axis. 
%In 2D,
num_rays = N*(N+2)/2; %Up to 168
num_rays_per_octant = num_rays/4;
direction_cosines = zeros(num_rays,3);
point_weights = zeros(num_rays,1);
level_weights = zeros(N/2,1); 

%This is the subset of [-1,1] that the direction_cosines are taken from in each
%dimension. Symmetry in each octant => only need one array
%Antisymmetric about zero
ordinates = zeros(1,N/2); 

%Choice of first cosine is arbitary. In analogy to Gaussian quadrature:
ordinates(1) = sqrt(1/(3*(N-1))); 
%In 3D space,
delta_cos =(1-3*ordinates(1)^2)* 2/(N-2) 

for i=2:N/2
    ordinates(i) = sqrt(ordinates(1)^2 + (i-1)*delta_cos);
end

%Derive weights for integration. See Bruls99 appendix
W = sqrt( 4/(3*(N-1)) + ([1:N/2] - 1)*2/(N-1));
%Compute level weights
level_weights(1) = W(1);
for i=2:N/2-1
    level_weights(i) = W(i) - W(i-1);
end

%following ATHENA, we correct the last level weight to renormalize...not
%sure why this happens
level_weights(N/2) = 1.0 - sum(level_weights); 
%permuation matrix
pl = zeros(N/2,3);
pmat = zeros(N/2,N/2); 
plab = zeros(num_rays,1);

%To select valid rays for the quadrature, we follow the index selection
%rule:
%sigma = i + j + k, indices of the direction cosines must equal ntheata/2 +2
%for normalization in 3D space
ray_count = 0;
np = 0;
ip = 0;
for i=1:N/2
   for j=1:(N/2+1-i)
       k = N/2+2-i-j;
       %how do we order the rays in each octant?
       ray_count = ray_count +1;
       direction_cosines(ray_count,1) = ordinates(i);
       direction_cosines(ray_count,2) = ordinates(j);
       direction_cosines(ray_count,3) = ordinates(k);
       if N <= 12 %system of equations is lin. ind. iff N<=12
       ip = permutation(i,j,k,pl,np);
       if ip == -1 %ray indices havent been loaded yet
           np = np+1;
           pl(np,1) = i;
           pl(np,2) = j;
           pl(np,3) = k;
           pmat(i,np) = pmat(i,np) + 1.0;
           plab(ray_count) = np;
       else %these indices have been loaded 
          pmat(i,ip) = pmat(i,ip) + 1.0;
          plab(ray_count) = ip; 
       end
       end
   end
end
assert(ray_count == num_rays_per_octant);

%Symmetry: reflect across second axis 
%(must be a better way of doing this)
direction_cosines(ray_count+1:2*ray_count,1) = -direction_cosines(1:ray_count,1);
direction_cosines(ray_count+1:2*ray_count,2) = direction_cosines(1:ray_count,2);
direction_cosines(ray_count+1:2*ray_count,3) = direction_cosines(1:ray_count,3);
ray_count = ray_count*2;
%Symmetry: reflect across first axis
direction_cosines(ray_count+1:2*ray_count,1) = direction_cosines(1:ray_count,1);
direction_cosines(ray_count+1:2*ray_count,2) = -direction_cosines(1:ray_count,2);
direction_cosines(ray_count+1:2*ray_count,3) = direction_cosines(1:ray_count,3);
ray_count = ray_count*2;
%Material symmetry: do not reflect across final axis

%test ray normalization
%some of these norms have a deviation from 1.0 in the 15th decimal place
assert(all((abs(sum((direction_cosines.*direction_cosines),2) - ones(num_rays,1)) < 1e-12)));

%test output in unit sphere
%quiver3(zeros(ray_count,1),zeros(ray_count,1),zeros(ray_count,1), ...
 %   direction_cosines(1:ray_count,1),direction_cosines(1:ray_count,2), ...
  %  direction_cosines(1:ray_count,3)); 

%Use Bruls method to calculate families of weights. Solve system of eqs
inv(pmat(1:N/2-1,1:N/2-1))
wpf = pmat(1:N/2-1,1:N/2-1)\level_weights(1:N/2-1);
for i=1:num_rays_per_octant %each quadrant, normalize %Is this ordered correctly?
    point_weights(i) = wpf(plab(i))*0.125;
    point_weights(i+num_rays_per_octant) = wpf(plab(i))*0.125; 
    point_weights(i+2*num_rays_per_octant) = wpf(plab(i))*0.25; 
    point_weights(i+3*num_rays_per_octant) = wpf(plab(i))*0.25; 
end

end