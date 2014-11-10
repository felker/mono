function [J] = jacobFD(f,x,delx) 
% Calculates the Jacobian of the
% system of non-linear equations:
% f(x) = 0, through finite differences.
% The Jacobian is built by columns
[m n] = size(x); 
for j = 1:m 
 xx = x; 
 xx(j) = x(j) + delx; 
 J(:,j) = (f(xx)-f(x))/delx; 
end; 
% end function