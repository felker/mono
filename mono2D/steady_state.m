%Solve steady state integral equation, eq 17 in Jiang14 with dI/dt=0
%solve at one arbitrary spatial point, nx/2,ny/2
%np is number of discretization levels on [0,2pi]
function [] = steady_state()
np = 40; 
sigma_a = 100;
steady_temp = 0.983881934648789;%final temperature given
[direction_cosines,point_weights] = lgwt(np,-1,1);
mu_x = direction_cosines;
mu_z = 0.881917103688197;
mu_y = sqrt(1-direction_cosines.^2-mu_z^2);
C =10;
vx = 0.3*C;
%Use secant method in multidimensions

function output = func2minimize(I) %from R^np --> R^np
    J = point_weights'*I; 
    Kxx = (direction_cosines(:,1).^2.*point_weights)'*I; 
    output = (C*sigma_a*(steady_temp^4/(4*pi) - I) + 3*mu_x*vx*steady_temp^4/(4*pi)*sigma_a + ...
        mu_x*vx*sigma_a.*I - sigma_a*vx^2/C*(J + Kxx));
end
% end function 
x0 = 0.1*ones(np,1);
dx = 0.000001;
[I_analytic] = secantm(x0,dx,@func2minimize); 
plot(I_analytic.*mu_x,I_analytic.*mu_y,'--b',I_analytic.*mu_x,-I_analytic.*mu_y,'--b');


np = 40; 
sigma_a = 100;
steady_temp = 0.983881934648789;;     %final temperature given
[direction_cosines,point_weights] = lgwt(np,-1,1);
mu_x = direction_cosines;
mu_z = 0.3333;
mu_y = sqrt(1-direction_cosines.^2-mu_z^2);
C =10;
vx = 0.3*C;
%Use secant method in multidimensions

x0 = 0.1*ones(np,1);
dx = 0.000001;
[I_analytic] = secantm(x0,dx,@func2minimize); 
plot(I_analytic.*mu_x,I_analytic.*mu_y,'--k',I_analytic.*mu_x,-I_analytic.*mu_y,'--k');
end
