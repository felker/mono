%Solve steady state integral equation, eq 17 in Jiang14 with dI/dt=0
%solve at one arbitrary spatial point, nx/2,ny/2
%np is number of discretization levels on [0,2pi]
function [] = steady_state()
sigma_a = 100;
steady_temp = 1.023395585899979;
%final temperature given
% [direction_cosines,point_weights] = lgwt(np,-radius,radius);
% mu_x = direction_cosines;
% mu_y = sqrt(radius^2-direction_cosines.^2);
[na, mu, point_weights, lw] = angular_quad2D(4); 
mu_x = mu(:,1);
mu_y = mu(:,2);
mu_remainder = sqrt(1-mu_x.^2 - mu_y.^2);
C =10;
vx = 2.955596777221964; %0.3*C;
%Use secant method in multidimensions
% np = np*2; 
% mu_x = [direction_cosines', direction_cosines']';
% mu_y = [sqrt(radius^2-direction_cosines'.^2), -sqrt(radius^2-direction_cosines'.^2)]';
% point_weights = [point_weights', point_weights']';
% point_weights = sum(point_weights)* point_weights;

function output = func2minimize(I) %from R^np --> R^np
    J = point_weights'*I; 
    Kxx = (mu_x.^2.*point_weights)'*I; 
    output = (C*sigma_a*(steady_temp^4/(4*pi) - I) + 3*mu_x*vx*steady_temp^4/(4*pi)*sigma_a + ...
        mu_x*vx*sigma_a.*I - sigma_a*vx^2/C*(J + Kxx));
end
% end function 
x0 = 0.1*ones(na,1);
dx = 0.000001;
[I_analytic] = secantm(x0,dx,@func2minimize);
%pack 
mu_z = 0.881917103688197; %into one buffer
out_count = 0;
for i=1:na
    if (abs(mu_remainder(i) - mu_z) < 1e-6)
        out_count = out_count +1; 
        x_out(out_count) = I_analytic(i).*mu_x(i);
        y_out(out_count) = I_analytic(i).*mu_y(i);
    end
end
plot(x_out, y_out,'--b'); 

mu_z = 0.3333333; %into one buffer
out_count = 0;
for i=1:na
    if (abs(mu_remainder(i) - mu_z) < 1e-6)
        out_count = out_count +1; 
        x_out2(out_count) = I_analytic(i).*mu_x(i);
        y_out2(out_count) = I_analytic(i).*mu_y(i);
    end
end
plot(x_out2, y_out2,'--k'); 

end
