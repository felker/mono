%Solve steady state integral equation, eq 17 in Jiang14 with dI/dt=0
%solve at one arbitrary spatial point, nx/2,ny/2
%np is number of discretization levels on [0,2pi]
function [] = steady_state()
sigma_a = 100;
steady_temp = 1;% 1.023395585899979;
np = 2000; %total number of rays 
mu_z = 0.3333333; %0.881917103688197; 
mu_remainder = sqrt(1-mu_z^2);
mu_x = mu_remainder*cos(linspace(0,pi,np/2));
mu_y = sqrt(mu_remainder^2-mu_x.^2); 
mu_x = [mu_x -mu_x]';
mu_y = [mu_y  -mu_y]';
%analytic solution is symmetric about x-axis
point_weights = ones(np,1)/np; 

C = 10;
vx = 3.0; %2.955596777221964; %0.3*C;
%error is with the coupled equation solver. this system is linear though!
%build A
A = zeros(np,np); 
for i=1:np %each row
    for j=1:np % each column
        if (j==i) %diagonal
            A(i,j) = sigma_a*(-C + mu_x(i)*vx - vx^2/C*point_weights(i)- vx^2/C*mu_x(i)^2*point_weights(i));
        else
            A(i,j) = sigma_a*(-vx^2/C*point_weights(j)- vx^2/C*mu_x(j)^2*point_weights(j));            
        end
    end
end
%build b
b = -1*sigma_a*steady_temp^4/(4*pi)*(C+3*vx*mu_x); 
I_analytic = A\b; 

% function output = func2minimize(I) %from R^np --> R^np
%     J = point_weights'*I; 
%     Kxx = (mu_x.^2.*point_weights)'*I; 
%     output = (C*sigma_a*(steady_temp^4/(4*pi) - I) + 3*mu_x*vx*steady_temp^4/(4*pi)*sigma_a + ...
%         mu_x*vx*sigma_a.*I - sigma_a*vx^2/C*(J + Kxx));
% end
% % end function 
% x0 = 0.1*ones(np,1);
% dx = 0.000001;
% [I_analytic] = secantm(x0,dx,@func2minimize);
x_out = mu_x.*I_analytic; 
y_out = mu_y.*I_analytic; 
out_count = 0;
% for i=1:np 
%         out_count = out_count +1; 
%         x_out(out_count) = I_analytic(i).*mu_x(i);
%         y_out(out_count) = I_analytic(i).*mu_y(i);
%     end
% end
plot(x_out, y_out,'--k'); 

end
