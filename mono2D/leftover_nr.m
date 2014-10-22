%left over multidimensional (2x2) NR solve in absorption term updates
%none of this worked correctly


           for m =1:MAX_ITER
       %         Build RHS vector representing the original nonlinear
       %        function, restricted to only the final two components
                F = zeros(2,1); 
                J = zeros(2,2);
                temp_step = temp(k,l);
                IN_step = intensity(k,l,na);
                
                                J(1,1) = A(na) + rho_a(k,l)*dt*vC(k,l)*(((A(na)*ones(na-1,1))./A(1:na-1))'*pw(1:na-1) + ...
                    pw(na)) + rho_a(k,l)*dt/C*(((vvnn(1:na-1).*pw(1:na-1))./A(1:na-1))'*A(na)*ones(na-1,1) + ...
                    pw(na)*vvnn(na));
                J(1,2) = -4*dt*rho_a(k,l)*C*(1+nv(k,l,na)/C)*temp_step^3/pi; 
                J(2,1) = -dt*(adiabatic -1)/(R*density(k,l))*(-8*pi*P*rho_a(k,l)*(((A(na)*squeeze(nv(k,l,1:na-1))./A(1:na-1))'*pw(1:na-1)...
                    + nv(k,l,na)*pw(na)) - c*vC(k,l)*(A(na)*ones(1,na-1)*((pw(1:na-1)./A(1:na-1))) + vvnn(na)*pw(na))) + ...
                    P*C*(1- vC(k,l))*rho_a(k,l)*4*pi*(pw(1:na-1)'*(A(na)*ones(na-1,1)./A(1:na-1)) + pw(na)));
                J(2,2) = 1 +4*(adiabatic -1)/(R*density(k,l))*dt*P*C*(1-vC(k,l))*rho_a(k,l)*temp_step^3; 
                
                
                                F(1) = IN_step - I_old(na) - dt*(C*rho_a(k,l)*(temp_step^4/(4*pi) - IN_step) + ...
                nv(k,l,na)*rho_a(k,l)*(3*temp_step^4/(4*pi) + IN_step) - rho_a(k,l)*C*vC(k,l)*(pw(1:na-1)'*D + ...
                pw(na)*IN_step) - rho_a(k,l)/C*((pw(1:na-1).*vvnn(1:na-1))'*D + pw(na)*vvnn(na)*IN_step));
            
                F(2) = temp_step - temp_old - (adiabatic -1)*dt/(R*density(k,l))*(-8*pi*P*rho_a(k,l)*( ...
                    (squeeze(nv(k,l,1:na-1)).*pw(1:na-1))'*D + nv(k,l,na)*pw(na)*IN_step - C*vC(k,l)*(D'*pw(1:na-1) + pw(na)*IN_step) - ...
                    1/C*(((vvnn(1:na-1).*pw(1:na-1))'*D + vvnn(na)*pw(na)*IN_step))) - P*C*(1-vC(k,l))*rho_a(k,l)*( ...
                    temp_step^4 - 4*pi*(pw(1:na-1)'*D + pw(na)*IN_step)));
                
                            %this item is from subtracting na-th line from lth line

                            D = ((A(na)*IN_step -I_old(na))*ones(na-1,1) + I_old(1:na-1))./A(1:na-1); 
          % end %done with NR solve at cell
          %  m;
                      %STATIC CASE: UPDATE OTHER RAYS
            %intensity(k,l,1:na-1) =  (I_old(r) -I_old(na))/(1 +dt*rho_a(k,l)*C) + IN_step; 
            %1D static case, update everything
%             intensity(k,l,1:na-1) =  (I_old(r) + C*rho_a(k,l)*temp_step^4/(4*pi))/(1 +dt*rho_a(k,l)*C);
%             temp(k,l) = temp_step;

  %ATTEMPT STATIC MIXED FRAME SOLUTION V=0 only
%                    D = zeros(na,1);
%                  for r=1:na 
%                      D(r) = (I_old(r) -I_old(na))/(1 +dt*rho_a(k,l)*C) + IN_step; 
%                  end;
%                 F(1) = (1 +dt*rho_a(k,l)*C)*IN_step - I_old(na) - C*rho_a(k,l)*temp_step^4/(4*pi);
%                 F(2) = temp_step- temp_old + dt*rho_a(k,l)*P*C/(density(k,l)*R)*(temp_step^4 - ...
%                     4*pi*D'*pw);
%                 J(1,1) = (1 +dt*rho_a(k,l)*C); 
%                 J(1,2) = -C*rho_a(k,l)/pi*temp_step^3; 
%                 J(2,1) = 1 + dt*rho_a(k,l)*P*C/(density(k,l)*R)*4*temp_step^3;
%                 J(2,2) = -4*pi*dt*rho_a(k,l)*P*C/(density(k,l)*R);
% 
%                 %Solve 2x2 system for I^n+1_na and T^n+1 steps
%                 delta_step = J\-F;
%                 IN_step = IN_step + delta_step(1);
%                 temp_step = temp_step + delta_step(2);
                %check convergence criteria. Maybe change this to delta
                %intensity?
%                 if max(abs(delta_step)) < delta_c
%                % if max(abs(delta_step./[IN_step temp_step]')) || (delta_step(1) == 0 && delta_step(2) == 0 && IN_step == 0 && temp_step ==0) < delta_c
%                     break;
%                 end  

             intensity(k,l,1:na-1) = ((A(na)*IN_step -I_old(na))*ones(na-1,1) + I_old(1:na-1))./A(1:na-1); 

                    %intensity(k,l,n) = (1-dt*(1+nv(k,l,n)/c)*rho_a(k,l))^-1*((IN_step - I_old(na))/dt + ...
                    %    rho_a(k,l)*IN_step + 3*temp_step^4/(4*pi*c)*rho_a(k,l)*(nv(k,l,n) - nv(k,l,na)) ...
                    %    - nv(k,l,na)/c*rho_a(k,l)*IN_step);           
                %end