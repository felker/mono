function [y_out] = time_series_output(y_out,time_out,rad_energy,rad_flux,rad_pressure,v,nx,ny,C,GasMomentum,GasKE,GasTE,photon_momentum)
%TIME_SERIES_OUTPUT Uncomment plots which are supposed to print out line
%plot values every interval timesteps, including the first one

%Need to constantly update function definition for what variables are needed depending on desired
%output
persistent out_count;

if isempty(out_count)
    out_count = 0;
end
out_count = out_count+1;

%Set to 2 if you want to ignore the initial condition step, t=0
start_index = 2; 

%Jiang14 Figure 2
y_out(out_count,1) = rad_flux(nx/2,ny/2,1);
%advective flux
y_out(out_count,2) = v(nx/2,ny/2,1)*(rad_energy(nx/2,ny/2) + rad_pressure(nx/2,ny/2,1,1))/C; 

y_out(out_count,3) = GasMomentum(nx/2,ny/2,1)+photon_momentum(nx/2,ny/2,1);
y_out(out_count,4) = GasMomentum(nx/2,ny/2,1);
y_out(out_count,5) = GasTE(nx/2,ny/2)+GasKE(nx/2,ny/2)+rad_energy(nx/2,ny/2);
y_out(out_count,6) = GasTE(nx/2,ny/2)+GasKE(nx/2,ny/2);

hi = subplot(3,1,1);
plot(time_out(start_index:out_count),y_out(start_index:out_count,1)','-',time_out(start_index:out_count),y_out(start_index:out_count,2)','-');
axis([0 0.1 0.25 0.5])
hi = subplot(3,1,2);
plot(time_out(start_index:out_count),y_out(start_index:out_count,3)','--',time_out(start_index:out_count),y_out(start_index:out_count,4)','-');
axis([0 0.1 2.9 3.15])
hi = subplot(3,1,3);
plot(time_out(start_index:out_count),y_out(start_index:out_count,5)','--',time_out(start_index:out_count),y_out(start_index:out_count,6)','-');


%pause(1.0);


     %Figure 1
%     y_out(out_count) = rad_energy(nx/2,ny/2);
%     y_out2(out_count) = temp(nx/2,ny/2)^4;
% % %     figure(1);
%       legend('E_r','T^4');

%      xlabel('time (s)');
%      ylabel('E_r,T^4');
%      title('\sigma_a = 1');


end

