function [uT mom vol uT_in mom_in vol_in] = boundaryFlux(u,v,rho,T,x_rho,y_rho,z,t)
% u,v,T should be in rho coordinates, but vertical doesn't need to be
% corrected (as long as we're on the boundary, where topography is flat)
u=convert2rho(u,'u');
v=convert2rho(v,'v');
rho=convert2rho(rho,'rho');
T=convert2rho(T,'rho')+200;
x_rho=convert2rho(x_rho,'rho');
y_rho=convert2rho(y_rho,'rho');

x=x_rho(1:end,1);
y=squeeze(y_rho(1,1:end));


north = -v(:,end,:,:).*T(:,end,:,:).*rho(:,end,:,:);
south = v(:,1,:,:).*T(:,1,:,:).*rho(:,1,:,:);
east = -u(end,:,:,:).*T(end,:,:,:);
west = u(1,:,:,:).*T(1,:,:,:);
%gives net flux as a function of time
uT = squeeze(trapz(z,trapz(x,north+south)+trapz(y,east+west)));
uT_in = squeeze(trapz(z,trapz(x,south)+trapz(y,west)));

north = -v(:,end,:,:).*rho(:,end,:,:);
south = v(:,1,:,:).*rho(:,1,:,:);
east = -u(end,:,:,:).*rho(end,:,:,:);
west = u(1,:,:,:).*rho(1,:,:,:);
%gives net flux as a function of time
vol = squeeze(trapz(z,trapz(x,north+south)+trapz(y,east+west)));
vol_in = squeeze(trapz(z,trapz(x,south)+trapz(y,west)));

north = -v(:,end,:,:).^2.*rho(:,end,:,:);
south = v(:,1,:,:).^2.*rho(:,1,:,:);
east = -u(end,:,:,:).^2.*rho(end,:,:,:);
west = u(1,:,:,:).^2.*rho(1,:,:,:);
%gives net flux as a function of time
mom = squeeze(trapz(z,trapz(x,north+south)+trapz(y,east+west)));
mom_in = squeeze(trapz(z,trapz(x,south)+trapz(y,west)));

% figure
% hold on
% plot(t,mom./mom_in)
% plot(t,vol./vol_in)
% plot(t,uT./uT_in)
% title('Flux Balance')
% ylabel('Proportion of error')
% xlabel('Timestep')
% legend('Momentum','Volume','Temperature')

