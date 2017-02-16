function [E Et] = energyBalance(u,v,rho,x_rho,y_rho,Hz,t)
% u,v,T should be in rho coordinates, but vertical doesn't need to be
% corrected (as long as we're on the boundary, where topography is flat)
Hz = convert2rho(shiftdim(Hz,1),'rho');%problem is here
u=convert2rho(u,'u');
v=convert2rho(v,'v');
rho=convert2rho(rho,'rho');
x_rho=convert2rho(x_rho,'rho');
y_rho=convert2rho(y_rho,'rho');

x=x_rho(1:end,1);
y=squeeze(y_rho(1,1:end));

u2=rho.*u.^2/2;
v2=rho.*v.^2/2;

E=zeros(size(squeeze(u(:,:,1,:))));

volume = Hz(end,1,1)*x(end)*y(end);
E_baro = -volume*mean(mean(mean(mean(rho))))*(0.1*tanh(t/10)).^2/2;

for i=1:length(x)
    for j=1:length(y)
        E(i,j,:) = trapz(squeeze(Hz(i,j,:)),u2(i,j,:,:))+trapz(squeeze(Hz(i,j,:)),v2(i,j,:,:));
    end
end

Et = squeeze(mean(diff(y))*mean(diff(x))*trapz(trapz(E)));

figure
plot(t,E_baro)
hold on
plot(t,Et)
xlabel('Time (days)')
ylabel('Kinetic Energy')
legend('Barotropic Prediction','Actual')



