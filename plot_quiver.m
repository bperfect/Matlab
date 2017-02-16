clear all
fname='/home/bperfect/ROMS/startup/ocean_his_coarse.nc';
Z1    = 17;
V1    = ncread(fname,'v',[1 1 Z1 1],[Inf Inf 1 Inf],[1 1 1 1]);
U1    = ncread(fname,'u',[1 1 Z1 1],[Inf Inf 1 Inf],[1 1 1 1]);
V1    = squeeze(V1);
U1    = squeeze(U1);

x_rho=ncread(fname,'x_rho');
y_rho=ncread(fname,'y_rho');

V1(:,end+1,:)=V1(:,end,:);
U1(end+1,:,:)=U1(end,:,:);

V1 = double(V1);
U1 = double(U1);

L1=2;
L2=2;
fac=100000;
for i=1:1:2530
quiver(x_rho(1:L1:end,1:L2:end),y_rho(1:L1:end,1:L2:end),fac*squeeze(U1(1:L1:end,1:L2:end,i)),fac*squeeze(V1(1:L1:end,1:L2:end,i)),0,'k');
axis equal
xlim([0 150000]);
ylim([0 100000]);
pause(0.25); 
end