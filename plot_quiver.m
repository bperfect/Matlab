clear all
fname='/home/bperfect/seamounts/datasets/run2/ocean_his_fine.nc';
Z1    = 3;
V1    = ncread(fname,'v',[1 1 Z1 311],[Inf Inf Inf 1],[1 1 1 1]);
U1    = ncread(fname,'u',[1 1 Z1 311],[Inf Inf Inf 1],[1 1 1 1]);
V1    = squeeze(V1);
U1    = squeeze(U1);

x_rho=ncread(fname,'x_rho');
y_rho=ncread(fname,'y_rho');

V1(:,end+1,:)=V1(:,end,:);
U1(end+1,:,:)=U1(end,:,:);

V1 = double(V1);
U1 = double(U1);

L1=4;
L2=4;
fac=20000;
for i=1:1:35
%subplot(2,4,i)
quiver(x_rho(1:L1:end,1:L2:end),y_rho(1:L1:end,1:L2:end),fac*squeeze(U1(1:L1:end,1:L2:end,i)),fac*squeeze(V1(1:L1:end,1:L2:end,i)),0,'k');
axis equal
axis off
% xlim([0 150000]);
% ylim([0 100000]);
pause(0.25); 
end