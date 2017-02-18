clear all
fname='/home/bperfect/seamounts/datasets/run2/ocean_his_fine.nc';
Z1    = 3;
V1    = ncread(fname,'v',[1 1 Z1 311],[Inf Inf Inf 1],[1 1 1 1]);
U1    = ncread(fname,'u',[1 1 Z1 311],[Inf Inf Inf 1],[1 1 1 1]);
V1    = squeeze(V1);
U1    = squeeze(U1);

h=ncread(fname,'h');
s_w=ncread(fname,'s_w');
s_rho=ncread(fname,'s_rho');
y_rho=ncread(fname,'y_rho');
x_rho=ncread(fname,'x_rho');
%% New variables
t=ncread(fname,'ocean_time');
t=t/3600/24; %convert from seconds to days
hc=ncread(fname,'hc');
Cs_w=ncread(fname,'Cs_w');
Cs_r=ncread(fname,'Cs_r');
x_rho=ncread(fname,'x_rho');
y_rho=ncread(fname,'y_rho');
zeta = ncread(fname,'zeta',[1 1 1],[Inf Inf Inf],[1 1 1]);

%% Obtain vertical coordinates in physical dimensions
[Hz,z_w,z_r]=get_depth_Hz_ROMS(h,hc,s_rho,s_w,Cs_r,Cs_w,zeta(:,:,1),2);

%% Choose time indices and physical vertical coordinates (meters) to interpolate
z_query = -50:-100:-2950;

z_r = convert2rho(z_r,'rho');
u=convert2rho(U1,'u');
v=convert2rho(V1,'v');
dims = size(u);

for row = 1:dims(1)
    for col = 1:dims(2)
        u_final(row,col,:) = interp1(z_r(:,row,col),squeeze(u(row,col,:)),z_query','spline',0);
        v_final(row,col,:) = interp1(z_r(:,row,col),squeeze(v(row,col,:)),z_query','spline',0);
    end
end

L1=4;
L2=4;
fac=20000;
for i=1:length(z_query)
quiver(x_rho(1:L1:end,1:L2:end),y_rho(1:L1:end,1:L2:end),fac*squeeze(U1(1:L1:end,1:L2:end,i)),fac*squeeze(V1(1:L1:end,1:L2:end,i)),0,'k');
axis equal
axis off
title(['Z= ' int2str(z_query(i)) 'meters'])
name = ['quiver' int2str(i+100) '.png'];
saveas(gcf,name)
pause(0.25); 
end