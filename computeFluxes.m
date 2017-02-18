%% Obtain flow field in physical coordinates
reread=1;
cd '/home/bperfect/seamounts/datasets/run2';
cd '/home/bperfect/ROMS/startup';
hisname = 'ocean_his_coarse.nc';
avgname = 'ocean_avg.nc';
if reread==1;
di=1;    
%% field variables
u=ncread(hisname,'u',[1 1 1 1],[Inf Inf Inf Inf],[1 1 1 di]);
v=ncread(hisname,'v',[1 1 1 1],[Inf Inf Inf Inf],[1 1 1 di]);
rho=ncread(hisname,'rho',[1 1 1 1],[Inf Inf Inf Inf],[1 1 1 di]);
zeta = ncread(hisname,'zeta',[1 1 1],[Inf Inf Inf],[1 1 1]);
temp=ncread(hisname,'temp',[1 1 1 1],[Inf Inf Inf Inf],[1 1 1 di]);
%rvort=nc_read(avgname,'rvorticity');
%pvort=nc_read(avgname,'pvorticity');

%% vertical grid variables
h=ncread(hisname,'h');
s_w=ncread(hisname,'s_w');
s_rho=ncread(hisname,'s_rho');

%% y-coordinate grid variables
y_rho=ncread(hisname,'y_rho');
%%  x-coordinate grid variables
x_rho=ncread(hisname,'x_rho');
%% New variables
t=ncread(hisname,'ocean_time',[1],[Inf],[di]);
t=t/3600/24; %convert from seconds to days
hc=ncread(hisname,'hc');
Cs_w=ncread(hisname,'Cs_w');
Cs_r=ncread(hisname,'Cs_r');

end
%% Obtain vertical coordinates in physical dimensions
[Hz,z_w,z_r]=get_depth_Hz_ROMS(h,hc,s_rho,s_w,Cs_r,Cs_w,zeta(:,:,1),2);


%% boundary flux balance
[uT mom vol uT_in mom_in vol_in]=boundaryFlux(u,v,rho,temp,x_rho,y_rho,Cs_r,t);
%% energy balance
[E Et]=energyBalance(u,v,rho,x_rho,y_rho,z_r,t);


%% Re-load variables for the fine case
hisname = 'ocean_his_fine.nc';
if reread==1;
di=6;    
%% field variables
u=ncread(hisname,'u',[1 1 1 1],[Inf Inf Inf Inf],[1 1 1 di]);
v=ncread(hisname,'v',[1 1 1 1],[Inf Inf Inf Inf],[1 1 1 di]);
rho=ncread(hisname,'rho',[1 1 1 1],[Inf Inf Inf Inf],[1 1 1 di]);
zeta = ncread(hisname,'zeta',[1 1 1],[Inf Inf Inf],[1 1 1]);
temp=ncread(hisname,'temp',[1 1 1 1],[Inf Inf Inf Inf],[1 1 1 di]);
%rvort=nc_read(avgname,'rvorticity');
%pvort=nc_read(avgname,'pvorticity');

%% vertical grid variables
h=ncread(hisname,'h');
s_w=ncread(hisname,'s_w');
s_rho=ncread(hisname,'s_rho');

%% y-coordinate grid variables
y_rho=ncread(hisname,'y_rho');
%%  x-coordinate grid variables
x_rho=ncread(hisname,'x_rho');
%% New variables
t2=ncread(hisname,'ocean_time',[1],[Inf],[di]);
t2=t2/3600/24; %convert from seconds to days
hc=ncread(hisname,'hc');
Cs_w=ncread(hisname,'Cs_w');
Cs_r=ncread(hisname,'Cs_r');

end
%% Obtain vertical coordinates in physical dimensions
[Hz,z_w,z_r]=get_depth_Hz_ROMS(h,hc,s_rho,s_w,Cs_r,Cs_w,zeta(:,:,1),2);


%% boundary flux balance
[uT2 mom2 vol2 uT_in2 mom_in2 vol_in2]=boundaryFlux(u,v,rho,temp,x_rho,y_rho,Cs_r,t2);
%% energy balance
[E2 Et2]=energyBalance(u,v,rho,x_rho,y_rho,z_r,t);

%t2=t2-4200000+t(end); %correction for run1 timing issues
figure
hold on
plot(t,Et)
plot(t2,Et2)
xlabel('Time (days)')
ylabel('Kinetic Energy')
legend('Coarse Startup','Fine')
title('Horizontal Kinetic Energy')

figure
hold on
plot(t,mom./mom_in)
plot(t,vol./vol_in)
plot(t,uT./uT_in)
plot(t2,mom2./mom_in2)
plot(t2,vol2./vol_in2)
plot(t2,uT2./uT_in2)
title('Flux Balance')
ylabel('Proportion of error')
xlabel('Timestep')
legend('Coarse Momentum','Coarse Volume','Coarse Temperature','Fine Momentum','Fine Volume','Fine Temperature')


