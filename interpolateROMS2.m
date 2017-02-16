%% Obtain flow field in physical coordinates
reread=1;
cd '/home/bperfect/ROMS/sweep_fourier/f1e-5N1e-3/';
hisname = 'ocean_his.nc';
avgname = 'ocean_avg.nc';
if reread==1;
    
%field variables
%u=nc_read(hisname,'u');
%v=nc_read(hisname,'v');
%w=nc_read(filename,'w');
%rho=nc_read(filename,'rho');
zeta = nc_read(hisname,'zeta');
%temp=nc_read(filename, 'temp');
rvort=nc_read(avgname,'rvorticity');

%vertical grid variables
h=nc_read(hisname,'h');
s_w=nc_read(hisname,'s_w');
s_rho=nc_read(hisname,'s_rho');

%y-coordinate grid variables
y_rho=nc_read(hisname,'y_rho');
y_u=nc_read(hisname,'y_u');
y_v=nc_read(hisname,'y_v');
y_psi=nc_read(hisname,'y_psi');
%x-coordinate grid variables
x_rho=nc_read(hisname,'x_rho');
x_v=nc_read(hisname,'x_v');
x_psi=nc_read(hisname,'x_psi');
x_u=nc_read(hisname,'x_u');
%%% New variables

hc=nc_read(hisname,'hc');
Cs_w=nc_read(hisname,'Cs_w');
Cs_r=nc_read(hisname,'Cs_r');
end
[Hz,z_w,z_r]=get_depth_Hz_ROMS(h,hc,s_rho,s_w,Cs_r,Cs_w,zeta(:,:,1),2);

x_grid = x_rho(2:end-1,2:end-1);
y_grid = y_rho(2:end-1,2:end-1);
z_grid = reshape(z_r(:,2:end-1,2:end-1),size(x_grid))
%% x_grid,y_grid,z_grid specify interior rho coordinates

% a=size(u);
% u_grid = zeros(a(1)-1,a(2)-2,a(3),a(4));
% a=size(v);
% v_grid = zeros(a(1)-2,a(2)-1,a(3),a(4));
% w_grid = w(2:end-1,2:end-1,:,:);

%fix u-point variables onto rho coords
% use linear interpolation

% for i = 1:length(x_grid(:,1))
%     u_grid(i,:,:,:) = 1/2*(u(i+1,2:end-1,:,:)+u(i,2:end-1,:,:));
% end
% %fix v-point variables onto rho coords
% for i = 1:length(y_grid(1,:))
%     v_grid(:,i,:,:) = 1/2*(v(2:end-1,i+1,:,:)+v(2:end-1,i,:,:));
% end

z_query = -100:-100:-3000
time=3333;
for row =1:91
    for col = 1:61
        u_final(row,col,:,time) = interp1(z_r(:,row,col),squeeze(rvort(row,col,:,time)),z_query','spline',0);

% for time = length(u_grid(1,1,1,:))
%     for row = 1:length(u_grid(:,1,1,1))
%         for col = 1:length(u_grid(1,:,1,1))
%             %spline interpolation, and for values outside the domain,
%             %use a given value (currently 0)
%             u_final(row,col,:,time) = interp1(z_r(:,row,col),squeeze(u_grid(row,col,:,time)),z_query','spline',0);
%             v_final(row,col,:,time) = interp1(z_r(:,row,col),squeeze(v_grid(row,col,:,time)),z_query','spline',0);
%         end
%     end
%     time
% end

        
            

u_in = scatteredInterpolant(x_u(:),y_u(:),h_u_3d(:),u_time(:));
u_rho = u_in(rho_grid_x,rho_grid_y,rho_grid_z);

pause()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

