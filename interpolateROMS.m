%% Obtain flow field in physical coordinates
reread=1;
cd '/home/bperfect/ROMS/sweep_his_files';
filename = 'sweep_f_1e-5_N_1e-3_hyak.nc';
if reread==1;
%field variables
u=nc_read(filename,'u');
v=nc_read(filename,'v');
w=nc_read(filename,'w');
rho=nc_read(filename,'rho');
zeta = nc_read(filename,'zeta');
temp=nc_read(filename, 'temp');

%vertical grid variables
h=nc_read(filename,'h');
s_w=nc_read(filename,'s_w');
s_rho=nc_read(filename,'s_rho');

%y-coordinate grid variables
y_rho=nc_read(filename,'y_rho');
y_u=nc_read(filename,'y_u');
y_v=nc_read(filename,'y_v');
y_psi=nc_read(filename,'y_psi');
%x-coordinate grid variables
x_rho=nc_read(filename,'x_rho');
x_v=nc_read(filename,'x_v');
x_psi=nc_read(filename,'x_psi');
x_u=nc_read(filename,'x_u');

hc=nc_read(filename,'hc');
Cs_w=nc_read(filename,'Cs_w');
Cs_r=nc_read(filename,'Cs_r');
[Hz,z_w,z_r]=get_depth_Hz_ROMS(h,hc,s_rho,s_w,Cs_r,Cs_w,zeta(:,:,1),2);
pause()

%% create the master grid that everything will be expressed in terms of
hLevels = length(s_rho);
rho_grid_z = zeros(length(x_rho(:,1)),length(x_rho(1,:)),length(s_rho));
rho_grid_x = repmat(x_rho,1,1,hLevels);
rho_grid_y = repmat(y_rho,1,1,hLevels);

%% Generate our z-dimension variables
% Express height function as a 3D tensor in u and v coordinates for later interpolation
h_in = griddedInterpolant(x_rho,y_rho,h); %interpolate h (given in rho coords)
h_u_2d = h_in(x_u(:,:,1),y_u(:,:,1));
h_v_2d = h_in(x_v(:,:,1),y_v(:,:,1));

h_u_3d = zeros(length(h_u_2d(:,1)),length(h_u_2d(1,:)),hLevels);
h_v_3d = zeros(length(h_v_2d(:,1)),length(h_v_2d(1,:)),hLevels);
h_w_3d = zeros(length(h(:,1)),length(h(1,:)),length(s_w));

%give our z-coordinate grids vertical variation
for i=1:hLevels
    h_u_3d(:,:,i) = h_u_2d*s_rho(i);
    h_v_3d(:,:,i) = h_v_2d*s_rho(i);
    rho_grid_z(:,:,i) = h*s_rho(i);
end

for i = 1:length(s_w)
    h_w_3d(:,:,i) = h*s_w(i);
end

%% Generate our x and y-dimension variables
x_u = repmat(x_u,1,1,hLevels);
y_u = repmat(y_u,1,1,hLevels);
x_v = repmat(x_v,1,1,hLevels);
y_v = repmat(y_v,1,1,hLevels);
x_w = repmat(x_rho,1,1,length(s_w));
y_w = repmat(y_rho,1,1,length(s_w));

end %if reread=1
%% Grab our flow field at a given time
frames=248;
F1(frames) = struct('cdata',[],'colormap',[]);
fig=figure
for timeIndex = 1:248
%timeIndex = length(u(1,1,1,:));
u_time = u(:,:,:,timeIndex);
v_time = v(:,:,:,timeIndex);
w_time = w(:,:,:,timeIndex);

%% Generate interpolation functions and interpolate
%u interpolation
u_in = scatteredInterpolant(x_u(:),y_u(:),h_u_3d(:),u_time(:));
u_rho = u_in(rho_grid_x,rho_grid_y,rho_grid_z);
%v interpolation
v_in = scatteredInterpolant(x_v(:),y_v(:),h_v_3d(:),v_time(:));
v_rho = v_in(rho_grid_x,rho_grid_y,rho_grid_z);
%w interpolation
w_in = scatteredInterpolant(x_w(:),y_w(:),h_w_3d(:),w_time(:));
w_rho = w_in(rho_grid_x,rho_grid_y,rho_grid_z);

% %Confirm that everything is in the same dimensions
% size(u_rho)
% size(v_rho)
% size(w_rho)
% size(rho)

rho_t = rho(:,:,:,timeIndex);
rho_0 = rho(:,:,:,1);
rho_perturb = rho_t-rho_0;

% surf(squeeze(rho_grid_x(:,50,:)), squeeze(rho_grid_z(:,50,:)), squeeze(rho_perturb(:,50,:)))
% view(0,90)

surf(squeeze(rho_grid_x(:,71,:)), squeeze(rho_grid_z(:,71,:)), squeeze(u_rho(:,71,:)),'edgecolor','none')
view(0,90)
caxis([.08 .12])
F1(i)=getframe(fig);
end

vv = VideoWriter('newfile.avi');
open(vv)
writeVideo(vv,F1)
close(vv)

%% Richardson number plots
% ri = zeros(size(u_rho));
% g=9.8;rho_o = 1026;g_prime = g*.04/rho_o;
% for i=1:hLevels-1
%     ri(:,:,i) = -g/rho_o*(rho_t(:,:,i)-rho_t(:,:,i+1)).*...
%         (rho_grid_z(:,:,i)-rho_grid_z(:,:,i+1))./...
%         (sqrt(v_rho(:,:,i).^2+u_rho(:,:,i).^2)-sqrt(v_rho(:,:,i+1).^2+u_rho(:,:,i+1).^2)).^2;
% end
% ri(ri>5) = 5;
% ri(ri<-1) = -1;
% for i=1:hLevels-1
% imagesc(ri(:,:,i))
% colorbar
% caxis([-1 5])
% pause()
% end

%Notes: matlab isosurface requires a 3d grid of values 
