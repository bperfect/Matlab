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


%% Obtain vertical coordinates in physical dimensions
[Hz,z_w,z_r]=get_depth_Hz_ROMS(h,hc,s_rho,s_w,Cs_r,Cs_w,zeta(:,:,1),2);

%% Choose time indices and physical vertical coordinates (meters) to interpolate
z_query = -100:-400:-2900;
t_query = 1:200:1700;

%% Convert to rho coordinates
rvort_test = convert2rho(rvort(:,:,:,t_query),'psi');
%pvort_test = convert2rho(pvort(:,:,:,t_query),'psi');
rho_test = convert2rho(rho(:,:,:,t_query),'rho');
z_test = convert2rho(permute(z_r,[2 3 1]),'rho');
dims = size(rvort_test);

for time = 1:dims(4)
    for row = 1:dims(1)
        for col = 1:dims(2)
            %spline interpolation, and for values outside the domain, use a given value (currently 0)
            %u_final(row,col,:,time) = interp1(z_r(:,row,col),squeeze(u_grid(row,col,:,time)),z_query','spline',0);
            %v_final(row,col,:,time) = interp1(z_r(:,row,col),squeeze(v_grid(row,col,:,time)),z_query','spline',0);
            rvort_final(row,col,:,time) = interp1(squeeze(z_test(row,col,:)),squeeze(rvort_test(row,col,:,time)),z_query,'spline',0);
            %pvort_final(row,col,:,time) = interp1(squeeze(z_test(row,col,:)),squeeze(pvort_test(row,col,:,time)),z_query,'spline',0);
            rho_final(row,col,:,time) = interp1(squeeze(z_test(row,col,:)),squeeze(rho_test(row,col,:,time)),z_query,'spline',0);
        end
    end
    time
end

dx = y_rho(1,2)-y_rho(1,1);
dy = x_rho(2,1)-x_rho(1,1);

figure
for t = 2:length(t_query)
    subplot(3,3,t)
    hold on
    for i=2:length(z_query)-1
        [psd,radspec,logspec,kx,ky,klin] = psd2d(rvort_final(:,:,i,t)',dx,0);
        plot( klin(1:(length(klin)-1)) , logspec ) ; xlabel ('k') ; ylabel ('energy'); title ('Radially averaged power spectrum')
    end
    title(['rvort time index = ' int2str(t_query(t))])
    legend(cellstr(num2str(z_query'))')
end

% figure
% for t = 1:length(t_query)
%     subplot(2,3,t)
%     hold on
%     for i=1:length(z_query)
%         [psd,radspec,logspec,kx,ky,klin] = psd2d(pvort_final(:,:,i,t),dx,0);
%         plot( klin(1:(length(klin)-1)) , logspec ) ; xlabel ('k') ; ylabel ('energy'); title ('Radially averaged power spectrum')
%     end
%     title(['pvort time index = ' int2str(t_query(t))])
%     legend(cellstr(num2str(z_query'))')
% end

figure
for t = 1:length(t_query)
    subplot(3,3,t)
    hold on
    for i=1:length(z_query)
        [psd,radspec,logspec,kx,ky,klin] = psd2d(rho_final(:,:,i,t),dx,0);
        plot( klin(1:(length(klin)-1)) , logspec ) ; xlabel ('k') ; ylabel ('energy'); title ('Radially averaged power spectrum')
    end
    title(['rho time index = ' int2str(t_query(t))])
    legend(cellstr(num2str(z_query'))')
end