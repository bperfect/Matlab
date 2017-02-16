function [Hz,z_w,z_r]=get_depth_Hz_ROMS(h,hc,s_rho,s_w,Cs_r,Cs_w,zeta,Vtransform);
% function get_depth_Hz_ROMS(h,hc,s_rho,s_w,Cs_r,Cs_w,zeta);
% This function computes the vertical grid metrics, z_r,z_w and Hz for a
% ROMS grid as a function of time:
% INPUT: 
% h    = Water depth (+ve in meters)
% hc   = Critical depth controlling the vertical stretchin (in m)
% s_rho= S coordinates at RHO points
% s_w  = S coordinates at w points
% Cs_r = S coordinate stretching curves at rho points
% Cs_w = S coordinate stretching curves at w points
% zeta = sea surface elevation (m)
% Vtransform= 1 or 2 

N1         = length(s_rho);
N2         = length(s_w);
[Mp,Lp]    = size(h);
j          = h==0;
h(j)       = eps;
z_w(1,:,:) = -h;

cff_r      = hc*s_rho;
cff_w      = hc*s_w;
cff1_r     = Cs_r; %originally Cs_r
cff1_w     = Cs_w; %originally Cs_w
hinv       = 1./(hc+h);
hwater     = h;
count=0;

for i=1:1:N1
    count = count+1;
    cff2_r= (cff_r(i)+cff1_r(i)*hwater).*hinv;
    cff2_w= (cff_w(i+1)+cff1_w(i+1)*hwater).*hinv;
    z_w(i+1,:,:) = zeta + (zeta+hwater).*cff2_w;
    z_r(i,:,:)   = zeta + (zeta+hwater).*cff2_r;
    Hz(i,:,:)    = z_w(i+1,:,:)-z_w(i,:,:);
        
end
