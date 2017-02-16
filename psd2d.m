 function [psd2d,radspec,logspec,kx,ky,klin] = psd2d(a,delta,verbosity)
%
% 
% [psd2d,amp,radspec,kx,ky]=psd2d(a,delta)
% Compute the Power spectral density for a 2d signal (psd2d) and the
% radially averaged power spectrum (radspec, logspec). The spectrum is
% averaged over a semicircular region on  the [kx+, ky]  plane.
% the matrix is extended automatically with a symmetric extension 
% (20% extension)
% (C) Fulvio Porzio, 2016, fulvio301@gmail.com 
% -------Input-------
% a     = Matrix
% delta = sampling interval in metres
% -------Output-------
% psd2d   = 2d power spectral density grid
% radspec = Radially averaged power spectrum 
% logspec = Logarithmic radially averaged spectrum
% kj      = angular wavenumber in radians/metric unit defined by user
    
    % -------Sets the matrix for fft-------
    [matr, window] = hanning2d (a, 30);
    percentage     = tie((length(a))*1.2);
    [matr, t]      = bordowav (matr,'sp0',percentage);
    [m, n]         = size (matr);
    
    
    % -------Sets the wavenumber vectors-------
    kx1 = mod ( 1/2+ ( 0:(m-1) )/m,1 )-1/2;   
    kx  = fftshift ( kx1.* ( 2*pi/delta ) );

    ky1 = mod ( 1/2+ ( 0:(n-1) )/n,1 )-1/2; 
    ky  = fftshift ( ky1.* ( 2*pi/delta ) );
    
    % -------Compute Power Spectrum-------
    amp   = fftshift ( abs (fft2(matr) ) );
    psd2d = (1./ ( 2.* (delta.*m) .* (delta.*n) ) ).* ( ( abs(amp) ).^2);
    
    % -------Compute Radial Spectrum-------
    [~, jx]     = find (kx>=0); 
    jy          = 1 : length(ky);
    kxpos       = kx (jx); 
    kypos       = ky (jy);
    [p,q]       = meshgrid (kxpos, kypos);
    kgrd        = sqrt ( (p.^2) + (q.^2) );
    zed         = zeros ( size(a, 1),size(a, 2) );
    
    for jk = 1:length(jx); jj = 1:length(jy);
            
        zed(jy(jj), jx(jk))  = 1 ;
        
    end
    
    psd2dposzed = psd2d.*zed ;
    
    psd2dposzed ( : , all(~psd2dposzed,1) ) = [];
    psd2dposzed ( all(~psd2dposzed,2) , : ) = [];
    

    
    klin      =  kxpos ;
    dk        =  (klin(2)/2);
    klin      =  klin - dk; 
    radspec   =  ones (length(klin),1 ) ;
    
    for i = 1 : (length(klin) - 1)
    
    rng         =  kgrd >= klin (i) & kgrd < klin (i+1);
    [ax, ay]    =  find (rng);
    rng         =  rng.*psd2dposzed;
    aa          =  find (rng);
    
        for jj = 1:length(ax);
        
        aa(jj) = psd2dposzed ( ax(jj), ay(jj) );
        
        end
    
    radspec (i) = nanmean (aa);
    
    end

radspec ( length(klin) ) = psd2dposzed ( size(psd2dposzed, 1), size(psd2dposzed, 2));    
    
      
    % -------plot results-------
    klin  = klin + dk;
    xplot = 0:delta:(size(a,1)*delta)-delta;
    yplot = 0:delta:(size(a,2)*delta)-delta;
    
    if max(xplot) > 10^3;
        
        xplot = xplot./1000;
        yplot = yplot./1000;
        m     = 'km';
        
    elseif max(xplot) < 10^3;
        m     = 'm';
    end
    
    signal  = 'Relative Vorticity';%input ('kind of signal?','s');
    measure = 'm';%input ('unit measure?','s');
    logspec = log (radspec(1:(length(klin)-1)));
    matr    = matr(t(1):t(2), t(3):t(4));
    if verbosity ==1
    figure; 
    size(xplot)
    size(yplot)
    size(a)
    subplot(2,3,1);   pcolor (xplot, yplot, a')      ;  xlabel(m);      ylabel(m);      shading flat;            title(signal);                       colormap jet; axis vis3d; axis tight
    subplot(2,3,4);   mesh   (xplot, yplot, a')      ;  xlabel(m);      ylabel(m);      zlabel(measure);         title(signal);                       colormap jet; axis vis3d; axis tight
    subplot(2,3,3);   pcolor (kx,ky,log (psd2d) )   ;  xlabel('Kx');   ylabel('Ky');   shading flat;            title('Power Spectral Density');     colormap jet; axis vis3d; axis tight
    subplot(2,3,6);   mesh   (kx,ky,log (psd2d) )   ;  xlabel('Kx');   ylabel('Ky');   zlabel('Energy');        title('Power Spectral Density');     colormap jet; axis vis3d; axis tight
    subplot(2,3,2);   mesh   (xplot, yplot, window') ;  xlabel(m);      ylabel(m);                               title('Hanning Window');             colormap jet; axis vis3d; axis tight
    subplot(2,3,5);   mesh   (xplot, yplot, matr')   ;  xlabel(m);      ylabel(m);      zlabel(measure);         title('Windowed signal');            colormap jet; axis vis3d; axis tight

    
    figure;
 
    plot( klin(1:(length(klin)-1)) , logspec, 'o') ; grid on; hold on; plot( klin(1:(length(klin)-1)) , logspec ) ; xlabel ('k') ; ylabel ('energy'); title ('Radially averaged power spectrum'); 
    end
 end
 
function [hh, hxy] = hanning2d(a, per)
%	[hh,hxy] = hanning2d(a,per);
%   Taper to zero the edges of a map according to Hanning window function
%   hxy = hanning2d(a,per);
%   per = 1; su tutta l'area (hanning forte)
%   per debole(100%) tapering m=n/10;
%
[n1, n2] = size(a);
m        = (n1/100) * per ;
m        = fix(m);
if n1 == m
    hx  = hanning(n1);
    hy  = hanning(n2)';
    hxy = hx*hy;
else
% hanning debole
    h                       = hanning(m);
    finestr                 = ones(n1, 1);
    finestr(1:m/2,1)        = h(1:m/2, 1);
    finestr(n1-m/2+1:n1,1)  = h(m/2+1:m, 1);
    hx                      = finestr;
    finestr                 = ones(n2, 1);
    finestr(1:m/2,1)        = h(1:m/2, 1);
    finestr(n2-m/2+1:n2,1)  = h(m/2+1:m, 1);
    hy                      = finestr';
    hxy                     = hx*hy;
end
hh = a.*hxy;
end
 
 function [val] = tie(num) 
 
    val = 2.*round (num./2) ;
 
 end
 
 function [grid4, t] = bordowav(a, mode, twopower)
% ------Function-----
%
% [grid4,t]=bordowav(a,mode,twopower);
% -------Output------
%
%   grid4    = extended matrix
%   t        = de-extending indexes
%
% -------Input-------
%
%   a        = matrix to extend
%   twopower = power of two or other dimension
%   mode     =  'zpd'    Zero extension'
%               'sp0'    Smooth extension of order 0
%               'spd'    Smooth extension of order 1
%               'sym'    Symmetric extension
%               'ppd'    Periodized extension (1)
%               'per'    Periodized extension (2)
%
if  mod(twopower-size(a,1),2) ~= 0;
    s1  = (twopower-size(a,1)+1)/2;
    t1  = s1+1;t2=t1+size(a,1)-1;
else s1 =(twopower-size(a,1))/2;
    
    t1=s1+1;t2=t1+size(a,1)-1;
    
end

if  mod(twopower-size(a,2),2)~= 0;
    s2  = (twopower-size(a,2)+1)/2;
    t3  = s2+1;t4=t3+size(a,2)-1;
    
else s2 = (twopower-size(a,2))/2;
     t3 = s2+1;t4=t3+size(a,2)-1;
end

grid4 = wextend(2,mode,a,[s1,s2]);
grid4 = grid4(1:twopower,1:twopower);
t     = [t1,t2,t3,t4];
 end