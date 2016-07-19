%!  Calculate the depth of the bottom.
%  do i=is,ie ; do j=js,je
%    x=(G%geoLonT(i,j)-G%west_lon)/G%len_lon
%    y=(G%geoLatT(i,j)-G%south_lat)/G%len_lat
%!  This sets topography that has a reentrant channel to the south.
%    D(i,j) = -D0 * ( y*(1.0 + 0.6*cos(4.0*PI*x)) &
%                   + 0.75*exp(-6.0*y) &
%                   + 0.05*cos(10.0*PI*x) - 0.7 )
%    if (D(i,j) > max_depth) D(i,j) = max_depth
%    if (D(i,j) < min_depth) D(i,j) = 0.
%  enddo ; enddo
clc
clear all;

dx=1/45;dy=1/30;
xx=[dx/2:dx:1]; yy=[dy/2:dy:1]; D0=1 ; max_depth = 1.1; min_depth = 0;

[y,x]=meshgrid(yy,xx);

D = D0;
D = D-(1.2*spike(x,0.2)+1.2*spike(x-1,0.2)).*spike(min(0,y-.3),0.2)...% South America
      -1.2*spike(x-0.5,0.2).*spike(min(0,y-.55),0.2)...%Africa
      -1.1*spike(y-1,0.12)-1.1*spike(y,0.12)...% The great northern wall and Antarctica
      -1.2*(spike(x,0.12)+spike(x-1,0.12)).*spike(max(0,y-.06),0.12)...% Antarctic Peninsula
      -0.1*(cosbell(x,0.1)+cosbell(x-1,0.1))... % Drake Passage ridge 
      -0.5*cosbell(x-0.16,0.05).*cosbell(y-0.18,0.13).^0.4...%Scotia Arc East
      -0.4*cosbell(x-0.09,0.08).^0.4.*cosbell(y-0.26,0.05)...%Scotia Arc North
      -0.4*cosbell(x-0.08,0.08).^0.4.*cosbell(y-0.1,0.05)...%Scotia Arc South
      -0.05*cos(14*pi*x).*sin(14*pi*y)... %roughness
      -0.05*cos(20*pi*x).*cos(20*pi*y);   %roughness

% I we want to make sure there is a wall at South America at coarse res,
% we nned something like this:
jmin=find(yy>0.3,1,'first');
D(1,jmin:end)=0;  

D(D>max_depth)=max_depth;
D(D<min_depth)=min_depth;

% Zonal wind stress:
tau_max=0.2;
taux=tau_max*(cosbell(y-0.28,0.28).*(y<=0.28)...
     +(0.35+0.65*cos(pi*(y-0.28)/0.52)).*(y>0.28).*(y<=0.8)...
     -0.3*cosbell(y-0.8,0.2).*(y>0.8)...
     );
 
len_lon=90 ;len_lat= 50;west_lon=0; south_lat=-65;

lat=y*len_lat+south_lat;
lon=x*len_lon+west_lon;

f=4*pi/86400*sind(lat);

figure(1)
clf
[h,c]=contourf(lon,lat,-D,[-1.15:0.05:0.5]); hold on;
set(c,'edgecolor','none')
contour(lon,lat,f./D,[-3e-4:2e-5:0],'k')
colorbar

figure(2)
clf
plot(lat(1,:),taux(1,:));hold on;
plot(0.5*(lat(1,1:end-1)+lat(1,2:end)),diff(taux(1,:))*30);
