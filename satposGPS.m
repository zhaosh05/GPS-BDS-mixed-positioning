function [satp,tcorr] = satposGPS(t,eph)
% satellite position of GPS

GM = 3.986005e14;  % earth's universal gravitational parameter m^3/s^2
Omegae_dot = 7.2921151467e-5; % earth rotation rate, rad/s
F = -4.442807633e-10; % Constant, [sec/(meter)^(1/2)]

svprn   =   eph(1);
af2	  =   eph(2);
M0	     =   eph(3);
roota   =   eph(4);
deltan  =   eph(5);
ecc	  =   eph(6);
omega   =   eph(7);
cuc	  =   eph(8);
cus	  =   eph(9);
crc	  =  eph(10);
crs	  =  eph(11);
i0	     =  eph(12);
idot    =  eph(13);
cic	  =  eph(14);
cis	  =  eph(15);
Omega0  =  eph(16);
Omegadot=  eph(17);
tgd = eph(18);
af0	  =  eph(19);
af1	  =  eph(20);
toc	  =  eph(21);
toe = eph(21);

% Procedure for coordinate calculation
A = roota*roota;
tk = check_t(t-toc);
tcorr = (af2*tk + af1)*tk + af0;

n0 = sqrt(GM/A^3);
n = n0+deltan;
M = M0+n*tk;
M = rem(M+2*pi,2*pi);
E = M;
for i = 1:10
   E_old = E;
   sinE=sin(E);
   E = M+ecc*sinE;
   dE = rem(E-E_old,2*pi);
   if abs(dE) < 1.e-12
      break;
   end
end
E = rem(E+2*pi,2*pi);
cosE=cos(E);
dE = n/(1-ecc * cosE);
dtr=F * ecc * roota * sinE;
v = atan2(sqrt(1-ecc^2)*sinE, cosE-ecc);
phi = v+omega;
phi = rem(phi,2*pi);
dphi=sqrt(1-ecc^2)*dE/(1-ecc*cosE);
cos2phi=cos(2*phi);
sin2phi=sin(2*phi);
u = phi		         + cuc*cos2phi+cus*sin2phi;
du=(1+2*(cus * cos2phi-cuc * sin2phi))*dphi;
r = A*(1-ecc*cosE) + crc*cos2phi+crs*sin2phi;
dr= A * ecc *sin(E) *dE + ...
        2 *( crs * cos2phi - crc * sin2phi) ...
        * dphi;
i = i0+idot*tk	      + cic*cos2phi+cis*sin2phi;
di = 2 * (cis * cos2phi - cic * sin2phi) ...
        * dphi + idot;
Omega = Omega0+(Omegadot-Omegae_dot)*tk-Omegae_dot*toe;
dOmega = Omegadot - Omegae_dot;
Omega = rem(Omega+2*pi,2*pi);

cosu=cos(u);sinu=sin(u);
cosi=cos(i);sini=sin(i);
cosOmega=cos(Omega);sinOmega=sin(Omega);

xdash = cosu*r;
ydash = sinu*r;
satp(1,1) = xdash*cosOmega-ydash*cosi*sinOmega;
satp(2,1) = xdash*sinOmega+ydash*cosi*cosOmega;
satp(3,1) = ydash*sini;
dxdash = dr * cosu - r * sinu * du;
dydash = dr * sinu + r * cosu * du;

Vx = dxdash * cosOmega -  dydash * cosi * sinOmega ...
+ ydash * sinOmega * sini * di - (xdash * sinOmega + ...
ydash * cosi * cosOmega) * dOmega;

Vy = dxdash * sinOmega + dydash * cosi * cosOmega - ...
ydash * sini *cosOmega * di + (xdash * cosOmega - ...
ydash * cosi * sinOmega) * dOmega;

Vz = dydash * sini + ydash * cosi * di;
satp(4,1) = Vx;
satp(5,1) = Vy;
satp(6,1) = Vz;

tcorr=tcorr+dtr-tgd;
