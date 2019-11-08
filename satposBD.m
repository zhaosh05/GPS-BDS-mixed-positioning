function [satp,tcorr] = satposBD(t,eph)
%Calculation of BDS satellites' X,Y,Z coordinates at time t
% reference: BeiDou Navigation Satellite System Signal In Space Interface Control Document-Open Service Signal B1I (Version 1.0) 

GM = 3.986004418e14;  % earth's universal gravitational parameter m^3/s^2
Omegae_dot = 7.2921150e-5; % earth rotation rate, rad/s

F          = -4.442807633e-10;

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
tgd	  =  eph(18);
af0	  =  eph(19);
af1	  =  eph(20);
toe	  =  eph(21);
toc=toe;


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
dtr = F * ecc * roota * sinE;
v = atan2(sqrt(1-ecc^2)*sinE, cosE-ecc);
phi = v+omega;
phi = rem(phi,2*pi);
dphi=sqrt(1-ecc^2)*dE/(1-ecc*cosE);
cos2phi=cos(2*phi);
sin2phi=sin(2*phi);
u = phi		         + cuc*cos2phi+cus*sin2phi;
du=(1+2*(cus * cos2phi-cuc * sin2phi))*dphi;
r = A*(1-ecc*cosE) + crc*cos2phi+crs*sin2phi;
dr= A * ecc *sinE *dE + ...
        2 * (crs * cos2phi - crc * sin2phi) ...
        * dphi;
i = i0+idot*tk	      + cic*cos2phi+cis*sin2phi;
di = 2 * (cis * cos2phi - cic * sin2phi) ...
        * dphi + idot;
x1 = cos(u)*r;
y1 = sin(u)*r;
if svprn>5%IGSO/MEO
    Omega = Omega0+(Omegadot-Omegae_dot)*tk-Omegae_dot*toe;
    dOmega = Omegadot - Omegae_dot;
    Omega = rem(Omega+2*pi,2*pi);

    satp(1,1) = x1*cos(Omega)-y1*cos(i)*sin(Omega);
    satp(2,1) = x1*sin(Omega)+y1*cos(i)*cos(Omega);
    satp(3,1) = y1*sin(i);
    xdash = r * cos(u);
        ydash = r * sin(u);

        dxdash = dr * cos(u) - r * sin(u) * du;
        dydash = dr * sin(u) + r * cos(u) * du;

        Vx = dxdash * cos(Omega) -  dydash * cos(i) * sin(Omega) ...
            + ydash * sin(Omega) * sin(i) * di - (xdash * sin(Omega) + ...
            ydash * cos(i) * cos(Omega)) * dOmega;

        Vy = dxdash * sin(Omega) + dydash * cos(i) * cos(Omega) - ...
            ydash * sin(i) *cos(Omega) * di + (xdash * cos(Omega) - ...
            ydash * cos(i) * sin(Omega)) * dOmega;

        Vz = dydash * sin(i) + ydash * cos(i) * di;
    
else%GEO
    Omega = Omega0 + (Omegadot-0)*tk - Omegae_dot*toe;
    dOmega = Omegadot-0;
    %Reduce to between 0 and 360 deg
    Omega = rem(Omega+2*pi,2*pi);
    % Compute satellite coordinates
    x_gk = x1 * cos(Omega) - y1 * cos(i)*sin(Omega);
    y_gk = x1 * sin(Omega) + y1 * cos(i)*cos(Omega);
    z_gk = y1 * sin(i);
    ang0=-5/180*pi;

    Rx=[1 0 0; 0 cos(ang0) sin(ang0);0 -sin(ang0) cos(ang0)];
    Rz=[cos(Omegae_dot*tk) sin(Omegae_dot*tk) 0;-sin(Omegae_dot*tk) cos(Omegae_dot*tk) 0;0 0 1];
    satp(1:3,1)=Rz*Rx*[x_gk,y_gk,z_gk]';
    xdash = r * cos(u);
        ydash = r * sin(u);

        dxdash = dr * cos(u) - r * sin(u) * du;
        dydash = dr * sin(u) + r * cos(u) * du;

        Vx_gk = dxdash * cos(Omega) -  dydash * cos(i) * sin(Omega) ...
            + ydash * sin(Omega) * sin(i) * di - (xdash * sin(Omega) + ...
            ydash * cos(i) * cos(Omega)) * dOmega;

        Vy_gk = dxdash * sin(Omega) + dydash * cos(i) * cos(Omega) - ...
            ydash * sin(i) *cos(Omega) * di + (xdash * cos(Omega) - ...
            ydash * cos(i) * sin(Omega)) * dOmega;

        Vz_gk = dydash * sin(i) + ydash * cos(i) * di;
        
        Vx= Vx_gk*cos(Omegae_dot*tk)+Vy_gk*cos(ang0)*sin(Omegae_dot*tk)+Vz_gk*sin(Omegae_dot*tk)*sin(ang0)...
            +Omegae_dot*(-x_gk*sin(Omegae_dot*tk)+y_gk*cos(Omegae_dot*tk)*cos(ang0)+z_gk*cos(Omegae_dot*tk)*sin(ang0));
        Vy= -Vx_gk*sin(Omegae_dot*tk)+Vy_gk*cos(Omegae_dot*tk)*cos(ang0)+Vz_gk*cos(Omegae_dot*tk)*sin(ang0)...
            +Omegae_dot*(-x_gk*cos(Omegae_dot*tk)-y_gk*sin(Omegae_dot*tk)*cos(ang0)-z_gk*sin(Omegae_dot*tk)*sin(ang0));
        Vz=-Vy_gk*sin(ang0)+Vz_gk*cos(ang0);
        
end
satp(4,1) = Vx;
satp(5,1) = Vy;
satp(6,1) = Vz;
tcorr=tcorr+dtr-tgd;
