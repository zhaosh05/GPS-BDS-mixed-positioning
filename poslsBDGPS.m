function [pos, obsnewGPS, satsnewGPS, obsnewBDS, satsnewBDS] = poslsBDGPS(obsGPS,satsGPS,obsBDS,satsBDS,time,EphGPS,EphBDS,pos0,elthr,useConst)
% least square iteration to solve position with two constellations
pos=zeros(5,1);%x, y, z, dt to reference system, dtsys between two constellations
if pos0(1)~=0
    pos(1:4,1)=pos0(1:4);
end

v_light = 299792458;
dtr = pi/180;
mGPS = size(obsGPS,1);  % number of GPS svs
mBDS = size(obsBDS,1);  % number of BDS svs

% identify ephemerides columns in Eph
col_EphGPS=zeros(1,mGPS);
for t = 1:mGPS
    col_EphGPS(t) = find_eph(EphGPS,satsGPS(t),time);
end
includeind= col_EphGPS~=0;
satsnewGPS=satsGPS(includeind);
obsnewGPS=obsGPS(includeind,:);
col_EphGPS=col_EphGPS(includeind);
mGPS=size(obsnewGPS,1);
ElGPS = zeros(1,mGPS);
varpGPS=zeros(mGPS,1);
XsatGPS=zeros(3,mGPS);
tcorrGPS=zeros(1,mGPS);

col_EphBDS=zeros(1,mBDS);
for t = 1:mBDS
    col_EphBDS(t) = find_eph(EphBDS,satsBDS(t),time-14);%BDS leap second 14s
end
includeind= col_EphBDS~=0;
satsnewBDS=satsBDS(includeind);
obsnewBDS=obsBDS(includeind,:);
col_EphBDS=col_EphBDS(includeind);
mBDS=size(obsnewBDS,1);
ElBDS = zeros(1,mBDS);
varpBDS=zeros(mBDS,1);
XsatBDS=zeros(3,mBDS);
tcorrBDS=zeros(1,mBDS);

no_iterations = 10; 

for i = 1:mGPS
    k = col_EphGPS(i);
    tx_RAW = time - obsnewGPS(i)/v_light;
    t0c = EphGPS(21,k);
    dt = check_t(tx_RAW-t0c);
    for j=1:3
        dt=dt-(EphGPS(2,k)*dt + EphGPS(20,k))*dt + EphGPS(19,k);
    end
    satcorr = (EphGPS(2,k)*dt + EphGPS(20,k))*dt + EphGPS(19,k)-EphGPS(18,k);
    tx = tx_RAW-satcorr;
    [XsatGPS(1:6,i), tcorrGPS(i)] = satposGPS(tx, EphGPS(:,k));
end
for i = 1:mBDS
    k = col_EphBDS(i);
    tx_RAW = time-14 - obsnewBDS(i)/v_light;
    t0c = EphBDS(21,k);
    dt = check_t(tx_RAW-t0c);
    for j=1:3
        dt=dt-(EphBDS(2,k)*dt + EphBDS(20,k))*dt + EphBDS(19,k);
    end
    satcorr = (EphBDS(2,k)*dt + EphBDS(20,k))*dt + EphBDS(19,k)-EphBDS(18,k);
    tx = tx_RAW-satcorr;
    [XsatBDS(1:6,i), tcorrBDS(i)] = satposBD(tx, EphBDS(:,k));
end
if mGPS+mBDS<4
    fprintf('only %d pr observables, cannot compute position\n', length(obs))
    pos=zeros(4,1);
    El = 0;
    GDOP=0;
    obsnewGPS=obs;
    satsnewGPS=sats;

    return
end
   
AGPS=zeros(mGPS,5);
resGPS=zeros(mGPS,1);
ABDS=zeros(mBDS,5);
resBDS=zeros(mBDS,1);
for iter = 1:no_iterations
    for i=1:mGPS
        if iter == 1
            Rot_X = XsatGPS(1:3,i);
            varpGPS(i)=1;
        else
            rho2 = (XsatGPS(1,i)-pos(1))^2+(XsatGPS(2,i)-pos(2))^2+(XsatGPS(3,i)-pos(3))^2;
            traveltime = sqrt(rho2)/v_light;
            Rot_X = corrotation(traveltime,XsatGPS(1:3,i));

            dx=Rot_X-pos(1:3,:);
            [phi,lambda,h] = togeod(6378137,298.257223563,pos(1),pos(2),pos(3));
            cl = cos(lambda*dtr); sl = sin(lambda*dtr);
            cb = cos(phi*dtr); sb = sin(phi*dtr);
            F = [-sl -sb*cl cb*cl;
                  cl -sb*sl cb*sl;
                   0    cb   sb];
            local_vector = F'*dx;
            E = local_vector(1);
            N = local_vector(2);
            U = local_vector(3);
            hor_dis = sqrt(E^2+N^2);
            if hor_dis < 1.e-20
               az = 0;
               el = 90;
            else
               az = atan2(E,N)/dtr;
               el = atan2(U,hor_dis)/dtr;
            end
            ElGPS(i) = el;
        end
        distGPS=norm(Rot_X-pos(1:3));
        resGPS(i,1) = obsnewGPS(i)-distGPS-pos(4)+v_light*tcorrGPS(i);
        ri=1/distGPS;
        AGPS(i,1:5) =[(-(Rot_X(1)-pos(1)))*ri...
                (-(Rot_X(2)-pos(2)))*ri ...
                (-(Rot_X(3)-pos(3)))*ri 1 0];
    end %i
    for i=1:mBDS
        if iter == 1
%                 traveltime = 0.092;
            Rot_X = XsatBDS(1:3,i);
            varpBDS(i)=1;
        else
            rho2 = (XsatBDS(1,i)-pos(1))^2+(XsatBDS(2,i)-pos(2))^2+(XsatBDS(3,i)-pos(3))^2;
            traveltime = sqrt(rho2)/v_light;
            Rot_X = corrotation(traveltime,XsatBDS(1:3,i));

            dx=Rot_X-pos(1:3,:);
            [phi,lambda,h] = togeod(6378137,298.257223563,pos(1),pos(2),pos(3));
            cl = cos(lambda*dtr); sl = sin(lambda*dtr);
            cb = cos(phi*dtr); sb = sin(phi*dtr);
            F = [-sl -sb*cl cb*cl;
                  cl -sb*sl cb*sl;
                   0    cb   sb];
            local_vector = F'*dx;
            E = local_vector(1);
            N = local_vector(2);
            U = local_vector(3);
            hor_dis = sqrt(E^2+N^2);
            if hor_dis < 1.e-20
               az = 0;
               el = 90;
            else
               az = atan2(E,N)/dtr;
               el = atan2(U,hor_dis)/dtr;
            end
            ElBDS(i) = el;
        end

        distBDS=norm(Rot_X-pos(1:3));
        resBDS(i,1) = obsnewBDS(i)-distBDS-pos(5)+v_light*tcorrBDS(i);
        ri=1/distBDS;
        ABDS(i,1:5) =[(-(Rot_X(1)-pos(1)))*ri...
                (-(Rot_X(2)-pos(2)))*ri ...
                (-(Rot_X(3)-pos(3)))*ri 0 1];
    end %i
    varp=[varpGPS;varpBDS];
    A=[AGPS;ABDS];
    res=[resGPS;resBDS];

    P=diag(1./varp);
    AA=(A'*P*A)^-1;
    x=AA*A'*P*res;
    pos = pos+x;
    if sum(x.*x)<1e-6
        break
    end
end % iter

