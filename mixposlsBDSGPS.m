function [pos, obsnewGPS, satsnewGPS, obsnewBDS, satsnewBDS,GDOP] = mixposlsBDSGPS(obsGPS,satsGPS,obsBDS,satsBDS,time,EphGPS,EphBDS,pos0,elthr,useConst,ecount)
%  full and fractional measurements mixed positioning

pos=zeros(5,1);
GDOP=0;
if pos0(1)==0
    pos=zeros(5,1);
else
    pos(1:4)=pos0(1:4);
end

v_light = 299792458;
lightms=v_light*1e-3;
dtr = pi/180;
mGPS = size(obsGPS,1);  % number of svs
mBDS = size(obsBDS,1);  % number of svs

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
nGEO=length(find(satsnewBDS<=5));
if nGEO<4
    pos=zeros(5,1);
    disp('insufficient GEO after eph!')
    return
end

nnoGEO=length(find(satsnewBDS>5))+mGPS;

pos=[pos(1:4);zeros(nnoGEO,1)];

% preliminary guess for receiver position and receiver clock offset
% pos = zeros(4,1);
no_iterations = 10; 
if mGPS+mBDS<4
    fprintf('only %d pr observables, ecount=%d\n', length(obs),ecount)
    pos=zeros(4,1);
    El = 0;
    GDOP=0;
    obsnewGPS=obs;
    satsnewGPS=sats;
    return
end

AGPS=zeros(mGPS,4+nnoGEO);
resGPS=zeros(mGPS,1);
ABDS=zeros(mBDS,4+nnoGEO);
resBDS=zeros(mBDS,1);

for iter = 1:no_iterations
    jj=0;
    for i=1:mBDS
        k = col_EphBDS(i);
        if satsnewBDS(i)>5
            jj=jj+1;
            tx_RAW = time-14 -(pos(4+jj)/v_light + obsnewBDS(i)/v_light); % 14s to compensate BDST and GPST
        else
            tx_RAW = time-14 - obsnewBDS(i)/v_light;
        end
        t0c = EphBDS(21,k);
        dt = check_t(tx_RAW-t0c);
        for j=1:3
            dt=dt-(EphBDS(2,k)*dt + EphBDS(20,k))*dt + EphBDS(19,k);
        end
        satcorr = (EphBDS(2,k)*dt + EphBDS(20,k))*dt + EphBDS(19,k)-EphBDS(18,k);
        tx = tx_RAW-satcorr;
        [XsatBDS(1:6,i), tcorrBDS(i)] = satposBD(tx, EphBDS(:,k));
        if iter == 1
            Rot_X = XsatBDS(1:3,i);
            rho2=norm(Rot_X-pos(1:3));
            varpBDS(i)=1;
            los=(Rot_X-pos(1:3))/rho2;

        else
            rho2 = (XsatBDS(1,i)-pos(1))^2+(XsatBDS(2,i)-pos(2))^2+(XsatBDS(3,i)-pos(3))^2;
            traveltime = sqrt(rho2)/v_light;
            Rot_X = corrotation(traveltime,XsatBDS(1:3,i));
            rho2=norm((Rot_X-pos(1:3)));
            los=(Rot_X-pos(1:3))/rho2;
            vlos=XsatBDS(4:6,i)'*los;

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
               el = 90;
            else
               el = atan2(U,hor_dis)/dtr;
            end
            El(i) = el;
        end

        ri=1/rho2;
        if satsnewBDS(i)<=5
            ABDS(i,1:4+nnoGEO) =[(-(Rot_X(1)-pos(1)))*ri...
                (-(Rot_X(2)-pos(2)))*ri ...
                (-(Rot_X(3)-pos(3)))*ri 1 zeros(1,nnoGEO)];
            resBDS(i,1) = obsnewBDS(i)-rho2-pos(4)+v_light*tcorrBDS(i);
        else
            tmp=zeros(1,nnoGEO);
            tmp(jj)=-1;%-lightms;
            ABDS(i,1:4+nnoGEO) =[(-(Rot_X(1)-pos(1)))*ri...
                (-(Rot_X(2)-pos(2)))*ri ...
                (-(Rot_X(3)-pos(3)))*ri 1 tmp];
            resBDS(i,1) = obsnewBDS(i)+pos(4+jj)-rho2-pos(4)+v_light*tcorrBDS(i);
        end
    end %i

    jj=0;
    for i=1:mGPS
        k = col_EphGPS(i);
        jj=jj+1;
        tx_RAW = time -(pos(4+mBDS-nGEO+jj)/v_light + obsnewGPS(i)/v_light);

        t0c = EphGPS(21,k);
        dt = check_t(tx_RAW-t0c);
        for j=1:3
            dt=dt-(EphGPS(2,k)*dt + EphGPS(20,k))*dt + EphGPS(19,k);
        end
        satcorr = (EphGPS(2,k)*dt + EphGPS(20,k))*dt + EphGPS(19,k)-EphGPS(18,k);
        tx = tx_RAW-satcorr;
        [XsatGPS(1:6,i), tcorrGPS(i)] = satposGPS(tx, EphGPS(:,k));
        if iter == 1
            Rot_X = XsatGPS(1:3,i);
            rho2=norm(Rot_X-pos(1:3));
            varpGPS(i)=1;
            los=(Rot_X-pos(1:3))/rho2;
       else
            rho2 = (XsatGPS(1,i)-pos(1))^2+(XsatGPS(2,i)-pos(2))^2+(XsatGPS(3,i)-pos(3))^2;
            traveltime = sqrt(rho2)/v_light;
            Rot_X = corrotation(traveltime,XsatGPS(1:3,i));
            rho2=norm((Rot_X-pos(1:3)));
            los=(Rot_X-pos(1:3))/rho2;
            vlos=XsatGPS(4:6,i)'*los;

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
               el = 90;
            else
               el = atan2(U,hor_dis)/dtr;
            end
            El(i) = el;
        end

        ri=1/rho2;
        
        tmp=zeros(1,nnoGEO);
        tmp(mBDS-nGEO+jj)=-1;
        AGPS(i,1:4+nnoGEO) =[(-(Rot_X(1)-pos(1)))*ri...
            (-(Rot_X(2)-pos(2)))*ri ...
            (-(Rot_X(3)-pos(3)))*ri 1 tmp];
        resGPS(i,1) = obsnewGPS(i)+pos(4+mBDS-nGEO+jj)-rho2-pos(4)+v_light*tcorrGPS(i);
    end %i
    varp=[varpBDS;varpGPS];
    A=[ABDS;AGPS];
    res=[resBDS;resGPS];
    P=diag(1./varp);
    
    AA0=A'*P*A;
    GDOP=sqrt(sum(1./eig(A(1:nGEO,1:4)'*A(1:nGEO,1:4))));
    if sum(GDOP)<2993.4 %eigenvalue test
        AA=AA0^-1;
    else
        pos=zeros(5,1);
        fprintf('eigenvalue too close to 0! GDOP=%.1f,ecount=%d\n',GDOP,ecount);
        return
    end
    
    x=AA*A'*P*res;

    pos(1:4) = pos(1:4)+x(1:4);
    pos(5:end)=pos(5:end)+round(x(5:end)/lightms)*lightms;
    if (sum(x(1:4).*x(1:4)))<1e-5
        obsBDSfull=[obsnewBDS(satsnewBDS<=5);pos(5:4+mBDS-nGEO)+obsnewBDS(satsnewBDS>5)];%recover full pseudorange
        obsGPSfull=pos(5+mBDS-nGEO:end)+obsnewGPS;
        % calculate final position using recovered full measurements
        posfull = poslsBDGPS(obsGPSfull,satsnewGPS,obsBDSfull,satsnewBDS,time,EphGPS,EphBDS,pos(1:4),elthr,'BDSGPS');
        break
    end
    
end % iter

pos=posfull(1:5);
