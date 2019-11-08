function [EphGPS,EphBDS,EphGLN,EphGAL,EphQZS]=readnav(navfile)
% read rinex ephemeris message file
EphGPS=0;
EphBDS=0;
EphGLN=0;
EphGAL=0;
EphQZS=0;

fide = fopen(navfile);
head_lines = 0;
talfa0=0;
talfa1=0;
talfa2=0;
talfa3=0;
tbeta0=0;
tbeta1=0;
tbeta2=0;
tbeta3=0;

while 1  % skip header
   head_lines = head_lines+1;
   line = fgetl(fide);
   temp=strfind(line,'RINEX VERSION / TYPE');
   if ~isempty(temp)
       [ver,remain]=strtok(line);
       ver=str2double(ver);
   end
   if ver>=3
       temp=strfind(line,'IONOSPHERIC CORR');
       if ~isempty(temp)
           if strcmp(line(1:4),'GPSA')==1
               talfa0=str2double(line(6:17));
               talfa1=str2double(line(18:29));
               talfa2=str2double(line(30:41));
               talfa3=str2double(line(42:51));
           elseif strcmp(line(1:4),'GPSB')==1
               tbeta0=str2double(line(6:17));
               tbeta1=str2double(line(18:29));
               tbeta2=str2double(line(30:41));
               tbeta3=str2double(line(42:51));
           end
       end
   elseif ver>=2
       temp=strfind(line,'ION ALPHA');
       if ~isempty(temp)
           talfa0=str2double(line(3:14));
           talfa1=str2double(line(15:26));
           talfa2=str2double(line(27:38));
           talfa3=str2double(line(39:50));
       end
       temp=strfind(line,'ION BETA');
       if ~isempty(temp)
           tbeta0=str2double(line(3:14));
           tbeta1=str2double(line(15:26));
           tbeta2=str2double(line(27:38));
           tbeta3=str2double(line(39:50));
       end
   end
   answer = strfind(line,'END OF HEADER');
   if ~isempty(answer), break;	end
end
% head_lines
offset=ftell(fide);
noeph = -1;
while 1
   noeph = noeph+1;
   line = fgetl(fide);
   if line == -1, break;  end
end
noeph = round(noeph/8);
%index to symbolize which line of observations are for the specific system
GPSind=zeros(1,noeph);
BDSind=zeros(1,noeph);
GALind=zeros(1,noeph);
GLNind=zeros(1,noeph);
QZSind=zeros(1,noeph);
fseek(fide,offset,'bof');

svprn	 = zeros(1,noeph);
weekno	 = zeros(1,noeph);
t0c	 = zeros(1,noeph);
tgd	 = zeros(1,noeph);
aodc	 = zeros(1,noeph);
toe	 = zeros(1,noeph);
af2	 = zeros(1,noeph);
af1	 = zeros(1,noeph);
af0	 = zeros(1,noeph);
aode	 = zeros(1,noeph);
deltan	 = zeros(1,noeph);
M0	 = zeros(1,noeph);
ecc	 = zeros(1,noeph);
roota	 = zeros(1,noeph);
toe	 = zeros(1,noeph);
cic	 = zeros(1,noeph);
crc	 = zeros(1,noeph);
cis	 = zeros(1,noeph);
crs	 = zeros(1,noeph);
cuc	 = zeros(1,noeph);
cus	 = zeros(1,noeph);
Omega0	 = zeros(1,noeph);
omega	 = zeros(1,noeph);
i0	 = zeros(1,noeph);
Omegadot = zeros(1,noeph);
idot	 = zeros(1,noeph);
accuracy = zeros(1,noeph);
health	 = zeros(1,noeph);
fit	 = zeros(1,noeph);
alfa0 = zeros(1,noeph);
alfa1 = zeros(1,noeph);
alfa2 = zeros(1,noeph);
alfa3 = zeros(1,noeph);
beta0 =  zeros(1,noeph);
beta1 =  zeros(1,noeph);
beta2 =  zeros(1,noeph);
beta3 =  zeros(1,noeph);

switch ver
    case {2.10 2.11}% only applicable to gps nav
        for i = 1:noeph
            line = fgetl(fide);	  %
            svprn(i) = str2double(line(1:2));
            year = line(3:6);
            month = line(7:9);
            day = line(10:12);
            hour = line(13:15);
            minute = line(16:18);
            second = line(19:22);
            af0(i) = str2double(line(23:41));
            af1(i) = str2double(line(42:60));
            af2(i) = str2double(line(61:79));
            line = fgetl(fide);	  %
            IODE = line(4:22);
            crs(i) = str2double(line(23:41));
            deltan(i) = str2double(line(42:60));
            M0(i) = str2double(line(61:79));
            line = fgetl(fide);	  %
            cuc(i) = str2double(line(4:22));
            ecc(i) = str2double(line(23:41));
            cus(i) = str2double(line(42:60));
            roota(i) = str2double(line(61:79));
            line=fgetl(fide);
            toe(i) = str2double(line(4:22));
            cic(i) = str2double(line(23:41));
            Omega0(i) = str2double(line(42:60));
            cis(i) = str2double(line(61:79));
            line = fgetl(fide);	    %
            i0(i) =  str2double(line(4:22));
            crc(i) = str2double(line(23:41));
            omega(i) = str2double(line(42:60));
            Omegadot(i) = str2double(line(61:79));
            line = fgetl(fide);	    %
            idot(i) = str2double(line(4:22));
            codes = str2double(line(23:41));
            weekno = str2double(line(42:60));
            L2flag = str2double(line(61:79));
            line = fgetl(fide);	    %
            svaccur = str2double(line(4:22));
            svhealth = str2double(line(23:41));
            tgd(i) = str2double(line(42:60));
            iodc = line(61:79);
            line = fgetl(fide);	    %
            tom(i) = str2double(line(4:22));
            alfa0(i) = talfa0;
            alfa1(i) = talfa1;
            alfa2(i) = talfa2;
            alfa3(i) = talfa3;
            beta0(i) = tbeta0;
            beta1(i) = tbeta1;
            beta2(i) = tbeta2;
            beta3(i) = tbeta3;
        end
        eph(1,:)  = svprn;
        eph(2,:)  = af2;
        eph(3,:)  = M0;
        eph(4,:)  = roota;
        eph(5,:)  = deltan;
        eph(6,:)  = ecc;
        eph(7,:)  = omega;
        eph(8,:)  = cuc;
        eph(9,:)  = cus;
        eph(10,:) = crc;
        eph(11,:) = crs;
        eph(12,:) = i0;
        eph(13,:) = idot;
        eph(14,:) = cic;
        eph(15,:) = cis;
        eph(16,:) = Omega0;
        eph(17,:) = Omegadot;
        eph(18,:) = tgd;%toe;modified by zsh
        eph(19,:) = af0;
        eph(20,:) = af1;
        eph(21,:) = toe;

        eph(22,:) =alfa0;
        eph(23,:) =alfa1;
        eph(24,:) =alfa2;
        eph(25,:) =alfa3;
        eph(26,:) =beta0;
        eph(27,:) =beta1;
        eph(28,:) =beta2;
        eph(29,:) =beta3;
        
        EphGPS=eph;

    case {3.01 3.02 3.03}
        for i = 1:noeph
            line = fgetl(fide);	  %
            switch line(1)
                case 'G'
                    GPSind(i)=1;
                    svprn(i) = str2double(line(2:3));
                    year = line(5:8);
                    month = line(10:11);
                    day = line(13:14);
                    hour = line(16:17);
                    minute = line(19:20);
                    second = line(22:23);
                    af0(i) = str2double(line(24:42));
                    af1(i) = str2double(line(43:61));
                    af2(i) = str2double(line(62:80));

                    line = fgetl(fide);	  %broadcast orbit 1
                    IODE = line(5:23);
                    crs(i) = str2double(line(24:42));
                    deltan(i) = str2double(line(43:61));
                    M0(i) = str2double(line(62:80));
                    line = fgetl(fide);	  %broadcast orbit 2
                    cuc(i) = str2double(line(5:23));
                    ecc(i) = str2double(line(24:42));
                    cus(i) = str2double(line(43:61));
                    roota(i) = str2double(line(62:80));
                    line=fgetl(fide); %broadcast orbit 3
                    toe(i) = str2double(line(5:23));
                    cic(i) = str2double(line(24:42));
                    Omega0(i) = str2double(line(43:61));
                    cis(i) = str2double(line(62:80));
                    line = fgetl(fide);	    %broadcast orbit 4
                    i0(i) =  str2double(line(5:23));
                    crc(i) = str2double(line(24:42));
                    omega(i) = str2double(line(43:61));
                    Omegadot(i) = str2double(line(62:80));
                    line = fgetl(fide);	    %broadcast orbit 5
                    idot(i) = str2double(line(5:23));
                    codes = str2double(line(24:42));
                    weekno = str2double(line(43:61));
                    L2flag = str2double(line(62:80));
                    line = fgetl(fide);	    %broadcast orbit 6
                    svaccur = str2double(line(5:23));
                    svhealth = str2double(line(24:42));
                    
                    %remove unhealthy sat
                    if svhealth ==1
                        GPSind(i)=0;
                    end
                    tgd(i) = str2double(line(43:61));
                    iodc = line(62:80);
                    line = fgetl(fide);	    %broadcast orbit 7
                    tom(i) = str2double(line(5:23));

                    alfa0(i) = talfa0;
                    alfa1(i) = talfa1;
                    alfa2(i) = talfa2;
                    alfa3(i) = talfa3;
                    beta0(i) = tbeta0;
                    beta1(i) = tbeta1;
                    beta2(i) = tbeta2;
                    beta3(i) = tbeta3;
                case 'C'
                    BDSind(i)=1;
                    svprn(i) = str2double(line(2:3));
                    year = line(5:8);
                    month = line(10:11);
                    day = line(13:14);
                    hour = line(16:17);
                    minute = line(19:20);
                    second = line(22:23);
                    af0(i) = str2double(line(24:42));
                    af1(i) = str2double(line(43:61));
                    af2(i) = str2double(line(62:80));

                    line = fgetl(fide);	  %broadcast orbit 1
                    IODE = line(5:23);
                    crs(i) = str2double(line(24:42));
                    deltan(i) = str2double(line(43:61));
                    M0(i) = str2double(line(62:80));
                    line = fgetl(fide);	  %broadcast orbit 2
                    cuc(i) = str2double(line(5:23));
                    ecc(i) = str2double(line(24:42));
                    cus(i) = str2double(line(43:61));
                    roota(i) = str2double(line(62:80));
                    line=fgetl(fide); %broadcast orbit 3
                    toe(i) = str2double(line(5:23));
                    cic(i) = str2double(line(24:42));
                    Omega0(i) = str2double(line(43:61));
                    cis(i) = -str2double(line(62:80));
                    line = fgetl(fide);	    %broadcast orbit 4
                    i0(i) =  str2double(line(5:23));
                    crc(i) = str2double(line(24:42));
                    omega(i) = str2double(line(43:61));
                    Omegadot(i) = str2double(line(62:80));
                    line = fgetl(fide);	    %broadcast orbit 5
                    idot(i) = str2double(line(5:23));
                    spare = str2double(line(24:42));
                    weekno = str2double(line(43:61));
                    spare = str2double(line(62:80));
                    line = fgetl(fide);	    %broadcast orbit 6
                    svaccur = str2double(line(5:23));
                    svhealth = str2double(line(24:42));
                    %remove unhealthy sat
                    if svhealth ==1
                        BDSind(i)=0;
                    end
                    tgd(i) = str2double(line(43:61));
                    tgd2 = line(62:80);
                    line = fgetl(fide);	    %broadcast orbit 7
                    tom(i) = str2double(line(5:23));

                    alfa0(i) = talfa0;
                    alfa1(i) = talfa1;
                    alfa2(i) = talfa2;
                    alfa3(i) = talfa3;
                    beta0(i) = tbeta0;
                    beta1(i) = tbeta1;
                    beta2(i) = tbeta2;
                    beta3(i) = tbeta3;
                case 'E'
                    GALind(i)=1;
                case 'R'
                    GLNind(i)=1;
                case 'J'
                    QZSind(i)=1;
                otherwise
                    disp('unidentified system in readnav')
            end
            
        %    spare = line(23:41);
        %    spare = line(42:60);
        %    spare = line(61:79);
        end
        eph(1,:)  = svprn;
        eph(2,:)  = af2;
        eph(3,:)  = M0;
        eph(4,:)  = roota;
        eph(5,:)  = deltan;
        eph(6,:)  = ecc;
        eph(7,:)  = omega;
        eph(8,:)  = cuc;
        eph(9,:)  = cus;
        eph(10,:) = crc;
        eph(11,:) = crs;
        eph(12,:) = i0;
        eph(13,:) = idot;
        eph(14,:) = cic;
        eph(15,:) = cis;
        eph(16,:) = Omega0;
        eph(17,:) = Omegadot;
        eph(18,:) = tgd;%toe;modified by zsh
        eph(19,:) = af0;
        eph(20,:) = af1;
        eph(21,:) = toe;

        eph(22,:) =alfa0;
        eph(23,:) =alfa1;
        eph(24,:) =alfa2;
        eph(25,:) =alfa3;
        eph(26,:) =beta0;
        eph(27,:) =beta1;
        eph(28,:) =beta2;
        eph(29,:) =beta3;

        EphGPS=eph(:,GPSind==1);
        EphBDS=eph(:,BDSind==1);
    otherwise
        disp('unidentified rinex version')
end


status = fclose(fide);
disp('Ephemerides decoded successfully')
fclose all;