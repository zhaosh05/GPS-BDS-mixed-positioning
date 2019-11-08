function [time,dayofyear,hour,minute,second, dt, satGPS, satGLN, satBDS, satGAL, satQZS, satSBA, eof, obsGPS,LLI_L1,obsGLN,obsGAL,obsBDS,LLI_B1,obsQZS,obsSBA]...
    = readobs(fid,ver,ObsGPS_types,ObsGLN_types,ObsGAL_types,ObsBDS_types,ObsQZS_types,ObsSBA_types)
% read rinex file and output time, obs and other measurements
time = 0;
dayofyear=0;
hour=0;
minute=0;
second=0;
dt = 0;
satGPS = zeros(32,1);%GPS obs sat table
satGLN = zeros(32,1);%glonass
satGAL = zeros(32,1);%galileo
satBDS = zeros(32,1);%beidou
satQZS = zeros(32,1);%qzss
satSBA = zeros(32,1);%sbas
NoSv=0;%total sv
NoGPS=0;%GPS sat
NoGLN=0;%glonass
NoGAL=0;%galilieo
NoBDS=0;%beidou
NoQZS=0;%qzss
NoSBA=0;%sbas
eof = 0;
switch ver
    case 2.0
        while 1
            lin = fgets(fid); %
            answer = strfind(lin,'END OF HEADER');
            if ~isempty(answer)
                lin = fgetl(fid);
            end
            if (feof(fid) == 1)
                eof = 1;
                break
            end
           if (length(lin)>=30)
            if (  (strcmp(lin(28:30),' 0 ') == 1))% (strcmp(lin(2),'0') == 1)  &     
              ll = length(lin)-2;
              if ll > 60, ll = 60; end
              linp = lin(1:ll);        
              %fprintf('%60s\n',linp);
              [year, lin] = strtok(lin);
%               year;
              [month, lin] = strtok(lin);
              [day, lin] = strtok(lin);
%               month;
%               day;
              [hour, lin] = strtok(lin);
              %hour
              [minute, lin] = strtok(lin);
              %minute
              [second, lin] = strtok(lin);
              %second
              [OK_flag, lin] = strtok(lin);
              h = str2double(hour)+str2double(minute)/60+str2double(second)/3600;
              jd = julday(str2double(year)+2000, str2double(month), str2double(day), h);
              [week, sow] = toGPStime(jd);

              time = sow;
              [NoSv, lin] = strtok(lin,'G');

              for k = 1:str2double(NoSv)
                 [sat, lin] = strtok(lin,'G');
                 prn(k) = str2double(sat);
              end

              sats = prn(:);
              dT = strtok(lin);
              if isempty(dT) == 0
                 dt = str2double(dT);
              end
              break

            end
          end
        end

        obsGPS=[];
        obsGLN=[];
        obsGAL=[];
        obsBDS=[];
        obsQZS=[];
        obsSBA=[];
    case {3.01 3.02}
        NoObsGPS_types = size(ObsGPS_types,2)/3;
        NoObsGLN_types = size(ObsGLN_types,2)/3;
        NoObsGAL_types = size(ObsGAL_types,2)/3;
        NoObsBDS_types = size(ObsBDS_types,2)/3;
        NoObsQZS_types = size(ObsQZS_types,2)/3;
        NoObsSBA_types = size(ObsSBA_types,2)/3;
        obsGPS=zeros(32,NoObsGPS_types);
        obsGLN=zeros(32,NoObsGLN_types);
        obsGAL=zeros(32,NoObsGAL_types);
        obsBDS=zeros(32,NoObsBDS_types);
        obsQZS=zeros(32,NoObsQZS_types);
        obsSBA=zeros(32,NoObsSBA_types);
        LLI_L1=zeros(32,1);
        LLI_B1=zeros(32,1);
        while 1
            lin = fgets(fid);
            answer = strfind(lin,'END OF HEADER');
            if ~isempty(answer)
                lin = fgetl(fid);
            end
            if (feof(fid) == 1)
                eof = 1;
                break
            end
            if (strcmp(lin(1),'>') == 1 && length(lin)>=32)%a complete epoch record line
                if ((strcmp(lin(32),'0') == 1))
                    year = str2double(lin(3:6));
                    month = str2double(lin(8:9));
                    day = str2double(lin(11:12));
                    hour = str2double(lin(14:15));
                    minute = str2double(lin(17:18));
                    second = str2double(lin(19:29));
                    OK_flag = str2double(lin(32));
                    NoSv=str2double(lin(33:35));
                    h = hour+minute/60+second/3600;
                    jd = julday(year, month, day, h);
                    dayofyear = julday(year,month,day,0)-julday(year,1,1,0)+1;
                    [week, sow] = toGPStime(jd);%week and second of week
                    time = sow;
                    if length(lin)>=56 % contains dt
                        dt = str2double(lin(42:end));
                    end
                    for i=1:NoSv
                        lin=fgets(fid);
                        switch lin(1)
                            case 'G'
                                NoGPS=NoGPS+1;
                                satGPS(NoGPS)=str2double(lin(2:3));
                                tt1=lin(4:end);%remove the sat # in the line head
                                for j=1:NoObsGPS_types
                                    if length(tt1)>=j*16
                                        obsGPS(NoGPS,j)=str2double(tt1((j-1)*16+1:(j-1)*16+14));
                                        if strcmp(ObsGPS_types((j-1)*3+1:j*3),'L1C')==1
                                            LLI_L1(NoGPS)=str2double(tt1((j-1)*16+15));
                                            if isnan(LLI_L1(NoGPS))
                                                LLI_L1(NoGPS)=0;
                                            end
                                        end
                                        %convert the blank measurement to zero
                                        if isnan(obsGPS(NoGPS,j))
                                            obsGPS(NoGPS,j)=0;
                                        end
                                    end
                                end
                            case 'R'%glonass
                                NoGLN=NoGLN+1;
                                satGLN(NoGLN)=str2double(lin(2:3));
                            case 'C'%bds
                                NoBDS=NoBDS+1;
                                satBDS(NoBDS)=str2double(lin(2:3));
                                tt1=lin(4:end);%remove the sat # in the line
                                for j=1:NoObsBDS_types
                                    if length(tt1)>=j*16
                                        %to cope with some uncompleted lines
                                        %every obs has 14bytes+2bytes including 3 after dot
                                        obsBDS(NoBDS,j)=str2double(tt1((j-1)*16+1:(j-1)*16+14));
                                        if strcmp(ObsBDS_types((j-1)*3+1:j*3),'L2I')==1
                                            LLI_B1(NoBDS)=str2double(tt1((j-1)*16+15));
                                            if isnan(LLI_B1(NoBDS))
                                                LLI_B1(NoBDS)=0;
                                            end
                                        end
                                        %convert the blank measurement to zero
                                        if isnan(obsBDS(NoBDS,j))
                                            obsBDS(NoBDS,j)=0;
                                        end
                                    end
                                end

                            case 'J'%qzss
                            case 'E'%galileo
                            case 'S'%sbas
                            otherwise
                                disp('error in obs data')
                        end
                    end
                end
              break
            end
        end
    otherwise
        disp('unidentified version')
end
