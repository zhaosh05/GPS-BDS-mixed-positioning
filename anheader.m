function [ObsGPS_types,ObsGLN_types,ObsGAL_types,ObsBDS_types,ObsQZS_types,ObsSBA_types, ant_delta, ver, filetype, ifound_types, eof] = anheader(file)
% analyze observation file header
% reference: RINEX 3.02

fid = fopen(file,'rt');
eof = 0;
ifound_types = 0;
Obs_types = [];
ant_delta = [];
ver=[];
filetype=[];
ObsGPS_types=[];
ObsGLN_types=[];
ObsGAL_types=[];
ObsBDS_types=[];
ObsQZS_types=[];
ObsSBA_types=[];
NobsGPS=0;
NobsGLN=0;
NobsGAL=0;
NobsBDS=0;
NobsQZS=0;
NobsSBA=0;
while 1			   % Gobbling the header
   line = fgetl(fid);
   answer = strfind(line,'END OF HEADER');
   if  ~isempty(answer), break; end
   if (line == -1), eof = 1; break; end
   answer = strfind(line,'RINEX VERSION / TYPE');
   if ~isempty(answer)
      ver=str2double(line(1:9));
      filetype = line(41);
      
   end
   
   answer = strfind(line,'ANTENNA: DELTA H/E/N');
   if ~isempty(answer)
      for k = 1:3
         [delta, line] = strtok(line);
         del = str2double(delta);
         ant_delta = [ant_delta del];
      end
   end
   switch ver
       case {2.10}
           answer = strfind(line,'# / TYPES OF OBSERV');
           if ~isempty(answer)
              [NObs, line] = strtok(line);
              NoObs = str2double(NObs);
              for k = 1:NoObs
                 [ot, line] = strtok(line);
                 ObsGPS_types = [Obs_types ot];
              end
              ifound_types = 1;
           end
       case{3.01 3.02}
           answer = strfind(line,'SYS / # / OBS TYPES');
           if ~isempty(answer)
               ifound_types = 1;
               switch line(1)
                   case 'G'
                       type=line(1);
                       NobsGPS=str2double(line(4:6));
                       line=line(8:58);
                       ncount=0;
                       while ~isempty(line)% sometimes, there are multiple lines
%                        for k = 1:NobsGPS
                           [ot, line] = strtok(line);
                           ObsGPS_types = [ObsGPS_types ot];
                           ncount=ncount+1; % record the numbers in the first line
                       end

                   case 'R'
                       type=line(1);
                   case 'E'
                       type=line(1);
                   case 'J'
                       type=line(1);
                   case 'C'
                       type=line(1);
                       NobsBDS=str2double(line(4:6));
                       line=line(8:58);
                       ncount=0;
                       while ~isempty(line)% sometimes, there are multiple lines
                           [ot, line] = strtok(line);
                           ObsBDS_types = [ObsBDS_types ot];
                           ncount=ncount+1;
                       end
                   case 'S'
                       type=line(1);
                   case 'M'
                       type=line(1);
                   case ' '
                       % a continuing line for a certain type
                       switch type
                           case 'G'
                               line=line(8:58);
                               for k = 1:NobsGPS-ncount
                                   [ot, line] = strtok(line);
                                   ObsGPS_types = [ObsGPS_types ot];
                               end
                           case 'R'
                           case 'E'
                           case 'J'
                           case 'C'
                               line=line(8:58);
                               for k = 1:NobsBDS-ncount 
                                   [ot, line] = strtok(line);
                                   ObsBDS_types = [ObsBDS_types ot];
                               end
                           case 'S'
                           case 'M'
                           otherwise
                               disp('unidentified obs type in observation header')
                       end
                   otherwise
                       disp('unidentified obs type in observation header')
               end
           end
       otherwise
           disp('unidentified version')
   end
end
fclose(fid);

