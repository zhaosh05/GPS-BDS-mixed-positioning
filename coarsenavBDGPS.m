%observation file and navigation message file should be RINEX 3.02 or 3.03
% obsfile = '.\ANMG20181001-1004.00o';
obsfile = '.\TOW200AUS20181001-1004.00o';%observation file name, can be changed
navfile = '.\brdm20181001-1004.18p';%navigation message file

[EphGPS,EphBDS,EphGLN,EphGAL,EphQZS]=readnav(navfile);

useConst='BDSGPS';%
v_light = 299792458;   % m/s
lightms=v_light*1e-3;

tic
time1 = -1*1.e50;
time2 = -1*1.e50;
[ObsGPS_types1,ObsGLN_types1,ObsGAL_types1,ObsBDS_types1,ObsQZS_types1,ObsSBA_type1, ant_delta1, ver1, filetype1, ifound_types1, eof11] = anheader(obsfile);

if ((ifound_types1 == 0) || (eof11 == 1))
   error('Basic information is missing in RINEX file')
end

fid1 = fopen(obsfile);

ecount = 0;%the count of epochs processed
stopcount=180;%5758;8640;%

posbase=zeros(5,stopcount);
posbase1=zeros(5,stopcount);
posbase2=zeros(5,stopcount);
timex=zeros(1,stopcount);
usedBDSGEO=zeros(5,stopcount);
GDOPBDSGEO=zeros(1,stopcount);

posb=zeros(4,1);
posb1=zeros(4,1);

while 1
    ifile_to_read = obsfile;
    [time1,dayofyear1,hour1,minute1,second1,dt1, satGPS1, satGLN1, satBDS1, satGAL1, satQZS1, satSBA1, eof1, obsGPS1,LLI1_L1,obsGLN1,obsGAL1,obsBDS1,LLI1_B1,obsQZS1,obsSBA1]=...
        readobs(fid1,ver1,ObsGPS_types1,ObsGLN_types1,ObsGAL_types1,ObsBDS_types1,ObsQZS_types1,ObsSBA_type1);
    if (eof1 == 1), break; end
    %GPS L1 observations
    NoSv1L1 = sum(satGPS1>0);
    
    if NoSv1L1>0
        noGPSflg=0;
        obs1L1=obsGPS1(1:NoSv1L1,:);
        sats1L1=satGPS1(1:NoSv1L1,:);
        LLI1L1=LLI1_L1(1:NoSv1L1);
        Obs_types1L1=ObsGPS_types1;
    else
        noGPSflg=1;
        if contains(useConst,'GPS')
            disp('no GPS obs in this epoch')
            continue
        end
    end
    %BDS B1 observations
    NoSv1B1 = sum(satBDS1>0);
    if NoSv1B1>0
        noBDSflg=0;
        obs1B1=obsBDS1(1:NoSv1B1,:);
        sats1B1=satBDS1(1:NoSv1B1,:);
        LLI1B1=LLI1_B1(1:NoSv1B1);
        Obs_types1B1=ObsBDS_types1;
    else
        noBDSflg=1;
        if contains(useConst,'BDS')
            disp('no BDS obs in this epoch')
            continue
        end
    end
    switch ver1
        case 2.10
        case {3.01 3.02}
            if ~noGPSflg
                indPR1L1=(strfind(Obs_types1L1,'C1C')+2)/3;
                indCP1L1=(strfind(Obs_types1L1,'L1C')+2)/3;
                indSNR1L1=(strfind(Obs_types1L1,'S1C')+2)/3;%SNR in rinex file
            end
            if ~noBDSflg
                indPR1B1=(strfind(Obs_types1B1,'C1I')+2)/3;%BDS B1
                indCP1B1=(strfind(Obs_types1B1,'L1I')+2)/3;
                indSNR1B1=(strfind(Obs_types1B1,'S1I')+2)/3;%SNR in rinex file
            end
        otherwise
            disp('cannot identify rinex version for base')
            return
    end
       
    ecount=ecount+1;
    timex(ecount)=time1;

    if ~noGPSflg
        indnonzero1=find(obs1L1(:,indPR1L1)~=0);
        if length(indnonzero1)<length(obs1L1(:,1))
            fprintf('defective GPS obs in rinex file, ecount=%d\n',ecount);
        end
        obs1L1=obs1L1(indnonzero1,:);
        sats1L1=sats1L1(indnonzero1,:);
        LLI1L1=LLI1L1(indnonzero1,:);
    end
    if ~noBDSflg
        indnonzero1=find(obs1B1(:,indPR1B1)~=0);
        if length(indnonzero1)<length(obs1B1(:,1))
            fprintf('defective BDS obs in rinex file, ecount=%d\n',ecount);
        end
        obs1B1=obs1B1(indnonzero1,:);
        sats1B1=sats1B1(indnonzero1,:);
        LLI1B1=LLI1B1(indnonzero1,:);
        
        temp=sortrows([sats1B1,obs1B1,LLI1B1],1);
        sats1B1=temp(:,1);
        obs1B1=temp(:,2:2+length(obs1B1(1,:))-1);
        LLI1B1=temp(:,3);
    end
%     posb(1:3)=[-2170120, 4385077,4078183]';%Weiqing
%     posb(1:3)=[ -2956913.1070  5075881.9300  2476417.1120]';%ANMG
    posb(1:3)=[ -5054582.6800  3275504.5720 -2091539.8860]';%TOW2 theoretical position as initial guess

    %% positioning using full measurements
    switch useConst
        case 'GPS'

        case 'BDS'

        case 'BDSGPS'
            %conventional dual-constellation positioning
            obs1L1new=obs1L1(:,indPR1L1);
            sats1L1new=sats1L1(:,:);
            obs1B1new=obs1B1(:,indPR1B1);
            sats1B1new=sats1B1(:,:);
            [posbBDSGPS, obs1newGPS, sats1newGPS, obs1newBDS, sats1newBDS] = poslsBDGPS(obs1L1new, sats1L1new,obs1B1new, sats1B1new, time1, EphGPS,EphBDS,posb,0,useConst);
            posbase(:,ecount)=posbBDSGPS;
    end


    %% reduce all pseudoranges to its 1ms fractional parts to test mixed positioning
    switch useConst
        case 'GPS'

        case 'BDS'
        
        case 'BDSGPS'
            posb(1:3)=zeros(1,3);
            obs1L1new=obs1L1(:,indPR1L1);
            sats1L1new=sats1L1(:,:);
            obsfracGPS=obs1L1new-round(obs1L1new/lightms)*lightms;
            
            obs1B1new=obs1B1(:,indPR1B1);
            sats1B1new=sats1B1(:,:);
            if length(find(sats1B1new<=5))<4 %number of GEO less than 4, exit
                fprintf('insufficient GEO, ecount=%d\n',ecount);
            else
                obs1B1new(sats1B1new>5)=obs1B1new(sats1B1new>5)-round(obs1B1new(sats1B1new>5)/lightms)*lightms;
            
                [posbBDSGPSGEO, obs1newGPS, sats1newGPS, obs1newBDS, sats1newBDS,GDOPBDSGEO(ecount)] = mixposlsBDSGPS(obsfracGPS, sats1L1new,obs1B1new, sats1B1new, time1, EphGPS,EphBDS,posbBDSGPS,0,useConst,ecount);
                posbase1(:,ecount)=posbBDSGPSGEO;
                usedBDSGEO(sats1newBDS(sats1newBDS<=5),ecount)=1;
            end
    end
    
    if ecount>=stopcount%
        break
    end

end
save('BDresult','posbase1','posbase','timex','usedBDSGEO','ecount','obsfile','navfile','GDOPBDSGEO')
toc
xlab=(timex-timex(1))/3600;%1/60:1/60:ecount/60;

xnonzero=find(posbase1(1,1:ecount)~=0);
figure

subplot(3,2,1)
stdxc=std(posbase1(1,xnonzero));
str1=['STD:',num2str(stdxc,'%.2f')];
plot(xlab(xnonzero),posbase1(1,xnonzero)-mean(posbase1(1,xnonzero)),'k.-')
grid on
axis tight
ylabel('X/m')
legend(str1)
subplot(3,2,3)
stdxc=std(posbase1(2,xnonzero));
str1=['STD:',num2str(stdxc,'%.2f')];
plot(xlab(xnonzero),posbase1(2,xnonzero)-mean(posbase1(2,xnonzero)),'r.-')
grid on
axis tight
ylabel('Y/m')
legend(str1)
subplot(3,2,5)
stdxc=std(posbase1(3,xnonzero));
str1=['STD:',num2str(stdxc,'%.2f')];
plot(xlab(xnonzero),posbase1(3,xnonzero)-mean(posbase1(3,xnonzero)),'b.-')
grid on
axis tight
ylabel('Z/m')
legend(str1)
xlabel('Time/hour')

subplot(3,2,2)
stdx=std(posbase(1,xnonzero));
str1=['STD:',num2str(stdx,'%.2f')];
plot(xlab(xnonzero),posbase(1,xnonzero)-mean(posbase(1,xnonzero)),'k.-')
grid on
axis tight
%         ylabel('X/m')
legend(str1)
subplot(3,2,4)
stdx=std(posbase(2,xnonzero));
str1=['STD:',num2str(stdx,'%.2f')];
plot(xlab(xnonzero),posbase(2,xnonzero)-mean(posbase(2,xnonzero)),'r.-')
grid on
axis tight
%         ylabel('Y/m')
legend(str1)
subplot(3,2,6)
stdx=std(posbase(3,xnonzero));
str1=['STD:',num2str(stdx,'%.2f')];
plot(xlab(xnonzero),posbase(3,xnonzero)-mean(posbase(3,xnonzero)),'b.-')
grid on
axis tight
%         ylabel('Z/m')
legend(str1)
xlabel('Time/hour')

fclose all;
