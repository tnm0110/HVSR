clear
%close all

% This is the main script that runs on daylong sac files (nodal seismomter) organised by event day direcotry (e.g. 2020-04-15) to calculate 
% horizontal to vertical spectral ratio (HVSR) from ambient noise. During the processing steps, it revoves the transient short time spikes using 
% STA/LTA and DBSCAN algorithom. 
% 

tic
%% stettin I/0
disp('setting I/O ...')
% STA LTA value in seconds 
sta=0.5; lta=60;
% STA LTA threshold range
r_min = 0.5; r_max = 1.50;
sdir='/space/mibhk/muse/noise/data2';   % input directory
addpath(sdir)
cd(sdir)
unix('ls -d 2020-04-15* > ndir.txt');

fid=fopen('ndir.txt','r');
ndir = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
ndir=ndir{1};
%%
% directory loop 
for m=1:length(ndir)
    dd=ndir{m};
    cd (dd)
    disp(['calculting day ', dd, ' ...']);
    %[~, files] = system('ls *SAC');
    unix('ls *X1.*SAC > sfile.txt');    % list all the SAC files inside the directory 
    fid=fopen('sfile.txt','r');
    files=textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    files=files{1};
    files = char(files);
    stn=files(:,22:24); 
    stn = unique(stn, 'rows');

% station loop
for k=1:size(stn,1)
%for k=1:2   % run only first five staitons 
    disp(['doing station ', stn(k,:), '... '])
    ssdir=pwd;
    % Use the dir function to list files with the .SAC extension
    mm=['*.' stn(k,:) '.*.SAC'];
    fileList = dir(fullfile(ssdir, mm));
    % Extract the names of the matching files
    fileNames = {fileList.name};

% inside staiton loop
%% remove the transient and split into hours
for i=1:length(fileNames)
    sfile=fileNames{i};
    cmp=sfile(26:28);
    yy=sfile(1:4); jd=sfile(6:8); s_name=sfile(22:24);
    sac=rdsac(sfile);
    f=1/sac.HEADER.DELTA;

    if strcmp(cmp,'DPZ')
        disp(['calculating noise window :: ==> station :: ', s_name,...
            ' ==> ', cmp '.... '])
        [amp_noise_dpz]=sta_lta_f(sac,sta,lta,r_min,r_max);
    elseif strcmp(cmp,'DP1')
        disp(['calculating noise window :: ==> station :: ', s_name,...
            ' ==> ', cmp '.... '])
        [amp_noise_dp1]=sta_lta_f(sac,sta,lta,r_min,r_max);
    else
        disp(['calculating noise window :: ==> station :: ', s_name,...
            ' ==> ', cmp '.... '])
        [amp_noise_dp2]=sta_lta_f(sac,sta,lta,r_min,r_max);
    end
end

%% 
    ntps=min([length(amp_noise_dp1),length(amp_noise_dp2),...
        length(amp_noise_dpz)]);
    wind_len=ceil(3600*f);
    wind_n=ceil(ntps/wind_len)-1;
    
    % saving the hour window
    disp(['saving the hour window ', s_name, ' ==> ', cmp '.... '])
    wind_hr=zeros(wind_len,3,wind_n);
    for j=1:wind_n
        st=(j-1)*wind_len +1;
        ed=j*wind_len;
        wind_cmp(:,1)=amp_noise_dpz(st:ed);
        wind_cmp(:,2)=amp_noise_dp1(st:ed);
        wind_cmp(:,3)=amp_noise_dp2(st:ed);
        wind_hr(:,:,j)=wind_cmp;
        
        wind_name=[yy '_' jd '_' s_name '.mat' ];
        %save(wind_name,"wind_hr"); %saving the window

% test plot
% figure(1)
% subplot(3,1,1)
% plot((1:size(wind_hr,1))/f,wind_hr(:,1),'LineWidth',0.8)
% 
% subplot(3,1,2)
% plot((1:size(wind_hr,1))/f,wind_hr(:,2),'LineWidth',0.8)
% 
% subplot(3,1,3)
% plot((1:size(wind_hr,1))/f,wind_hr(:,3),'LineWidth',0.8)
    end  
toc

%% calculating HVSR
disp(['calculating hvsr ', s_name, ' ==> ', cmp '.... '])
fs=f/2; % sampling frequency
hvsr_name=['hsvr_', yy '_' jd '_' s_name '.mat' ];
hvsr=zeros(size(wind_hr,1),3,size(wind_hr,3));

for kk=1:size(wind_hr,3)
    dpz=wind_hr(:,1,kk);
    dp1=wind_hr(:,2,kk);
    dp2=wind_hr(:,3,kk);
    
    % doing fourier transform
    f_dpz=fft(dpz);
    f_dp1=fft(dp1);
    f_dp2=fft(dp2);

    N=length(dpz);
    freq = (0:N-1) * (fs/N); % Frequency vector

    % calculate the fourier amplitude of the spectrum
    f_dpz_mag = abs(f_dpz) / N;
    f_dp1_mag = abs(f_dp1) / N;
    f_dp2_mag = abs(f_dp2) / N;

    % getting average horizonatal components
    %amp_hz=(amp_dp1+amp_dp2);
    %f_hz_mag=(abs(f_dp1_mag) + abs(f_dp2_mag));
    f_hz_mag=sqrt((f_dp1_mag).^2 + (f_dp2_mag).^2);
    
    % Apply moving average for smoothing
    wsize=1000;
    amp_dpz_sm = movmean(f_dpz_mag, wsize);
    amp_dp1_sm = movmean(f_dp1_mag, wsize);
    amp_dp2_sm = movmean(f_dp2_mag, wsize);
    amp_hz_sm = movmean(f_hz_mag, wsize);


% Apply median filtering to remove the spike

% windowSize = 50;  % Choose an appropriate window size
% amp_hz_sm = medfilt1(amp_hz_sm, windowSize);
% amp_dp1_sm = medfilt1(amp_dp1_sm, windowSize);
% amp_dp2_sm = medfilt1(amp_dp2_sm, windowSize);
% amp_dpz_sm = medfilt1(amp_dpz_sm, windowSize);

    %calculate spectral ratio
    hvsr_d1=amp_dp1_sm./amp_dpz_sm;
    hvsr_hz=amp_hz_sm./amp_dpz_sm;
    hvsr_d2=amp_dp2_sm./amp_dpz_sm;
    hvsr(:,:,kk)=[hvsr_hz hvsr_d1 hvsr_d2];
    
end
save(hvsr_name,"hvsr")

end     % station loop ends
cd ../
end     % dir loop ends

disp('All Finished')

% %Plot the spectral amplitude vs. frequency
% figure(1);
% for ff=1:size(hvsr,3)
% semilogx(freq, hvsr(:,1,ff));
% hold on
% end
% xlabel('Frequency (Hz)');
% ylabel('HVSR');    
% title('Spectral Amplitude (Horizontal) vs. Frequency');
% grid on;
% xlim([0.1, 10]);
% 
% figure(2);
% for ff=1:size(hvsr,3)
% semilogx(freq, hvsr(:,2,ff));
% hold on
% end
% xlabel('Frequency (Hz)');
% ylabel('HVSR');    
% title('Spectral Amplitude (DP1) vs. Frequency');
% grid on;
% xlim([0.1, 10]);
% 
% figure(3);
% for ff=1:size(hvsr,3)
% semilogx(freq, hvsr(:,3,ff));
% hold on
% end
% xlabel('Frequency (Hz)');
% ylabel('HVSR');    
% title('Spectral Amplitude (DP2) vs. Frequency');
% grid on;
% xlim([0.1, 10]);
% 
% % xlim([0.1, 10]);

% 
% figure(2)
% semilogx(freq,hvsr_hz,'r')
% hold on
% semilogx(freq,hvsr_d1,'b')
% hold on
% semilogx(freq,spec_hvsr)
% 
% xlabel('Frequency (Hz)');
% ylabel('Amplitude Ratio');    
% grid on;
% xlim([0.1, 20]);
% 
