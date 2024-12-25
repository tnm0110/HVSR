function [amp_noise]=sta_lta_f(sac,sta,lta,r_min,r_max)

%   This function takes the sac file and run STA/LTA algorithom to get trasient spikes then run density based clustering 
%   algoritho (DBSCAN) to isolate trasinent noise cluster to remove from the waveform 
% 
%   Input :: sac   ==> day long sac file
%            sta   ==> stort-time average window (s)
% 	     lta   ==> long-time average windo (s)
% 	     r_min ==> minimum STA/LTA threshold  [ needed if you want to use STA/LTA only ]
%	     r_max ==> maximum STA/LTA threshold
%   Output :: amp_noise ==> noise window devoid of any transient to calculate HVSR in the next step 
% 
%   			Md Mohimanul Islam
%	    	   University of Missouri Columbia 
%	   	     Last Modified : 12/02/2023

%
%

%file='2020.092.10.00.00.X1.205.DPZ.SAC';
%sf=rdsac(file);

amp=sac.d;
dt=sac.HEADER.DELTA;
f=1/dt;

% STA LTA value in seconds 
%sta=0.5; 
%lta=60;

% STA LTA threshold range
%r_min = 0.5;
%r_max = 1.50;


sta_n=sta*f;
lta_n=lta*f;

% reamove mean and trend
amp_m=amp - mean(amp);
amp_tr=detrend(amp_m);

% remove the constant amplitdue time window
tc=1e-2;  % tolerance limit 
[amp_2]=remove_cons_amp_2(amp_tr,tc);



%% apply a bandpass butterworth filter
% % filter parameters
%disp('Filtering the waveform ...')
fc1 = 0.1;    % Lower cutoff frequency of the bandpass filter
fc2 = 10;     % Upper cutoff frequency of the bandpass filter
order = 3;    % Filter order
fn=f/2;
% Design the bandpass filter using Butterworth design
[b, a] = butter(order, [fc1/fn, fc2/fn], 'bandpass');

% Apply the bandpass filter using filtfilt
amp_filt = filtfilt(b, a, amp_2);

%% STA/LTA moving averagae
disp('running STA/LTA window ...')
amp_filt=amp_2;
% create a moving window to calculate STA/LTA
amp_sta=movmean(abs(amp_filt),sta_n);
amp_lta=movmean(abs(amp_filt),lta_n);

r=abs(amp_sta./amp_lta); % STA/LTA ratio

% mask the wavefrom by STA/LTA threshold range 

mask = (r >= r_min) & (r <= r_max);

% Extract portions within threshold range
amp_3 = amp_filt(mask);


%% find the clustering of the points and remove them as whole
disp('removing the trigggered clusters ...')
[rr, ~] = find(r > r_max);

eps=150;     % Max distance between points to assign as the same cluster
minPts=50;   % minumum point needed to be a cluster
bsize=2000;  % batach size for seqeuntial dbscan

%[cluster_sp]=get_cluster(rr,eps,minPts); % for total, tooks a lot of time
[cluster_sp]=get_cluster_2(rr,eps,minPts,bsize); % for sequential search

% % epand the trasinent cluster range by some factor
% ef=1.25; %expension factor
% exp_lim=ceil((cluster_sp(2:end-1,2)-cluster_sp(2:end-1,1))*ef);
% %cluster_sp(2:end-1,1)=cluster_sp(2:end-1,1)- exp_lim;
% cluster_sp(2:end-1,2)=cluster_sp(2:end-1,2) + exp_lim;

% remove the coherent trasinent from the waveform
amp_filt=amp_filt'; % if not using filtered signal 

% Create a logical mask
if size(cluster_sp,1) > 1
    msk2 = true(size(amp_filt));

    % Set mask elements to false for the specified index ranges
    for i = 1:size(cluster_sp, 1)
        msk2(cluster_sp(i, 1):cluster_sp(i, 2)) = false;
    end

    % remivng the transient applying the mask
    amp_4 = amp_filt(msk2);
else
    amp_4=amp_filt;
end

% do the STA/LTA again

%% create a new moving window to calculate STA/LTA
disp('running second STA/LTA ...')
amp_sta2=movmean(abs(amp_4),sta_n);
amp_lta2=movmean(abs(amp_4),lta_n);

r=abs(amp_sta2./amp_lta2); % STA/LTA ratio

% Define threshold range
r_min = 0.5;
r_max = 1.5;
mask = (r >= r_min) & (r <= r_max);

% Extract portions within threshold range
amp_5 = amp_4(mask);
[~, rr2] = find(r > r_max);
rr2=rr2';

%% second cluster
disp('removing second triggerd cluster ...')
tic
eps=150;
minPts=50;

%[cluster_sp2]=get_cluster(rr2,eps,minPts)
[cluster_sp2]=get_cluster_2(rr2,eps,minPts,bsize);

%  % epand the trasinent cluster range by some factor
% ef=1.25; %expansion factor
% exp_lim=ceil((cluster_sp2(2:end-1,2)-cluster_sp2(2:end-1,1))*ef);
% cluster_sp2(2:end-1,1)=cluster_sp2(2:end-1,1)- exp_lim;
% cluster_sp2(2:end-1,2)=cluster_sp2(2:end-1,2) + exp_lim;

% remove the coherent trasinent from the waveform
amp_filt2=amp_4;

% Create a logical mask
if size(cluster_sp2,1) > 1
    msk2 = true(size(amp_filt2));
    % Set mask elements to false for the specified index ranges
    for i = 1:size(cluster_sp2, 1)
        msk2(cluster_sp2(i, 1):cluster_sp2(i, 2)) = false;
    end

% Create the final noise waveform by applying the mask
amp_noise = amp_filt2(msk2);
else
    amp_noise = amp_filt2;
end
end
