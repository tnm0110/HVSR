
function [amp_filt]=remove_cons_amp(amp,tc)

% fuction to remove the constant time series signal


%load aa.mat
% tc=1e-2;
% amp=amp_tr;



% Find segments with constant amplitude
diff_amp = abs(diff(amp));

% Find indices where amplitude changes beyond the tolerance
ind_nc = find(diff_amp > tc);


% Initialize variables to store non-constant segments
filtered_time = [];
filtered_data = [];


% Extract non-constant segments

prev_idx = 1;
for i = 1:length(ind_nc)
    segment_idx = ind_nc(i);
    %filtered_time = [filtered_time, time(prev_idx:segment_idx)];
    filtered_data = [filtered_data, amp(prev_idx:segment_idx)'];
    prev_idx = segment_idx + 1;
end

amp_filt=filtered_data';

% figure(1)
% subplot(2,1,1)
% plot(amp_tr)
% xlim([0 length(amp)])
% subplot(2,1,2)
% plot(filtered_data)
% xlim([0 length(amp)])
end

