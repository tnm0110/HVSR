function [cluster_sp]=get_cluster_2(rr,eps,minPts,bsize)

% this function get the trigger point cluster and get the range 
% of the culster using DBSACN algorithom. It splits the data into chunks to inrease the speed
% Input :: rr   ==> wavefoms
%	   eps  ==> Max distance between points to be cosidered as  the same cluster (in data points, e.g. 150)
%        minPts ==> Min number of points required to form a dense regions to remove (e.g. 150)
%         bsize ==> batach size for seqeuntial dbscan (e.g. 20000)
%


amp=rr;

% Fit DBSCAN clustering
%X = amp;
% disp('inside cluster ==> running dbscan ...')
% tic
% [clabel, ~] = dbscan(amp, eps, minPts);
% clabel_amp=[clabel amp];  % cluster label with amphlitude
% %tst=idx;
% 
% toc

%disp('inside cluster ==> getting the cluster value in sequential ..')

% sequential dbscan

n_batcth = ceil(length(amp) / bsize);

Clabel_amp=zeros(length(amp),2); % final clutering level with amp

last_label=0;
for batch_idx=1:n_batcth
    
    st_idx = (batch_idx - 1) * bsize + 1;
    en_idx = min(st_idx + bsize - 1, length(amp));
    batchData = amp(st_idx:en_idx); % current batch data
    % run the DBSACN 
    [batchLabels, ~] = dbscan(batchData, eps, minPts);   
    blabel_amp=[batchLabels batchData];

    % Find indices where the first column is not -1
    add_ind = blabel_amp(:, 1) ~= -1;
    blabel_amp(add_ind, 1) = blabel_amp(add_ind, 1) + last_label;
    last_label=max(blabel_amp(:,1));

    Clabel_amp(st_idx:en_idx,:)=blabel_amp;

end


zz=Clabel_amp(:,1);
ii=zz > -1;
Clabel_cc=Clabel_amp(ii,:);   % levels that are only in the cluster
uniqueValues = unique(Clabel_cc(:, 1)); 
numUnique = numel(uniqueValues);

cluster_sp=zeros(length(numUnique),2);
for i = 1:numUnique
    kk = Clabel_cc(:, 1) == uniqueValues(i);
    mat2=Clabel_cc(kk,:);
    cluster_sp(i,:)=[mat2(1,2) mat2(end,2)];
end

%cluser_sp=mat3;
end
