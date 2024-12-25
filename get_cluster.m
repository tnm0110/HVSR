function [cluser_sp]=get_cluster(amp,eps,minPts)

% this function get the trigger point cluster and get the range 
% of the culster 
%

%load clr.mat    % test data
% epsilon = 150; % Max distance between points in the same cluster
% minPts = 50;   % Min number of points required to form a dense region
% amp=rr;

%data = [2, 3, 4, 4.5, 10, 11, 12, 12.5, 40, 41, 45 100 110 230];

% Fit DBSCAN clustering
%X = amp;
[idx, ~] = dbscan(amp, eps, minPts);
idx=[idx amp];

% Display cluster indices
% disp("Cluster Indices:");
% disp(idx);

% Display whether each point is noise or not
% disp("Is Noise:");
% disp(isNoise);
zz=idx(:,1);
ii=find(zz > -1);
idx=idx(ii,:);


uniqueValues = unique(idx(:, 1)); 
numUnique = numel(uniqueValues);
%matrix3D = zeros(size(idx, 1), size(idx, 2), numUnique);

for i = 1:numUnique
    kk = idx(:, 1) == uniqueValues(i);
    mat2=idx(kk,:);
    mat3(i,:)=[mat2(1,2) mat2(end,2)];
end
cluser_sp=mat3;
end