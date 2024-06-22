

clc
clear




class_1 = 1;
class_2 = 2;

n_subs = 94;


channel_selected = 1:64;
flag_precue = 0;

dataset = extraction_phisionet(channel_selected,flag_precue);



features = zeros(numel(channel_selected),n_subs*...
    numel(dataset(1).eeg)/numel(channel_selected));

samp_per_sub = numel(dataset(1).eeg)/numel(channel_selected);
labels = zeros(numel(channel_selected),size(features,2));
for i = 1:n_subs
    
    idx_strt = (i-1)*samp_per_sub + 1;
    idx_end = i*samp_per_sub;
    
    features(:,idx_strt:idx_end) = ...
        reshape(dataset(i).eeg, ...
        [numel(channel_selected),...
        samp_per_sub]);
    labels(:,idx_strt:idx_end) = ...
        reshape(repmat(dataset(i).label,1,640),...
        [numel(channel_selected), samp_per_sub]);
    
    disp(i)
end

labels = labels(1,:)';

%%

[index,featureScore] = feature_rank(features,labels);




