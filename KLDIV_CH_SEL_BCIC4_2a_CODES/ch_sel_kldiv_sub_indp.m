%% Channel Selection Algorithm Using KL Divergence
% To be submitted in IEEE TBME
% Author - Raghav Dev @ IIT Delhi, rdevm23@gmail.com
% Copyright (C), Raghav Dev, IIT Delhi
% last modified on - 19th Oct, 2023
% Completely deterministic. Replicable 100%.

%% Parameters and filters


clear, clc

all_channels = 1:22;
n_subjects = 9;
trials = 1:48;
tmin = 3;
tmax = 6;
path = "C:\Raghav\git_repos\ch_selection_TBME\channel_selection-main\datasets\bcic4_2a\";
TorE = ["T","E"];
dataset = struct;
for s = 1:n_subjects
    for t_e = 1:numel(TorE)
        subj = strcat(path,'A0',num2str(s));
        % correcting the runs for A04T
        if (s == 4)&& (t_e == 1)
            runs = 2:7;
        else
            runs = 4:9;
        end 
        [dataset(s,t_e).imagery, ...
            dataset(s,t_e).classes, ...
            dataset(s,t_e).f_samp]  = ...
            extraction(strcat(subj,TorE(t_e)),...
            runs,all_channels,trials,tmin,tmax);
    end
end

x = dataset;
dataset = struct;
for s = 1:n_subjects
    dataset(s).imagery = [x(s,1).imagery;x(s,2).imagery];
    dataset(s).classes = [x(s,1).classes;x(s,2).classes];
    dataset(s).f_samp = x(s,1).f_samp;
end

%%
x = dataset;
for s = 1:n_subjects
dataset(s).imagery = zeros(6336,750);
dataset(s).classes = zeros(6336,1);
end

for s = 1:n_subjects
    idx = 1;
    for r = 1:size(x(s).imagery,1)
        if (x(s).classes(r) == 1)||(x(s).classes(r) == 2)
            dataset(s).imagery(idx,1:750) = x(s).imagery(r,:);
            dataset(s).classes(idx,1) = x(s).classes(r,1);
            idx = idx + 1;
        end
    end
    disp(s)
end


%%
n_channels = numel(all_channels);
filters_ON = 1;

dataset_filt = struct;
% Filtering
n_subs = 9;
fs = 250;
if filters_ON == 1
for s = 1:n_subjects
    
    signal_pre_proc = dataset(s).imagery;
    class = dataset(s).classes;
    
    classify_obj = cls_classify_motorI(...
        signal_pre_proc,class,n_channels,fs);
    
    % Creating the filter
    if (s==1)
        classify_obj = classify_obj.create_filter();
        filter = classify_obj.filter;
    else
        classify_obj.filter = filter;
    end
    
    % Other Preprocessing
    classify_obj = classify_obj.filter_bank();
    
    dataset_filt(s).eeg = classify_obj.signal;
    dataset_filt(s).label = class;
    disp('filtering is done for one subject')
end
else
    dataset_filt = dataset;
end
  

%% Finding the Ranks

all_klds = find_the_dist_(n_channels,dataset_filt);
%%
all_ranks = rankings(n_channels,all_klds);

%% User Defined Functions
% 1. all_ranks = rankings(n_channels,all_klds)
% 2. all_klds = find_the_dist_(n_channels,dataset)
% 3. norm_sigs = normlization(dataset)
% 4. [counts,centers] = find_the_hist(signal)
% 5. expect = expectation(dist)

function all_ranks = rankings(n_channels,all_klds)
all_ranks = struct;

all_ranks.refs(1).ranks = zeros(n_channels,1);
all_ranks.refs(2).ranks = zeros(n_channels,1);
all_ranks.refs(3).ranks = zeros(n_channels,1);
klds_c3 = all_klds.klds.c3;
klds_cz = all_klds.klds.cz;
klds_c4 = all_klds.klds.c4;

for c = 1:n_channels
    all_ranks.refs(1).ranks(c,1) = expectation(klds_c3(c,:));
    all_ranks.refs(2).ranks(c,1) = expectation(klds_cz(c,:));
    all_ranks.refs(3).ranks(c,1) = expectation(klds_c4(c,:));
    disp(c)
end

end



function all_klds = find_the_dist_(n_channels,dataset)
% dataset struct of size 1*n_subs with eeg and labeles as field
n_subs = size(dataset,2);
n_trials = int32((size(dataset(1).eeg,1)/n_channels));

n_trials_subs = int32((size(dataset(1).eeg,1)/n_channels)*n_subs);
n_samples = size(dataset(1).eeg,2);
norm_sigs = normlization(dataset);
all_klds = struct;


signals = zeros(n_channels,n_trials_subs,n_samples);

for i = 1: n_subs
    signals_of_one_sub = norm_sigs(i).eeg_logged;
    signals_of_one_sub = reshape(signals_of_one_sub,[n_channels,n_trials,n_samples]);
    idx_strt = (i-1)*n_trials + 1;
    idx_end = ((i)*n_trials );
    signals(:,idx_strt:idx_end,:) = signals_of_one_sub;
end

% signals size = n_channels*n_trials*n_samples
cz_sig = reshape(signals(8,:,:),[n_trials_subs,n_samples]);
c3_sig = reshape(signals(10,:,:),[n_trials_subs,n_samples]);
c4_sig = reshape(signals(12,:,:),[n_trials_subs,n_samples]);


all_klds.klds.c3 = zeros(n_channels,n_samples);
all_klds.klds.cz = zeros(n_channels,n_samples);
all_klds.klds.c4 = zeros(n_channels,n_samples);

for j = 1:n_samples
    % size c3, c4, cz  = n_trials*n_samples
    [c3_count,centers] = find_the_hist(c3_sig(:,j));
    [cz_count,~] = find_the_hist(cz_sig(:,j));
    [c4_count,~] = find_the_hist(c4_sig(:,j));
    for k = 1:n_channels
        ck_signals = reshape(signals(k,:,j),[n_trials_subs,1]);
        [ck_count,~] = find_the_hist(ck_signals);
        kld_c3 = kldiv(centers,ck_count,c3_count,'js');
        kld_cz = kldiv(centers,ck_count,cz_count,'js');
        kld_c4 = kldiv(centers,ck_count,c4_count,'js');
        
        
        all_klds.klds.c3(k,j) = kld_c3;
        all_klds.klds.cz(k,j) = kld_cz;
        all_klds.klds.c4(k,j) = kld_c4;
    end
    disp(j)
end

end




function norm_sigs = normlization(dataset)
% this witll log the data
n_subs = size(dataset,2);
norm_sigs = struct;
for i = 1: n_subs
    eeg = dataset(i).eeg;
    %make_it_one = (max(eeg,[],2)-min(eeg,[],2));
    
    eeg = (eeg - min(eeg,[],2))./(max(eeg,[],2)-min(eeg,[],2));
    eeg_logged = log(1 + eeg);
    eeg_logged = (eeg_logged-min(eeg_logged,[],2))...
        ./(max(eeg_logged,[],2)-min(eeg_logged,[],2));
    norm_sigs(i).eeg_logged = eeg_logged;
end


end

function [counts,centers] = find_the_hist(signal)
% signal = 1*n_trials
% histos = 1*n_trials

bin_edge = 0:0.1:1;
[counts,centers] = hist(signal,bin_edge);
counts = counts + 1; % 1 IS ADDED TO AVOID THE NaN AND inf
counts = counts/sum(counts);

end


function expect = expectation(dist)
% n_bins = 10;
[counts,centers] = hist(dist); % Erlier I had used hist() and it has no bins it was auto.
expect = sum(counts.*centers);
end




