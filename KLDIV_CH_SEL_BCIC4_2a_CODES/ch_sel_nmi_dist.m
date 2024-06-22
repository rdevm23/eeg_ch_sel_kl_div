clear
close all
% clc


%% Uncomment to regen ch_sel_obj.mat
% load('dataset_filtered.mat')
% ch_sel_obj = cls_NMI_rank_bcic4_2a(dataset_filt);
% ch_sel_obj = ch_sel_obj.nmi();
% NMIs_BCIC4_2a = ch_sel_obj.all_NMIs;
% save('NMIs_BCIC4_2a.mat','NMIs_BCIC4_2a');
%% 

load("NMIs_BCIC4_2a_STRUCT.mat");
n_channels = 22;
n_sub = 9;
n_trial = 288;
all_ch_combs = nchoosek(1:n_channels,2);
n_combs = size(all_ch_combs,1);
NMI_Array = zeros(n_sub,n_combs,n_trial);

for s = 1:n_sub
    for t = 1:n_trial
        NMI_Array(s,:,t) = NMIs_BCIC4_2a(s,t).all_combs;
    end
end

%%

NMI_avg = mean(mean(NMI_Array,3),1);
NMI_avg= NMI_avg';

alpha_all = 0; % Need to be tuned



alpha = alpha_all;
ordered_combs = ORDERED_COMBS(alpha,NMI_avg,all_ch_combs);


%% User Defined Function

function ordered_combs = ORDERED_COMBS(alpha,NMI_avg,all_ch_combs)

dist_and_nmi_ranks = (1-alpha)*NMI_avg;

[~,idx] = sort(dist_and_nmi_ranks,'descend');
ordered = all_ch_combs(idx,:);
ordered_combs = zeros(numel(all_ch_combs),1);
count = 1;
for i = 1:size(all_ch_combs,1)
    ordered_combs(count) = ordered(i,1);
    ordered_combs(count+1) = ordered(i,2);
    count = count + 2;
end


[~,b] = unique(ordered_combs,'first');
ordered_combs = ordered_combs(sort(b));


end




function y = normalization(x)
min_x = min(x);
max_x = max(x);
y = (x-min_x)/(max_x-min_x);
end











