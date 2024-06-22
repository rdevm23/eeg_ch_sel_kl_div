
clc
clear 
close all


%% Standard Dev and res for the TABLE

% 

rand_all_res = load('RESULTS_RAND_PHISIONET.mat');
sub_indp = load('RESULTS_PHYSIONET_SUB_INDP.mat');
nmi_res = load('RESULTS_PHISIONET_NMIs_2.mat');
fishers = load('RESULTS_COMP_FISHERS_PHISIONET.mat');
csp_rank = load('RESULTS_COMP_CSP_RANK_PHYSIONET.mat');

rand_all_avg_sub = struct;

for c = 3:64
    acc = zeros(94,3);
    for r = 1:10
        acc = acc + rand_all_res.RESULTS_PHYSIONET(c,r).accuracy_results;
    end
    acc = acc/10;
    rand_all_avg_sub(c,1).acc = acc;
    
end


chs = [10, 25, 40, 55, 64];

%%
all_results = struct;
for i = 1:5
    all_results(i).acc = zeros(94,5);
end
%%

for i = 1:5
    
    all_results(i).acc(:,1) = ...
        rand_all_avg_sub(chs(i),1).acc(:,1);
    all_results(i).acc(:,2) = ...
        sub_indp.RESULTS_PHYSIONET(chs(i)).accuracy_results(:,1);
    all_results(i).acc(:,3) = ...
        nmi_res.RESULTS_PHYSIONET(chs(i)).accuracy_results(:,1);
    all_results(i).acc(:,4) = ...
        fishers.RESULTS_PHYSIONET(chs(i)).accuracy_results(:,1);
    all_results(i).acc(:,5) = ...
        csp_rank.RESULTS_PHYSIONET(chs(i)).accuracy_results(:,1);
    
end
%%


for i = 1:5
    for j = 1:5
        all_results(i).acc(:,j) = sort(all_results(i).acc(:,j),'descend');
    end
end


%%

avg_n_subs = [35, 50, 85, 94];


all_results_avg = struct;
for i = 1:5
    all_results_avg(i).acc = zeros(4,5);
    all_results_avg(i).sd = zeros(4,5);
end


for i = 1:5 % channels to be reported
    for j = 1:5 % number of methods
        for k = 1:4 % number of averages subjects reported
            subs = avg_n_subs(k);
            all_results_avg(i).acc(k,j) = ...
                mean(all_results(i).acc(1:subs,j));
            all_results_avg(i).sd(k,j) = ...
                std(all_results(i).acc(1:subs,j));
        end
    end
end




