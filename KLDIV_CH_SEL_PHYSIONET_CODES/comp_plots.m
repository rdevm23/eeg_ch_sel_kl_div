
%Comp_with_others_physionet.eps

% clear
% 
% 
% selected_channels_rand = zeros(118,10);
% for i = 1:10
%     rng(i*5647*i^2)
%     selected_channels_rand(:,i) = randperm(118,118);
% end
% 


clc
clear 
close all
figure(1)

linewidth = 2;

load('RESULTS_RAND_PHISIONET.mat')
x_all = zeros(62,10);
for i = 1:10
    
    hold on
    for j = 3:64
        x_all(j,i) = RESULTS_PHYSIONET(j,i).avg_results(1);
    end
    plot(x_all(3:64,i),'DisplayName',string(i),'LineWidth',linewidth)
end
legend
grid on

%%

sub_indp = load('RESULTS_PHYSIONET_SUB_INDP.mat');

rand_avg = mean(x_all,2);

rand_avg = rand_avg(3:64);

y_sub_indp = zeros(64,1);
for i = 3:64
    y_sub_indp(i) = sub_indp.RESULTS_PHYSIONET(i).avg_results(1);
end
y_sub_indp(1:2) = [];


plot(y_sub_indp,'k','DisplayName',"Sub Indp",'LineWidth',linewidth)
hold off
legend
grid on


fontzie = 12;
ylabel("Accuracy (%)",'FontSize',fontzie,'FontWeight','bold')
xlabel("Number of Selected Channels",'FontSize',fontzie,'FontWeight','bold')



figure(2)
plot(y_sub_indp,'k','DisplayName',"Sub Indp",'LineWidth',linewidth)
hold on
% plot(x,'DisplayName',"Rand Max",'LineWidth',linewidth)
plot(rand_avg,'DisplayName',"Rand Avg",'LineWidth',linewidth)
legend
grid on

fontzie = 12;
ylabel("Accuracy (%)",'FontSize',fontzie,'FontWeight','bold')
xlabel("Number of Selected Channels",'FontSize',fontzie,'FontWeight','bold')

%% Fishers

fishers = load('RESULTS_COMP_FISHERS_PHISIONET.mat');
fishers = fishers.RESULTS_PHYSIONET;


y_fishers = zeros(64,1);
for i = 3:64

    y_fishers(i) = fishers(i).avg_results(1);
    
end

y_fishers(1:2) = [];
plot(y_fishers,'DisplayName',"FS",'LineWidth',linewidth)



%% CSP 



csp_rank = load('RESULTS_COMP_CSP_RANK_PHYSIONET.mat');
csp_rank = csp_rank.RESULTS_PHYSIONET;



y_csp = zeros(64,1);
for i = 3:64

    y_csp(i) = csp_rank(i).avg_results(1);
    
end

y_csp(1:2) = [];
plot(y_csp,'DisplayName',"CSP Rank",'LineWidth',linewidth)

%% NMI

nmi_res = load('RESULTS_PHISIONET_NMIs_2.mat');

nmi_res = nmi_res.RESULTS_PHYSIONET;




y_nmi = zeros(64,1);
for i = 3:64

    y_nmi(i) = nmi_res(i,1).avg_results(1);
    
end

y_nmi(1:2) = [];
plot(y_nmi,'DisplayName',"NMI",'LineWidth',linewidth)



%%

n_chs = 3:64;
n_chs = n_chs';


ranks_acc_per_ch = y_sub_indp./n_chs;

five_per_max = (y_sub_indp(62))*(1- 1/20);



%% For the table



load('RESULTS_RAND_PHISIONET.mat')
rand_at_40 = zeros(10,3);
for i = 1:10
    
   rand_at_40(i,:) = RESULTS_PHYSIONET(40,i).avg_results;
end
rand_at_40_mean = mean(rand_at_40);


%%

aaaa = zeros(64,4);
% y_sub_indp = y_nmi;
for i = 1:62
    if y_sub_indp(i) > rand_avg(i)
        aaaa(i+2,1) = 1;
    end
    if y_sub_indp(i) > y_csp(i)
        aaaa(i+2,2) = 1;
    end
    if y_sub_indp(i) > y_fishers(i)
        aaaa(i+2,3) = 1;
    end
    if y_sub_indp(i) > y_nmi(i)
        aaaa(i+2,4) = 1;
    end
end
sum(aaaa)

%% 
aaaa = zeros(64,4);
% y_sub_indp = y_nmi;
for i = 1:62
        aaaa(i+2,1) = 1;
        aaaa(i+2,2) = 1;
        aaaa(i+2,3) = 1;
        aaaa(i+2,4) = 1;
end
sum(aaaa);
%% Standard Dev and res for the TABLE

% 

rand_all_res = load('RESULTS_RAND_PHISIONET.mat');
sub_indp = load('RESULTS_PHYSIONET_SUB_INDP.mat');
nmi_res = load('RESULTS_PHISIONET_NMIs_2.mat');
fishers = load('RESULTS_COMP_FISHERS_PHISIONET.mat');
csp_rank = load('RESULTS_COMP_CSP_RANK_PHYSIONET.mat');

all_results = struct;
all_results.avg = zeros(64,3,5);
all_results.sd = zeros(64,3,5);


rand_all_avg_sub = struct;

for c = 3:64
    acc = zeros(94,3);
    for r = 1:10
        acc = acc + rand_all_res.RESULTS_PHYSIONET(c,r).accuracy_results;
    end
    acc = acc/10;
    rand_all_avg_sub(c,1).acc = acc;
    
end

for c = 3:64
    all_results.avg(c,:,1) = mean(rand_all_avg_sub(c,1).acc);
    all_results.sd(c,:,1) = std(rand_all_avg_sub(c,1).acc);
    
    all_results.avg(c,:,5) = mean(sub_indp.RESULTS_PHYSIONET(c).accuracy_results);
    all_results.sd(c,:,5) = std(sub_indp.RESULTS_PHYSIONET(c).accuracy_results);
    
    all_results.avg(c,:,4) = mean(nmi_res.RESULTS_PHYSIONET(c).accuracy_results);
    all_results.sd(c,:,4) = std(nmi_res.RESULTS_PHYSIONET(c).accuracy_results);
    
    all_results.avg(c,:,3) = mean(fishers.RESULTS_PHYSIONET(c).accuracy_results);
    all_results.sd(c,:,3) = std(fishers.RESULTS_PHYSIONET(c).accuracy_results);
    
    all_results.avg(c,:,2) = mean(csp_rank.RESULTS_PHYSIONET(c).accuracy_results);
    all_results.sd(c,:,2) = std(csp_rank.RESULTS_PHYSIONET(c).accuracy_results);
    
end




% now the table: chs = 10,25,40,55


table_res = struct;

table_res.acc = zeros(5,12);
table_res.std = zeros(5,12);
chs = [10,25,40,55];

for ch = 1:numel(chs)
for m = 1:5
    idx = 3*(ch-1)+1;
    table_res.acc(m,idx:idx+2) = all_results.avg(chs(ch),:,m);
    table_res.std(m,idx:idx+2) = all_results.sd(chs(ch),:,m);
end
end


disp("done")



%% 


