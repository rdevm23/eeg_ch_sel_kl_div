clear
clc
close all

%                C:\Raghav\git_repos\ch_selection_TBME\CHANNEL_SEL_INFORMATION_THEORY_PROJECT\KLDIV_CH_SEL_BCIC3_4a_CODES
sub_indp = load("C:\Raghav\git_repos\ch_selection_TBME\CHANNEL_SEL_INFORMATION_THEORY_PROJECT\MATS_FILES_RESULTS_FULLPAPER\RESULTS_ALL_SUBJECT_WITH_SUB_INDP.mat");
sub_dp   = load("C:\Raghav\git_repos\ch_selection_TBME\CHANNEL_SEL_INFORMATION_THEORY_PROJECT\MATS_FILES_RESULTS_FULLPAPER\RESULTS_ALL_SUBJECT_WITH_ALL_CHS_OF_ALL_SUBS.mat");

% col = ["m","g","c","y","b","r","k"];
line_width = [1,1.5,2];
DisplayName = ["Sub 1","Sub 2","Sub 3","Sub 4","Sub 5","Avg Ranks","Sub Indp"];

x_lims_low = [55,75,50,60,75,70];
x_lims_high = [100,100,85,105,105,95];
fontzie = 10;

signal = zeros(1,118);
avg_sub_dep = zeros(1,116);
avg_sub_dep_on_other_sub = zeros(10,116);
avg_of_avg_sub_dep = zeros(1,116);
for filt = 1:2
    figure(filt)
    avg_sub_dep = zeros(1,116);
    avg_of_avg_sub_dep = zeros(1,116);

    for s = 1:6 % 6th is for avg accuracy over subjects
        avg_sub_dep_on_other_sub = zeros(10,116);
        kk = 1;
        subplot(2,3,s)
        hold on
        for idx_sub_dp = 1:7 % 6th if for avg rank
            for ch = 3:118
                if idx_sub_dp ~= 7
                    % disp("Sub Dep");
                    if s ~= 6
                        signal(ch) = sub_dp.RESULTS_BCIC3_4A(ch,filt,idx_sub_dp).accuracy_results(s,1);
                    else
                        signal(ch) = sub_dp.RESULTS_BCIC3_4A(ch,filt,idx_sub_dp).avg_results(1);
                    end
                else
                    if s ~= 6
                        signal(ch) = sub_indp.RESULTS_BCIC3_4A(ch,filt).accuracy_results(s,1);
                    else
                        signal(ch) = sub_indp.RESULTS_BCIC3_4A(ch,filt).avg_results(1);
                    end
                    % disp("Sub Indep");
                end
            end
            signal(1:2) = [];
            if (idx_sub_dp == s) && (idx_sub_dp <= 5)
                
                LineWidth = line_width(2);
                plot(signal(1:1:116),'LineWidth',...
                    LineWidth,'DisplayName',DisplayName(idx_sub_dp))
                
                avg_sub_dep = avg_sub_dep + 0.2*signal;
                
            elseif idx_sub_dp == 7
                LineWidth = line_width(3);
                plot(signal(1:1:116),"k",'LineWidth',...
                    LineWidth,'DisplayName',DisplayName(idx_sub_dp))
                
            elseif idx_sub_dp == 6
                LineWidth = line_width(2);
                plot(signal(1:1:116),"r",'LineWidth',...
                    LineWidth,'DisplayName',DisplayName(idx_sub_dp))
                
            else
                % disp(idx_sub_dp)
                % disp(s)
                avg_sub_dep_on_other_sub(idx_sub_dp,:) = signal;
            end
            
            
            
            xlim([2 120]); ylim([x_lims_low(s) x_lims_high(s)]); grid on;
            ylabel("Accuracy (%)",'FontSize',fontzie,'FontWeight','bold')
            xlabel("Number of Channels",'FontSize',fontzie,'FontWeight','bold')
            if s <= 5
                title(strcat("Subject ", string(s)));
            else
                title("Average Accuracy Over Subjects");
            end
            
        end
        
        if s ~= 6
            avg_sub_dep_on_other_sub_ = sum(avg_sub_dep_on_other_sub,1)/4;
            LineWidth = line_width(2);
            plot(avg_sub_dep_on_other_sub_(1:1:116),...
                "Color","#77AC30",'LineWidth',...
                LineWidth,'DisplayName',"Avg Sub Dep ~X")
            avg_of_avg_sub_dep = ...
                avg_of_avg_sub_dep...
                + 0.2*avg_sub_dep_on_other_sub_;
            
        elseif s == 6
            
            LineWidth = line_width(2);
            plot(avg_of_avg_sub_dep(1:1:116),...
                "Color","#77AC30",'LineWidth',...
                LineWidth,'DisplayName',"Avg Sub Dep ~X")
            
            
            LineWidth = line_width(2);
            plot(avg_sub_dep(1:1:116),"Color","#0072BD",'LineWidth',...
                LineWidth,'DisplayName',"Sub Dep")
        end
        hold off
    end
    legend;
end
%%

% avg_sub_dep_on_other_sub_ = sum(avg_sub_dep_on_other_sub,1)/4;
% 
% for i = 1:2
%     for j = 1:6
%         
%     end
% end




%% Filter Vs Without Filter


figure(3)


col = ["#D95319","#000000"];
line_width = [1,2.5];
DisplayName1 = ["Without Filter","Filtered"];
DisplayName2 = [" Avg Rank"," Subject Ind"," Sub 1"," Sub 2"];

x_lims_low = 65;
x_lims_high = 95;
s = 6;
signal_avg_rank = zeros(1,118);
signal_sub_indp = zeros(1,118);
signal_sub1 = zeros(1,118);
signal_sub2 = zeros(1,118);
for filt = 1:2
    for ch = 3:118
        signal_avg_rank(ch) = sub_dp.RESULTS_BCIC3_4A(ch,filt,6).avg_results(1);
        signal_sub_indp(ch) = sub_indp.RESULTS_BCIC3_4A(ch,filt).avg_results(1);
        signal_sub1(ch) = sub_dp.RESULTS_BCIC3_4A(ch,filt,1).avg_results(1);
        signal_sub2(ch) = sub_dp.RESULTS_BCIC3_4A(ch,filt,2).avg_results(1);

    end    
    signal_avg_rank(1:2) = [];
    signal_sub_indp(1:2) = [];
    signal_sub1(1:2) = [];
    signal_sub2(1:2) = [];
    
    subplot(2,2,1)
    hold on
    plot(signal_avg_rank,'MarkerFace',"#D95319",'LineWidth',line_width(filt),...
        'DisplayName',strcat(DisplayName1(filt),DisplayName2(1)))
    hold off
    
    xlim([2 120]); 
    ylim([x_lims_low x_lims_high]); 
    grid on;
    
    ylabel("Accuracy (%)",'FontSize',fontzie,'FontWeight','bold')
    xlabel("Number of Channels",'FontSize',fontzie,'FontWeight','bold')
    
    legend;

    subplot(2,2,2)
    hold on
    plot(signal_sub_indp,'MarkerFace',"#000000",'LineWidth',line_width(filt),...
        'DisplayName',strcat(DisplayName1(filt),DisplayName2(2)))
    
    hold off
    
    xlim([2 120]); 
    ylim([x_lims_low x_lims_high]); 
    grid on;
    ylabel("Accuracy (%)",'FontSize',fontzie,'FontWeight','bold')
    xlabel("Number of Channels",'FontSize',fontzie,'FontWeight','bold')
    legend;
    
    subplot(2,2,3)
    hold on
    plot(signal_sub1,'MarkerFace',"#4DBEEE",'LineWidth',line_width(filt),...
        'DisplayName',strcat(DisplayName1(filt),DisplayName2(3)))
    hold off
    
    xlim([2 120]); 
    ylim([x_lims_low x_lims_high]); 
    grid on;
    
    ylabel("Accuracy (%)",'FontSize',fontzie,'FontWeight','bold')
    xlabel("Number of Channels",'FontSize',fontzie,'FontWeight','bold')
    
    legend;

    subplot(2,2,4)
    hold on
    plot(signal_sub2,'MarkerFace',"#77AC30",'LineWidth',line_width(filt),...
        'DisplayName',strcat(DisplayName1(filt),DisplayName2(4)))
    
    hold off
    
    xlim([2 120]); 
    ylim([x_lims_low x_lims_high]); 
    grid on;
    
    ylabel("Accuracy (%)",'FontSize',fontzie,'FontWeight','bold')
    xlabel("Number of Channels",'FontSize',fontzie,'FontWeight','bold')
    
    sgtitle("Filter Vs Without Filter for Average Accuracy Over Subjects");
    legend;

end



















