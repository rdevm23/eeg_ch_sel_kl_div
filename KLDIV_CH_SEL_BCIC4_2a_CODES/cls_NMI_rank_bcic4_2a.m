classdef cls_NMI_rank_bcic4_2a
    
    properties (Access = public)
        % dataset_name: It is a flag will have 3 values
        %               a) 1 == for bcic3_4a
        %               b) 2 == for bcic4_2a
        %               c) 3 == for physionet
        dataset_filt
        electrode_positions % a struct of array; set in self = import_positions(self)
        all_comb_of_chs % shall be set in min_shunting()
        all_loss_shunting % shall be set in min_dispshunting()
        
        dist_ranks
        all_NMIs
        combs
    end
    
    %% ===================== PUBLIC METHOD =========================
    methods (Access = public)
        
        function self = cls_NMI_rank_bcic4_2a(dataset_filt)
            self.dataset_filt = dataset_filt;
        end
        
        function self = nmi(self)
            
            % no of total channels, will act as max int for the combinations
            n_channels = 22;
            n_subs = 9;
            n_trials = 288;
            n_samples = 750;
            
            all_ch_combs = nchoosek(1:n_channels,2);
            n_combs = size(all_ch_combs,1);
                
            % norm_dataset = normlization(self,self.dataset);
            NMIs = struct;
            for s = 1:n_subs
                % continue
                signal = self.dataset_filt(s).eeg;
                signal = reshape(signal,[n_channels,n_trials,n_samples]);
                
                for t = 1:n_trials
                    NMIs(s,t).all_combs = zeros(n_combs,1);
                    all_ch_entropy = zeros(n_channels,1);
                    for ch = 1: n_channels
                        ch_sig_per_ch = reshape(signal(ch,t,:),[1,n_samples]);
                        all_ch_entropy(ch) = ENTROPY(self,ch_sig_per_ch);
                    end
                    
                    
                    for i = 1:n_combs
                        
                        ch_x = all_ch_combs(i,1);
                        ch_y = all_ch_combs(i,2);
                        ch_sig_per_ch_x = reshape(signal(ch_x,t,:),[1,n_samples]);
                        ch_sig_per_ch_y = reshape(signal(ch_y,t,:),[1,n_samples]);
                        
                        joint_entropy = JOINT_ENTROPY(self,...
                            ch_sig_per_ch_x,ch_sig_per_ch_y);
                        shannon_entropy_x = all_ch_entropy(ch_x);
                        shannon_entropy_y = all_ch_entropy(ch_y);
                        
                        NMIs(s,t).all_combs(i) = NMI(self,joint_entropy,...
                            shannon_entropy_x,shannon_entropy_y);
                        
                    end
                    % disp(t)
                end
                disp(s)
            end
            self.all_NMIs = NMIs;
            self.combs = all_ch_combs;
        end
    end
    % ===============================================================
    
    %% ===================== PRIVATE METHOD =========================
    
    
    
    
    methods (Access = private)
        % Functions
        % 1. shannon_entropy = ENTROPY(~,ch_sig_per_ch)
        % 2. joint_entropy = JOINT_ENTROPY(~,...
        %       ch_sig_per_ch_x,ch_sig_per_ch_y)
        % 3. nmi = NMI(~,joint_entropy,...
        %       shannon_entropy_x,shannon_entropy_y)
        % 4. [d_max, d_min] = d_max_d_min(~,positions)
        
        
        
        function norm_sigs = normlization(~,dataset)
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
        
        function shannon_entropy = ENTROPY(~,ch_sig_per_ch)
            n_bins = 10;
            
            [N,~] = histcounts(ch_sig_per_ch,n_bins);
            % distribution.counts = N;
            % distribution.center = edges;
            
            shannon_entropy = sum((-1*N).*log(1 + N));
        end
        
        
        
        
        function joint_entropy = JOINT_ENTROPY(~,...
                ch_sig_per_ch_x,ch_sig_per_ch_y)
            n_bin = 10;
            
            
            [N,~,~] = histcounts2(ch_sig_per_ch_x,ch_sig_per_ch_y,n_bin);
            
            N = N(:);
            joint_entropy = sum((-1*N).*log(1 + N));
            
        end
        
        function nmi = NMI(~,joint_entropy,...
                shannon_entropy_x,shannon_entropy_y)
            mi = shannon_entropy_x + shannon_entropy_y - joint_entropy;
            
            nmi = mi/(shannon_entropy_x + shannon_entropy_y);
        end
        
        
        
        
        function [d_max, d_min] = d_max_d_min(~,positions)
            % positions: all the position of all the channels
            n_chs = size(positions,1);
            comb = nchoosek(1:n_chs,2); % all two combination of chs
            vect_i = positions(comb(:,1),:);
            vect_j = positions(comb(:,2),:);
            diff = 0.5*(sum((vect_i-vect_j).^2,2));
            d_max = max(diff);
            d_min = min(diff);
        end
        
        
    end
end















