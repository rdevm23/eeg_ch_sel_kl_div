

clear 
close all
% clc

n_iter_optimization = 10000;
n_selected_chs = 20;
self = cls_channel_selection_old_code(1);
self = self.import_positions();
self = self.min_shunting(n_iter_optimization,n_selected_chs);

[min_loss, idx_min] = min(self.all_loss_shunting);
comb = nchoosek(1:n_selected_chs,2);
lossy_length = size(comb,1);

k_d_max = 1.790569827*lossy_length;
k_d_min = 0.003978427*lossy_length;


fprintf('min loss is %d...\n',min_loss);
fprintf('k_d_max %d...\n',k_d_max);
fprintf('k_d_min %d...\n',k_d_min);

chs_with_min_shunt = self.all_comb_of_chs(idx_min,:)';

[max_loss, idx_max] = max(self.all_loss_shunting);
chs_with_max_shunt = self.all_comb_of_chs(idx_max,:)';


% 
% figure(1)
% plot(nfo.xpos,nfo.ypos,'o')
% grid on
% axis([-2 2 -1 1])
% text(nfo.xpos,nfo.ypos,nfo.clab')
% text(nfo.xpos,nfo.ypos,strcat('\newline',string(1:118)))
% 
% 
% 








% 
% 
% N = 118; % Max integer value the data can have.
% k = 100; % Number of rows you want in the final output
% p = 20;  % Number of columns
% % Make a k-by-p array of values from 1 to N
% % Make more rows than we need because we may throw out some if they're duplicates.
% rows = 5 * k;
% data = randi(N, rows, p);
% 
% % Now we have sample data, and we can begin.....
% % First go down row-by-row throwing out any lines with duplicated numbers like [1, 1, 4]
% goodRows = []; % Keep track of good rows.
% % Bad rows will have unique() be less than the number of elements in the row.
% for row = 1 : rows
%     thisRow = data(row, :);
%     if length(unique(thisRow)) == length(thisRow)
%         goodRows = [goodRows, row];
%     end
% end
% if ~isempty(goodRows)
%     data = data(goodRows,:);
% end
% % Sort row-by-row with columns going in ascending order
% [~, sortIndexes] = sort(data, 2);
% % Get the list of rows all scrambled up
% [uniqueRows, ia, ~] = unique(data, 'rows', 'stable');
% % ia is the rows from A that we kept (extracted and stored in uniqueRows).
% % Extract those same rows from sortIndexes so we know how to "unsort" the rows
% sortIndexes2 = sortIndexes(ia,:);
% % Extract the first N rows.
% sortedOutput = uniqueRows(1:N,:);
% % Now "unsort" them if you want to do that.
% output = zeros(size(sortedOutput, 1),p);
% for row = 1 : size(sortedOutput, 1)
%     output(row,:) = sortedOutput(row, sortIndexes2(row, :));
% end