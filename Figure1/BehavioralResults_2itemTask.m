%% behavior tally: Ocean/Lemon, correct rate by rank, by chosen item, by chosen rank
savedir = 'G:\PaperPrep';

Ocean_tags = {'0709','0710','0713','0714','0715','0720','0721','0722','0723','0726','0727','0729','0730'};
Ocean_behav_dir = 'G:\FW_data\NewJuly\ocean_data\preprocessed\trial_cuts';

Lemon_tags = {'0224','0227','0302','0304','0305','0306','0310','0311','0312','0313','0315','0316','0317'};
Lemon_behav_dir = 'G:\ZC_data\Preprocessed_New\trial_info';

monkey_tags = {'Ocean','Lemon'};
num_monkeys = length(monkey_tags);

num_ranks = 2;
num_items = 6;
num_days = 13;

rate_by_rank = zeros(num_monkeys,num_ranks,3); % dim 3: repeat/mirror/all
rate_by_chosen_item = zeros(num_monkeys,num_items,num_items,3); % dim 4: repeat/mirror/all
rate_by_chosen_rank = zeros(num_monkeys,num_ranks,num_ranks,3); % dim 4: repeat/mirror/all
dim_labels_rank = {'monkey','rank','repeat/mirror/all'};
dim_labels_chosen = {'monkey','rank/item_actual','rank/item_chosen','repeat/mirror/all'};

rate_by_rank_sepday = zeros(num_monkeys,num_ranks,3,num_days); % dim 3: repeat/mirror/all
rate_by_chosen_item_sepday = zeros(num_monkeys,num_items,num_items,3,num_days); % dim 4: repeat/mirror/all
rate_by_chosen_rank_sepday = zeros(num_monkeys,num_ranks,num_ranks,3,num_days); % dim 4: repeat/mirror/all
dim_labels_rank_sepday = {'monkey','rank','repeat/mirror/all','days'};
dim_labels_chosen_sepday = {'monkey','rank/item_actual','rank/item_chosen','repeat/mirror/all','days'};

for iM = 1:num_monkeys
    curr_monkey = monkey_tags{iM};
    
    eval(sprintf('num_files = length(%s_tags);',curr_monkey));
    eval(sprintf('tag_set = %s_tags;',curr_monkey));
    eval(sprintf('curr_dir = %s_behav_dir;',curr_monkey));
    for iFile = 1:num_files
        curr_tag = tag_set{iFile};
        monkey_str = lower(curr_monkey);
        load(sprintf('%s\\%s2021-%s-%s_TrialInfo_split.mat',...
            curr_dir,monkey_str,curr_tag(1:2),curr_tag(3:4)),...
            'Rule','TrueCorrect','Target','Response','CueFrameOnTime')
        
        usable_idx = find(~isnan(Rule) & ~isnan(CueFrameOnTime));
        
        % rate by rank
        repeat_R1_idx = intersect(intersect(usable_idx,find(Rule == 1)),find(~isnan(Response(:,1))));
        repeat_R2_idx = intersect(intersect(usable_idx,find(Rule == 1)),find(~isnan(Response(:,2))));
        mirror_R1_idx = intersect(intersect(usable_idx,find(Rule == 2)),find(~isnan(Response(:,2))));
        mirror_R2_idx = intersect(intersect(usable_idx,find(Rule == 2)),find(~isnan(Response(:,1))));
        T1_R_idx = intersect(repeat_R1_idx,find(Target(:,1) == Response(:,1)));
        T2_R_idx = intersect(repeat_R2_idx,find(Target(:,2) == Response(:,2)));
        T1_M_idx = intersect(mirror_R1_idx,find(Target(:,1) == Response(:,2)));
        T2_M_idx = intersect(mirror_R2_idx,find(Target(:,2) == Response(:,1)));
        rate_by_rank_sepday(iM,1,1,iFile) = length(T1_R_idx)/length(repeat_R1_idx);
        rate_by_rank_sepday(iM,2,1,iFile) = length(T2_R_idx)/length(repeat_R2_idx);
        rate_by_rank_sepday(iM,1,2,iFile) = length(T1_M_idx)/length(mirror_R1_idx);
        rate_by_rank_sepday(iM,2,2,iFile) = length(T2_M_idx)/length(mirror_R2_idx);
        rate_by_rank_sepday(iM,1,3,iFile) = (length(T1_R_idx)+length(T1_M_idx))/(length(repeat_R1_idx) + length(mirror_R1_idx));
        rate_by_rank_sepday(iM,2,3,iFile) = (length(T2_R_idx)+length(T2_M_idx))/(length(repeat_R2_idx) + length(mirror_R2_idx));
        
        % rate by chosen item
        for iT = 1:num_items
            repeat_R1_chosen_vec = zeros(1,num_items);
            repeat_R2_chosen_vec = zeros(1,num_items);
            mirror_R1_chosen_vec = zeros(1,num_items);
            mirror_R2_chosen_vec = zeros(1,num_items);
            
            for jT = 1:num_items
                repeat_R1_chosen_vec(jT) = length(intersect(repeat_R1_idx,find(Target(:,1)==iT & Response(:,1)==jT)));
                repeat_R2_chosen_vec(jT) = length(intersect(repeat_R2_idx,find(Target(:,2)==iT & Response(:,2)==jT)));
                mirror_R1_chosen_vec(jT) = length(intersect(mirror_R1_idx,find(Target(:,1)==iT & Response(:,2)==jT)));
                mirror_R2_chosen_vec(jT) = length(intersect(mirror_R2_idx,find(Target(:,2)==iT & Response(:,1)==jT)));
            end
            temp_repeat = repeat_R1_chosen_vec + repeat_R2_chosen_vec;
            temp_mirror = mirror_R1_chosen_vec + mirror_R2_chosen_vec;
            rate_by_chosen_item_sepday(iM,iT,:,1,iFile) = temp_repeat./sum(temp_repeat);
            rate_by_chosen_item_sepday(iM,iT,:,2,iFile) = temp_mirror./sum(temp_mirror);
            rate_by_chosen_item_sepday(iM,iT,:,3,iFile) = (temp_repeat+temp_mirror)./(sum(temp_repeat) + sum(temp_mirror));
        end
        
        % rate by chosen rank
        repeat_R1_chosen_T1 = length(intersect(repeat_R1_idx,find(Target(:,1) == Response(:,1))));
        repeat_R1_chosen_T2 = length(intersect(repeat_R1_idx,find(Target(:,1) == Response(:,2))));
        repeat_R2_chosen_T1 = length(intersect(repeat_R2_idx,find(Target(:,2) == Response(:,1))));
        repeat_R2_chosen_T2 = length(intersect(repeat_R2_idx,find(Target(:,2) == Response(:,2))));
        
        mirror_R1_chosen_T1 = length(intersect(mirror_R1_idx,find(Target(:,1) == Response(:,2))));
        mirror_R1_chosen_T2 = length(intersect(mirror_R1_idx,find(Target(:,2) == Response(:,2))));
        mirror_R2_chosen_T1 = length(intersect(mirror_R2_idx,find(Target(:,2) == Response(:,1))));
        mirror_R2_chosen_T2 = length(intersect(mirror_R2_idx,find(Target(:,1) == Response(:,1))));
        
        rate_by_chosen_rank_sepday(iM,1,1,1,iFile) = repeat_R1_chosen_T1/(repeat_R1_chosen_T1+repeat_R1_chosen_T2);
        rate_by_chosen_rank_sepday(iM,1,2,1,iFile) = repeat_R1_chosen_T2/(repeat_R1_chosen_T1+repeat_R1_chosen_T2);
        rate_by_chosen_rank_sepday(iM,2,1,1,iFile) = repeat_R2_chosen_T1/(repeat_R2_chosen_T1+repeat_R2_chosen_T2);
        rate_by_chosen_rank_sepday(iM,2,2,1,iFile) = repeat_R2_chosen_T2/(repeat_R2_chosen_T1+repeat_R2_chosen_T2);
        
        rate_by_chosen_rank_sepday(iM,1,1,2,iFile) = mirror_R2_chosen_T1/(mirror_R2_chosen_T1+mirror_R2_chosen_T2);
        rate_by_chosen_rank_sepday(iM,1,2,2,iFile) = mirror_R2_chosen_T2/(mirror_R2_chosen_T1+mirror_R2_chosen_T2);
        rate_by_chosen_rank_sepday(iM,2,1,2,iFile) = mirror_R1_chosen_T2/(mirror_R1_chosen_T1+mirror_R1_chosen_T2);
        rate_by_chosen_rank_sepday(iM,2,2,2,iFile) = mirror_R1_chosen_T1/(mirror_R1_chosen_T1+mirror_R1_chosen_T2);
        
        rate_by_chosen_rank_sepday(iM,1,1,3,iFile) = (repeat_R1_chosen_T1+mirror_R2_chosen_T1)/(repeat_R1_chosen_T1+mirror_R2_chosen_T1+repeat_R1_chosen_T2+mirror_R2_chosen_T2);
        rate_by_chosen_rank_sepday(iM,1,2,3,iFile) = (repeat_R1_chosen_T2+mirror_R2_chosen_T2)/(repeat_R1_chosen_T1+mirror_R2_chosen_T1+repeat_R1_chosen_T2+mirror_R2_chosen_T2);
        rate_by_chosen_rank_sepday(iM,2,1,3,iFile) = (repeat_R2_chosen_T1+mirror_R1_chosen_T2)/(repeat_R2_chosen_T2+mirror_R1_chosen_T1+repeat_R2_chosen_T1+mirror_R1_chosen_T2);
        rate_by_chosen_rank_sepday(iM,2,2,3,iFile) = (repeat_R2_chosen_T2+mirror_R1_chosen_T1)/(repeat_R2_chosen_T2+mirror_R1_chosen_T1+repeat_R2_chosen_T1+mirror_R1_chosen_T2);
    end
    
    % rate by rank
    for i = 1:num_ranks
        for j = 1:3
            rate_by_rank(iM,i,j) = mean(squeeze(rate_by_rank_sepday(iM,i,j,:)));
            rate_by_rank(iM,i,j) = mean(squeeze(rate_by_rank_sepday(iM,i,j,:)));
        end
    end
    
    % rate by chosen item
    for i = 1:num_items
        for j = 1:num_items
            for k = 1:3
                rate_by_chosen_item(iM,i,j,k) = mean(squeeze(rate_by_chosen_item_sepday(iM,i,j,k,:)));
            end
        end
    end
    
    % rate by chosen rank
    for i = 1:num_ranks
        for j = 1:num_ranks
            for k = 1:3
                rate_by_chosen_rank(iM,i,j,k) = mean(squeeze(rate_by_chosen_rank_sepday(iM,i,j,k,:)));
            end
        end
    end
    
end
save(sprintf('%s\\TwoMonk_Behav_2Seq.mat',savedir),'rate_*','dim_labels*')

%% Groot behavior tally: correct rate by rank, by chosen item, by chosen rank
savedir = 'G:\PaperPrep';

Groot_tags = {'0415','0507','0509','0513','0514','0515','0516','0520','0521','0522'};
Groot_behav_dir = 'E:\ZC_data\groot_data_switch\trial_cuts';

monkey_tags = {'Groot'};
num_monkeys = length(monkey_tags);

num_ranks = 2;
num_items = 6;
num_days = 10;

rate_by_rank = zeros(num_monkeys,num_ranks,3); % dim 3: repeat/mirror/all
rate_by_chosen_item = zeros(num_monkeys,num_items,num_items,3); % dim 4: repeat/mirror/all
rate_by_chosen_rank = zeros(num_monkeys,num_ranks,num_ranks,3); % dim 4: repeat/mirror/all
dim_labels_rank = {'monkey','rank','repeat/mirror/all'};
dim_labels_chosen = {'monkey','rank/item_actual','rank/item_chosen','repeat/mirror/all'};

rate_by_rank_sepday = zeros(num_monkeys,num_ranks,3,num_days); % dim 3: repeat/mirror/all
rate_by_chosen_item_sepday = zeros(num_monkeys,num_items,num_items,3,num_days); % dim 4: repeat/mirror/all
rate_by_chosen_rank_sepday = zeros(num_monkeys,num_ranks,num_ranks,3,num_days); % dim 4: repeat/mirror/all
dim_labels_rank_sepday = {'monkey','rank','repeat/mirror/all','days'};
dim_labels_chosen_sepday = {'monkey','rank/item_actual','rank/item_chosen','repeat/mirror/all','days'};

for iM = 1:num_monkeys
    curr_monkey = monkey_tags{iM};
    
    eval(sprintf('num_files = length(%s_tags);',curr_monkey));
    eval(sprintf('tag_set = %s_tags;',curr_monkey));
    eval(sprintf('curr_dir = %s_behav_dir;',curr_monkey));
    for iFile = 1:num_files
        curr_tag = tag_set{iFile};
        monkey_str = lower(curr_monkey);
        load(sprintf('%s\\%s2024-%s-%s_TrialInfo_split.mat',...
            curr_dir,monkey_str,curr_tag(1:2),curr_tag(3:4)),...
            'Rule','TrueCorrect','Target','Response','CueFrameOnTime')
        
        usable_idx = find(~isnan(Rule) & ~isnan(CueFrameOnTime));
        
        % rate by rank
        repeat_R1_idx = intersect(intersect(usable_idx,find(Rule == 1)),find(~isnan(Response(:,1))));
        repeat_R2_idx = intersect(intersect(usable_idx,find(Rule == 1)),find(~isnan(Response(:,2))));
        mirror_R1_idx = intersect(intersect(usable_idx,find(Rule == 2)),find(~isnan(Response(:,2))));
        mirror_R2_idx = intersect(intersect(usable_idx,find(Rule == 2)),find(~isnan(Response(:,1))));
        T1_R_idx = intersect(repeat_R1_idx,find(Target(:,1) == Response(:,1)));
        T2_R_idx = intersect(repeat_R2_idx,find(Target(:,2) == Response(:,2)));
        T1_M_idx = intersect(mirror_R1_idx,find(Target(:,1) == Response(:,2)));
        T2_M_idx = intersect(mirror_R2_idx,find(Target(:,2) == Response(:,1)));
        rate_by_rank_sepday(iM,1,1,iFile) = length(T1_R_idx)/length(repeat_R1_idx);
        rate_by_rank_sepday(iM,2,1,iFile) = length(T2_R_idx)/length(repeat_R2_idx);
        rate_by_rank_sepday(iM,1,2,iFile) = length(T1_M_idx)/length(mirror_R1_idx);
        rate_by_rank_sepday(iM,2,2,iFile) = length(T2_M_idx)/length(mirror_R2_idx);
        rate_by_rank_sepday(iM,1,3,iFile) = (length(T1_R_idx)+length(T1_M_idx))/(length(repeat_R1_idx) + length(mirror_R1_idx));
        rate_by_rank_sepday(iM,2,3,iFile) = (length(T2_R_idx)+length(T2_M_idx))/(length(repeat_R2_idx) + length(mirror_R2_idx));
        
        % rate by chosen item
        for iT = 1:num_items
            repeat_R1_chosen_vec = zeros(1,num_items);
            repeat_R2_chosen_vec = zeros(1,num_items);
            mirror_R1_chosen_vec = zeros(1,num_items);
            mirror_R2_chosen_vec = zeros(1,num_items);
            
            for jT = 1:num_items
                repeat_R1_chosen_vec(jT) = length(intersect(repeat_R1_idx,find(Target(:,1)==iT & Response(:,1)==jT)));
                repeat_R2_chosen_vec(jT) = length(intersect(repeat_R2_idx,find(Target(:,2)==iT & Response(:,2)==jT)));
                mirror_R1_chosen_vec(jT) = length(intersect(mirror_R1_idx,find(Target(:,1)==iT & Response(:,2)==jT)));
                mirror_R2_chosen_vec(jT) = length(intersect(mirror_R2_idx,find(Target(:,2)==iT & Response(:,1)==jT)));
            end
            temp_repeat = repeat_R1_chosen_vec + repeat_R2_chosen_vec;
            temp_mirror = mirror_R1_chosen_vec + mirror_R2_chosen_vec;
            rate_by_chosen_item_sepday(iM,iT,:,1,iFile) = temp_repeat./sum(temp_repeat);
            rate_by_chosen_item_sepday(iM,iT,:,2,iFile) = temp_mirror./sum(temp_mirror);
            rate_by_chosen_item_sepday(iM,iT,:,3,iFile) = (temp_repeat+temp_mirror)./(sum(temp_repeat) + sum(temp_mirror));
        end
        
        % rate by chosen rank
        repeat_R1_chosen_T1 = length(intersect(repeat_R1_idx,find(Target(:,1) == Response(:,1))));
        repeat_R1_chosen_T2 = length(intersect(repeat_R1_idx,find(Target(:,1) == Response(:,2))));
        repeat_R2_chosen_T1 = length(intersect(repeat_R2_idx,find(Target(:,2) == Response(:,1))));
        repeat_R2_chosen_T2 = length(intersect(repeat_R2_idx,find(Target(:,2) == Response(:,2))));
        
        mirror_R1_chosen_T1 = length(intersect(mirror_R1_idx,find(Target(:,1) == Response(:,2))));
        mirror_R1_chosen_T2 = length(intersect(mirror_R1_idx,find(Target(:,2) == Response(:,2))));
        mirror_R2_chosen_T1 = length(intersect(mirror_R2_idx,find(Target(:,2) == Response(:,1))));
        mirror_R2_chosen_T2 = length(intersect(mirror_R2_idx,find(Target(:,1) == Response(:,1))));
        
        rate_by_chosen_rank_sepday(iM,1,1,1,iFile) = repeat_R1_chosen_T1/(repeat_R1_chosen_T1+repeat_R1_chosen_T2);
        rate_by_chosen_rank_sepday(iM,1,2,1,iFile) = repeat_R1_chosen_T2/(repeat_R1_chosen_T1+repeat_R1_chosen_T2);
        rate_by_chosen_rank_sepday(iM,2,1,1,iFile) = repeat_R2_chosen_T1/(repeat_R2_chosen_T1+repeat_R2_chosen_T2);
        rate_by_chosen_rank_sepday(iM,2,2,1,iFile) = repeat_R2_chosen_T2/(repeat_R2_chosen_T1+repeat_R2_chosen_T2);
        
        rate_by_chosen_rank_sepday(iM,1,1,2,iFile) = mirror_R2_chosen_T1/(mirror_R2_chosen_T1+mirror_R2_chosen_T2);
        rate_by_chosen_rank_sepday(iM,1,2,2,iFile) = mirror_R2_chosen_T2/(mirror_R2_chosen_T1+mirror_R2_chosen_T2);
        rate_by_chosen_rank_sepday(iM,2,1,2,iFile) = mirror_R1_chosen_T2/(mirror_R1_chosen_T1+mirror_R1_chosen_T2);
        rate_by_chosen_rank_sepday(iM,2,2,2,iFile) = mirror_R1_chosen_T1/(mirror_R1_chosen_T1+mirror_R1_chosen_T2);
        
        rate_by_chosen_rank_sepday(iM,1,1,3,iFile) = (repeat_R1_chosen_T1+mirror_R2_chosen_T1)/(repeat_R1_chosen_T1+mirror_R2_chosen_T1+repeat_R1_chosen_T2+mirror_R2_chosen_T2);
        rate_by_chosen_rank_sepday(iM,1,2,3,iFile) = (repeat_R1_chosen_T2+mirror_R2_chosen_T2)/(repeat_R1_chosen_T1+mirror_R2_chosen_T1+repeat_R1_chosen_T2+mirror_R2_chosen_T2);
        rate_by_chosen_rank_sepday(iM,2,1,3,iFile) = (repeat_R2_chosen_T1+mirror_R1_chosen_T2)/(repeat_R2_chosen_T2+mirror_R1_chosen_T1+repeat_R2_chosen_T1+mirror_R1_chosen_T2);
        rate_by_chosen_rank_sepday(iM,2,2,3,iFile) = (repeat_R2_chosen_T2+mirror_R1_chosen_T1)/(repeat_R2_chosen_T2+mirror_R1_chosen_T1+repeat_R2_chosen_T1+mirror_R1_chosen_T2);
    end
    
    % rate by rank
    for i = 1:num_ranks
        for j = 1:3
            rate_by_rank(iM,i,j) = mean(squeeze(rate_by_rank_sepday(iM,i,j,:)));
            rate_by_rank(iM,i,j) = mean(squeeze(rate_by_rank_sepday(iM,i,j,:)));
        end
    end
    
    % rate by chosen item
    for i = 1:num_items
        for j = 1:num_items
            for k = 1:3
                rate_by_chosen_item(iM,i,j,k) = mean(squeeze(rate_by_chosen_item_sepday(iM,i,j,k,:)));
            end
        end
    end
    
    % rate by chosen rank
    for i = 1:num_ranks
        for j = 1:num_ranks
            for k = 1:3
                rate_by_chosen_rank(iM,i,j,k) = mean(squeeze(rate_by_chosen_rank_sepday(iM,i,j,k,:)));
            end
        end
    end
    
end
save(sprintf('%s\\Groot_Behav_2Seq.mat',savedir),'rate_*','dim_labels*')
