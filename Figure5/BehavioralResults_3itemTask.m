%% behavior tally: Ocean 3Seq, correct rate by chosen rank
savedir = 'G:\PaperPrep';
file_tags = {'0827','0830','0831','0901','0902','0906','0907','0908','0909','0910'};
behav_dir = 'G:\FW_data\NewJuly\ocean_data\sorting';

num_ranks = 3;
num_days = 10;
num_files = length(file_tags);

rate_by_chosen_rank_sepday = zeros(num_ranks,num_ranks,num_days); % dim 3: repeat/mirror/all
dim_labels_chosen_sepday = {'rank_actual','rank_chosen','days'};

for iFile = 1:num_files
    curr_tag = file_tags{iFile};
    
    load(sprintf('%s\\ocean2021-%s-%s\\TrialInfo.mat',behav_dir,curr_tag(1:2),curr_tag(3:4)),'TrialInfo')
    
    Rule = TrialInfo.Rule;
    CueFrameOnTime = TrialInfo.CueFrameOnTime;
    Target = TrialInfo.Target;
    Response = TrialInfo.Response;
    
    usable_idx = find(~isnan(Rule) & ~isnan(CueFrameOnTime) & ~isnan(Target(:,3)));
%     repeat_R1_idx = intersect(usable_idx,find(~isnan(Target(:,3)) & ~isnan(Response(:,1))));
%     repeat_R2_idx = intersect(usable_idx,find(~isnan(Target(:,3)) & ~isnan(Response(:,2))));
%     repeat_R3_idx = intersect(usable_idx,find(~isnan(Target(:,3)) & ~isnan(Response(:,3))));
    
    repeat_R1_chosen_T1 = length(intersect(usable_idx,find(Target(:,1) == Response(:,1))));
    repeat_R1_chosen_T2 = length(intersect(usable_idx,find(Target(:,1) == Response(:,2))));
    repeat_R1_chosen_T3 = length(intersect(usable_idx,find(Target(:,1) == Response(:,3))));
    
    repeat_R2_chosen_T1 = length(intersect(usable_idx,find(Target(:,2) == Response(:,1))));
    repeat_R2_chosen_T2 = length(intersect(usable_idx,find(Target(:,2) == Response(:,2))));
    repeat_R2_chosen_T3 = length(intersect(usable_idx,find(Target(:,2) == Response(:,3))));
    
    repeat_R3_chosen_T1 = length(intersect(usable_idx,find(Target(:,3) == Response(:,1))));
    repeat_R3_chosen_T2 = length(intersect(usable_idx,find(Target(:,3) == Response(:,2))));
    repeat_R3_chosen_T3 = length(intersect(usable_idx,find(Target(:,3) == Response(:,3))));
        
    rate_by_chosen_rank_sepday(1,1,iFile) = repeat_R1_chosen_T1/(repeat_R1_chosen_T1+repeat_R1_chosen_T2+repeat_R1_chosen_T3);
    rate_by_chosen_rank_sepday(1,2,iFile) = repeat_R1_chosen_T2/(repeat_R1_chosen_T1+repeat_R1_chosen_T2+repeat_R1_chosen_T3);
    rate_by_chosen_rank_sepday(1,3,iFile) = repeat_R1_chosen_T3/(repeat_R1_chosen_T1+repeat_R1_chosen_T2+repeat_R1_chosen_T3);
    
    rate_by_chosen_rank_sepday(2,1,iFile) = repeat_R2_chosen_T1/(repeat_R2_chosen_T1+repeat_R2_chosen_T2+repeat_R2_chosen_T3);
    rate_by_chosen_rank_sepday(2,2,iFile) = repeat_R2_chosen_T2/(repeat_R2_chosen_T1+repeat_R2_chosen_T2+repeat_R2_chosen_T3);
    rate_by_chosen_rank_sepday(2,3,iFile) = repeat_R2_chosen_T3/(repeat_R2_chosen_T1+repeat_R2_chosen_T2+repeat_R2_chosen_T3);
    
    rate_by_chosen_rank_sepday(3,1,iFile) = repeat_R3_chosen_T1/(repeat_R3_chosen_T1+repeat_R3_chosen_T2+repeat_R3_chosen_T3);
    rate_by_chosen_rank_sepday(3,2,iFile) = repeat_R3_chosen_T2/(repeat_R3_chosen_T1+repeat_R3_chosen_T2+repeat_R3_chosen_T3);
    rate_by_chosen_rank_sepday(3,3,iFile) = repeat_R3_chosen_T3/(repeat_R3_chosen_T1+repeat_R3_chosen_T2+repeat_R3_chosen_T3);
end
save(sprintf('%s\\Ocean_Behav_3Seq.mat',savedir),'rate_*','dim_labels*')

%% Groot behavior tally: Ocean 3Seq, correct rate by chosen rank
savedir = 'G:\PaperPrep\PostNatNeuroRevision';
file_tags = {'0311','0314','0315','0316','0319','0321','0322','0323','0325','0326','0327'};
behav_dir = 'C:\Users\DELL\Documents\ZC_data\groot_data\sorting';

num_ranks = 3;
num_days = 11;
num_files = length(file_tags);

rate_by_chosen_rank_sepday = zeros(num_ranks,num_ranks,num_days); % dim 3: repeat/mirror/all
dim_labels_chosen_sepday = {'rank_actual','rank_chosen','days'};

for iFile = 1:num_files
    curr_tag = file_tags{iFile};
    
    load(sprintf('%s\\groot2024-%s-%s\\TrialInfo.mat',behav_dir,curr_tag(1:2),curr_tag(3:4)),'TrialInfo')
    
    Rule = TrialInfo.Rule;
    CueFrameOnTime = TrialInfo.CueFrameOnTime;
    Target = TrialInfo.Target;
    Response = TrialInfo.Response;
    
    usable_idx = find(~isnan(Rule) & ~isnan(CueFrameOnTime) & ~isnan(Target(:,3)));
%     repeat_R1_idx = intersect(usable_idx,find(~isnan(Target(:,3)) & ~isnan(Response(:,1))));
%     repeat_R2_idx = intersect(usable_idx,find(~isnan(Target(:,3)) & ~isnan(Response(:,2))));
%     repeat_R3_idx = intersect(usable_idx,find(~isnan(Target(:,3)) & ~isnan(Response(:,3))));
    
    repeat_R1_chosen_T1 = length(intersect(usable_idx,find(Target(:,1) == Response(:,1))));
    repeat_R1_chosen_T2 = length(intersect(usable_idx,find(Target(:,1) == Response(:,2))));
    repeat_R1_chosen_T3 = length(intersect(usable_idx,find(Target(:,1) == Response(:,3))));
    
    repeat_R2_chosen_T1 = length(intersect(usable_idx,find(Target(:,2) == Response(:,1))));
    repeat_R2_chosen_T2 = length(intersect(usable_idx,find(Target(:,2) == Response(:,2))));
    repeat_R2_chosen_T3 = length(intersect(usable_idx,find(Target(:,2) == Response(:,3))));
    
    repeat_R3_chosen_T1 = length(intersect(usable_idx,find(Target(:,3) == Response(:,1))));
    repeat_R3_chosen_T2 = length(intersect(usable_idx,find(Target(:,3) == Response(:,2))));
    repeat_R3_chosen_T3 = length(intersect(usable_idx,find(Target(:,3) == Response(:,3))));
        
    rate_by_chosen_rank_sepday(1,1,iFile) = repeat_R1_chosen_T1/(repeat_R1_chosen_T1+repeat_R1_chosen_T2+repeat_R1_chosen_T3);
    rate_by_chosen_rank_sepday(1,2,iFile) = repeat_R1_chosen_T2/(repeat_R1_chosen_T1+repeat_R1_chosen_T2+repeat_R1_chosen_T3);
    rate_by_chosen_rank_sepday(1,3,iFile) = repeat_R1_chosen_T3/(repeat_R1_chosen_T1+repeat_R1_chosen_T2+repeat_R1_chosen_T3);
    
    rate_by_chosen_rank_sepday(2,1,iFile) = repeat_R2_chosen_T1/(repeat_R2_chosen_T1+repeat_R2_chosen_T2+repeat_R2_chosen_T3);
    rate_by_chosen_rank_sepday(2,2,iFile) = repeat_R2_chosen_T2/(repeat_R2_chosen_T1+repeat_R2_chosen_T2+repeat_R2_chosen_T3);
    rate_by_chosen_rank_sepday(2,3,iFile) = repeat_R2_chosen_T3/(repeat_R2_chosen_T1+repeat_R2_chosen_T2+repeat_R2_chosen_T3);
    
    rate_by_chosen_rank_sepday(3,1,iFile) = repeat_R3_chosen_T1/(repeat_R3_chosen_T1+repeat_R3_chosen_T2+repeat_R3_chosen_T3);
    rate_by_chosen_rank_sepday(3,2,iFile) = repeat_R3_chosen_T2/(repeat_R3_chosen_T1+repeat_R3_chosen_T2+repeat_R3_chosen_T3);
    rate_by_chosen_rank_sepday(3,3,iFile) = repeat_R3_chosen_T3/(repeat_R3_chosen_T1+repeat_R3_chosen_T2+repeat_R3_chosen_T3);
end
save(sprintf('%s\\Groot_Behav_3Seq.mat',savedir),'rate_*','dim_labels*')