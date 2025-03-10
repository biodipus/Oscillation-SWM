%% mapping SPK FW (Ocean), new colormap
addpath('C:\Users\DELL\Documents\scripts\export_fig-master');
load('G:\PaperPrep\SPK_TH_FW_array_dist\new_color_maps.mat');
load('C:\Users\DELL\Documents\FW_data\ocean_data\ocean_pos_labels_0915_2021.mat','RAS_sorted');
RAS_full = [RAS_sorted(1:143,:);[NaN,NaN,NaN];RAS_sorted(144:157,:);[NaN,NaN,NaN];[NaN,NaN,NaN];RAS_sorted(158:end,:)];

loaddir = 'C:\Users\DELL\Documents\FW_data\sharing\Unit_PEV\NewJuly';
spike_dir = 'G:\FW_data\NewJuly\ocean_data\sorting';
idx_dir = 'C:\Users\DELL\Documents\FW_data\sharing\Unit_PEV\NewJuly\spk_ch_info';

load('C:\Users\DELL\Documents\FW_data\ocean_data\preprocessed\new_area_codes_July2021.mat','good_chans','area_code');
F_idx = intersect(setdiff(good_chans,12),find(area_code==2 | area_code==3));
RAS_F = RAS_full(F_idx,:);

file_tags = {'0709','0710','0713','0714','0715','0720','0721','0722','0723','0726','0727','0729','0730'};
num_files = length(file_tags);
num_chans = length(F_idx);

% plot properties
fig_pos = [622,585,369,290];
y_bounds = [-2,28];
outer_size = 38;
inner_size = 38;
panel_str_set = {'SPK Ent','SPK Mem1','TH Ent', 'SPK Mem2'};
color_hold_set = nan(num_chans,4,num_files);
fig_dir = 'G:\PaperPrep\PostNatNeuroRevision\SPK_TH_FW_array_dist\figs_newCM';
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

for iFile = 1:num_files
    curr_tag = file_tags{iFile};
    load(sprintf('%s\\ocean2021-07-%s\\Neuron.mat',spike_dir,curr_tag(end-1:end)),'Spk_channel')
    load(sprintf('%s\\%s_Frontal_idx.mat',idx_dir,curr_tag),'chIn256');
    
    [~, uniqueIdx] = unique(Spk_channel);
    dupeIdx = ismember( Spk_channel, Spk_channel( setdiff( 1:numel(Spk_channel), uniqueIdx ) ) );
    dupes = unique(Spk_channel(dupeIdx)); dupe_locs = find(dupeIdx);

    Spk_chan_trim = setdiff(chIn256,dupes);
    [~,uF_sel] = intersect(F_idx,Spk_chan_trim);

    uF_idx_full = chIn256;
    [~,uF_sel2] = setdiff(uF_idx_full,dupes); % subselect Spk FW
    
    color_hold = nan(num_chans,4); % dim2: entry/mem1/TH-entry/mem2
    
    % Entry FW SPK
    load(sprintf('%s\\tdr_frontal_SPK\\tdr_frontal\\tdr_sens\\%s_weight_enc250.mat',loaddir,curr_tag),'enc_weight') % SPK FW
    SPK_Ent_FW = enc_weight(uF_sel2);
    color_hold(uF_sel,1) = SPK_Ent_FW;
    
    % Mem1 Post-Delay
    load(sprintf('%s\\tdr_frontal_SPK\\tdr_frontal\\tdr_mem_delay\\%s_weight_mem_delay.mat',loaddir,curr_tag),'r1_weight') % SPK FW
    SPK_S1_FW = r1_weight(uF_sel2);
    color_hold(uF_sel,2) = SPK_S1_FW;
    
    % Entry FW TH
    load(sprintf('%s\\tdr_frontal_LFP\\tdr_frontal\\tdr_sens\\%s_weight_enc250_4_8.mat',loaddir,curr_tag),'enc_weight') % TH FW
    TH_FW = enc_weight(uF_sel);
    color_hold(uF_sel,3) = TH_FW;
    
    % Mem2 Post-Delay
    load(sprintf('%s\\tdr_frontal_SPK\\tdr_frontal\\tdr_mem_delay\\%s_weight_mem_delay.mat',loaddir,curr_tag),'r2_weight') % SPK FW
    SPK_S2_FW = r2_weight(uF_sel2);
    color_hold(uF_sel,4) = SPK_S2_FW;
    color_hold_set(:,:,iFile) = color_hold;
end

color_hold_ave = squeeze(nanmean(color_hold_set,3));

for iP = 1:4
    h = figure; set(h,'Position',fig_pos)
        panel_str = panel_str_set{iP};
        scatter(-RAS_F(:,2),RAS_F(:,1),outer_size,[0.5,0.5,0.5]); ylim(y_bounds);
        hold on; scatter(-RAS_F(:,2),RAS_F(:,1),inner_size,color_hold_ave(:,iP),'filled','MarkerEdgeColor',[0.5,0.5,0.5]); 
        colorbar;
        title(panel_str);

    if iP == 1
        colormap(SPK_Ent_map);
    elseif iP == 2
        colormap(SPK_S1_map);
    elseif iP == 3
        colormap(TH_Ent_map);
    elseif iP == 4
        colormap(SPK_S2_map);
    end
    caxis([0.075,0.2]);
    set(gca,'FontSize',10);
    saveas(h,sprintf('%s\\Ocean13_SPK_TH_FW_array_dist_ZC_%s.fig',fig_dir,panel_str));
    export_fig(sprintf('%s\\Ocean13_SPK_TH_FW_array_dist_ZC_%s.pdf',fig_dir,panel_str),h);
    close(h);
end
save('G:\PaperPrep\PostNatNeuroRevision\SPK_TH_FW_array_dist\Ocean13_SPK_TH_FW_array_dist_hold_newCM.mat',...
    'F_idx','color_hold_set','RAS_F','file_tags');

%% mapping SPK FW (Ocean), rand sel peak shuffle + marginal dist
load('G:\PaperPrep\PostNatNeuroRevision\SPK_TH_FW_array_dist\Ocean13_SPK_TH_FW_array_dist_hold_newCM.mat',...
    'F_idx','color_hold_set','RAS_F','file_tags');

num_shuf = 1000;
num_chans = length(F_idx);
scope_size = 40;
peak_size = 4;

pos_hold_SPK_Ent = zeros(num_shuf,2);
pos_hold_SPK_Mem1 = zeros(num_shuf,2);
pos_hold_TH_Ent = zeros(num_shuf,2);
pos_hold_SPK_Mem2 = zeros(num_shuf,2);

% color_hold_ave = squeeze(nanmean(color_hold_set,3));

for iS = 1:num_shuf
    scope_idx = randperm(num_chans,scope_size);
    color_hold_ave = squeeze(nanmean(color_hold_set(scope_idx,:,:),3));
    FW_in_scope = color_hold_ave;
%     FW_in_scope = color_hold_ave(scope_idx,:);
    RAS_in_scope = RAS_F(scope_idx,:);
    
    [~,idx] = sort(FW_in_scope(:,1),'descend','MissingPlacement','last');
    pos_hold_SPK_Ent(iS,:) = mean(RAS_in_scope(idx(1:peak_size),1:2),1);
    
    [~,idx] = sort(FW_in_scope(:,2),'descend','MissingPlacement','last');
    pos_hold_SPK_Mem1(iS,:) = mean(RAS_in_scope(idx(1:peak_size),1:2),1);
    
    [~,idx] = sort(FW_in_scope(:,3),'descend','MissingPlacement','last');
    pos_hold_TH_Ent(iS,:) = mean(RAS_in_scope(idx(1:peak_size),1:2),1);
    
    [~,idx] = sort(FW_in_scope(:,4),'descend','MissingPlacement','last');
    pos_hold_SPK_Mem2(iS,:) = mean(RAS_in_scope(idx(1:peak_size),1:2),1);
end

SPK_Ent_pos = [-pos_hold_SPK_Ent(:,2),pos_hold_SPK_Ent(:,1)];
SPK_S1_pos = [-pos_hold_SPK_Mem1(:,2),pos_hold_SPK_Mem1(:,1)];
TH_Ent_pos = [-pos_hold_TH_Ent(:,2),pos_hold_TH_Ent(:,1)];
SPK_S2_pos = [-pos_hold_SPK_Mem2(:,2),pos_hold_SPK_Mem2(:,1)];

save('G:\PaperPrep\PostNatNeuroRevision\SPK_TH_FW_array_dist\Ocean13_scope40_peak4_Ent_Mem_pos.mat','*_pos');

%% mapping SPK FW (Groot), new colormap
addpath('C:\Users\DELL\Documents\scripts\export_fig-master');
load('G:\PaperPrep\SPK_TH_FW_array_dist\new_color_maps.mat');
load('E:\ZC_data\groot_data\groot_pos_labels.mat','RAS_sorted');

RAS_full = [RAS_sorted(1:143,:);[NaN,NaN,NaN];RAS_sorted(144:157,:);[NaN,NaN,NaN];[NaN,NaN,NaN];RAS_sorted(158:end,:)];

loaddir = 'E:\ZC_data\groot_data_switch\sharing';
spike_dir = 'E:\ZC_data\groot_data_switch\Groot_switch_sorting\sorting';
idx_dir = 'E:\ZC_data\groot_data_switch\sharing\groot_switch_ch_idx';

load('C:\Users\DELL\Documents\ZC_data\groot_data\new_area_codes_Mar2024.mat','good_chans','area_code');
F_idx = intersect(good_chans,find(area_code==2 | area_code==3));
RAS_F = RAS_full(F_idx,:);

file_tags = {'0415','0507','0509','0513','0514','0515','0516','0520','0521','0522'};
num_files = length(file_tags);
num_chans = length(F_idx);

% plot properties
fig_pos = [622,585,369,290];
y_bounds = [-22,-2];
outer_size = 38;
inner_size = 38;
panel_str_set = {'SPK Ent','SPK Mem1','TH Ent', 'SPK Mem2'};
color_hold_set = nan(num_chans,4,num_files);
fig_dir = 'G:\PaperPrep\PostNatNeuroRevision\SPK_TH_FW_array_dist\figs_newCM';
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

for iFile = 1:num_files
    curr_tag = file_tags{iFile};
    load(sprintf('%s\\groot2024-%s-%s\\Neuron.mat',spike_dir,curr_tag(end-3:end-2),curr_tag(end-1:end)),'Spk_channel')
    Spk_channel = Spk_channel(:,1);
    load(sprintf('%s\\%s_ch_idx.mat',idx_dir,curr_tag),'chIn256');
    
    [~, uniqueIdx] = unique(Spk_channel);
    dupeIdx = ismember( Spk_channel, Spk_channel( setdiff( 1:numel(Spk_channel), uniqueIdx ) ) );
    dupes = unique(Spk_channel(dupeIdx)); dupe_locs = find(dupeIdx);

    Spk_chan_trim = setdiff(chIn256,dupes);
    [~,uF_sel] = intersect(F_idx,Spk_chan_trim);

    uF_idx_full = chIn256;
    [~,uF_sel2] = setdiff(uF_idx_full,dupes); % subselect Spk FW
    
    color_hold = nan(num_chans,4); % dim2: entry/mem1/TH-entry/mem2
    
    % Entry FW SPK
    load(sprintf('%s\\groot_sw_weight_enc250\\%s-%s_weight_enc250.mat',loaddir,curr_tag(end-3:end-2),curr_tag(end-1:end)),'enc_weight');
    SPK_Ent_FW = enc_weight(uF_sel2);
    color_hold(uF_sel,1) = SPK_Ent_FW;
    
    % Mem1 Post-Delay
    load(sprintf('%s\\groot_swith_11days_weight_mem_delay(1)\\%s-%s_weight_mem_delay.mat',loaddir,curr_tag(end-3:end-2),curr_tag(end-1:end)),'r1_weight') % SPK FW
    SPK_S1_FW = r1_weight(uF_sel2);
    color_hold(uF_sel,2) = SPK_S1_FW;
    
    % Entry FW TH
    load(sprintf('%s\\groot_sw_weight_enc250_4_8\\%s-%s_weight_enc250_4_8.mat',loaddir,curr_tag(end-3:end-2),curr_tag(end-1:end)),'enc_weight');
    TH_FW = enc_weight(uF_sel);
    color_hold(uF_sel,3) = TH_FW;
    
    % Mem2 Post-Delay
    load(sprintf('%s\\groot_swith_11days_weight_mem_delay(1)\\%s-%s_weight_mem_delay.mat',loaddir,curr_tag(end-3:end-2),curr_tag(end-1:end)),'r2_weight') % SPK FW
    SPK_S2_FW = r2_weight(uF_sel2);
    color_hold(uF_sel,4) = SPK_S2_FW;
    color_hold_set(:,:,iFile) = color_hold;
end

color_hold_ave = squeeze(nanmean(color_hold_set,3));

for iP = 1:4
    h = figure; set(h,'Position',fig_pos)
        panel_str = panel_str_set{iP};
        scatter(-RAS_F(:,2),RAS_F(:,1),outer_size,[0.5,0.5,0.5]); ylim(y_bounds);
        hold on; scatter(-RAS_F(:,2),RAS_F(:,1),inner_size,color_hold_ave(:,iP),'filled','MarkerEdgeColor',[0.5,0.5,0.5]); 
        colorbar;
        title(panel_str);

    if iP == 1
        colormap(SPK_Ent_map);
    elseif iP == 2
        colormap(SPK_S1_map);
    elseif iP == 3
        colormap(TH_Ent_map);
    elseif iP == 4
        colormap(SPK_S2_map);
    end
    caxis([0.075,0.2]);
    set(gca,'FontSize',10);
    saveas(h,sprintf('%s\\Groot10_SPK_TH_FW_array_dist_ZC_%s.fig',fig_dir,panel_str));
    export_fig(sprintf('%s\\Groot10_SPK_TH_FW_array_dist_ZC_%s.pdf',fig_dir,panel_str),h);
    close(h);
end
save('G:\PaperPrep\PostNatNeuroRevision\SPK_TH_FW_array_dist\Groot10_SPK_TH_FW_array_dist_hold_newCM.mat',...
    'F_idx','color_hold_set','RAS_F','file_tags');

%% mapping SPK FW (Groot), rand sel peak shuffle + marginal dist
load('G:\PaperPrep\PostNatNeuroRevision\SPK_TH_FW_array_dist\Groot10_SPK_TH_FW_array_dist_hold_newCM.mat',...
    'F_idx','color_hold_set','RAS_F','file_tags');
num_shuf = 1000;
num_chans = length(F_idx);
scope_size = 40;
peak_size = 4;

pos_hold_SPK_Ent = zeros(num_shuf,2);
pos_hold_SPK_Mem1 = zeros(num_shuf,2);
pos_hold_TH_Ent = zeros(num_shuf,2);
pos_hold_SPK_Mem2 = zeros(num_shuf,2);

% color_hold_ave = squeeze(nanmean(color_hold_set,3));

for iS = 1:num_shuf
    scope_idx = randperm(num_chans,scope_size);
    color_hold_ave = squeeze(nanmean(color_hold_set(scope_idx,:,:),3));
    FW_in_scope = color_hold_ave;
%     FW_in_scope = color_hold_ave(scope_idx,:);
    RAS_in_scope = RAS_F(scope_idx,:);
    
    [~,idx] = sort(FW_in_scope(:,1),'descend','MissingPlacement','last');
    pos_hold_SPK_Ent(iS,:) = mean(RAS_in_scope(idx(1:peak_size),1:2),1);
    
    [~,idx] = sort(FW_in_scope(:,2),'descend','MissingPlacement','last');
    pos_hold_SPK_Mem1(iS,:) = mean(RAS_in_scope(idx(1:peak_size),1:2),1);
    
    [~,idx] = sort(FW_in_scope(:,3),'descend','MissingPlacement','last');
    pos_hold_TH_Ent(iS,:) = mean(RAS_in_scope(idx(1:peak_size),1:2),1);
    
    [~,idx] = sort(FW_in_scope(:,4),'descend','MissingPlacement','last');
    pos_hold_SPK_Mem2(iS,:) = mean(RAS_in_scope(idx(1:peak_size),1:2),1);
end

SPK_Ent_pos = [-pos_hold_SPK_Ent(:,2),pos_hold_SPK_Ent(:,1)];
SPK_S1_pos = [-pos_hold_SPK_Mem1(:,2),pos_hold_SPK_Mem1(:,1)];
TH_Ent_pos = [-pos_hold_TH_Ent(:,2),pos_hold_TH_Ent(:,1)];
SPK_S2_pos = [-pos_hold_SPK_Mem2(:,2),pos_hold_SPK_Mem2(:,1)];

save('G:\PaperPrep\PostNatNeuroRevision\SPK_TH_FW_array_dist\Groot10_scope40_peak4_Ent_Mem_pos.mat','*_pos');
