%% mapping SPK FW (Ocean), Seq3
load('C:\Users\DELL\Documents\FW_data\ocean_data\ocean_pos_labels_0915_2021.mat','RAS_sorted');
RAS_full = [RAS_sorted(1:143,:);[NaN,NaN,NaN];RAS_sorted(144:157,:);[NaN,NaN,NaN];[NaN,NaN,NaN];RAS_sorted(158:end,:)];
load('G:\PaperPrep\SPK_TH_FW_array_dist\new_color_maps.mat');

loaddir = 'C:\Users\DELL\Documents\FW_data\sharing\Unit_PEV\PostSept';
spike_dir = 'G:\FW_data\NewJuly\ocean_data\sorting';
idx_dir = 'C:\Users\DELL\Documents\FW_data\sharing\Unit_PEV\PostSept\tdr_frontal_spk_L3\tdr_sens';

load('C:\Users\DELL\Documents\FW_data\ocean_data\preprocessed\new_area_codes_July2021.mat','good_chans','area_code');
F_idx = intersect(setdiff(good_chans,12),find(area_code==2 | area_code==3));
RAS_F = RAS_full(F_idx,:);

file_tags = {'0827','0830','0831','0901','0902','0906','0907','0908','0909','0910'};
num_files = length(file_tags);
num_chans = length(F_idx);

% plot properties
fig_pos = [622,585,369,290];
y_bounds = [-2,28];
outer_size = 38;
inner_size = 38;
panel_str_set = {'SPK Mem1','SPK Ent','SPK Mem2','TH Ent','SPK Mem3'};
color_hold_set = nan(num_chans,5,num_files);
fig_dir = 'G:\PaperPrep\PostNatNeuroRevision\SPK_TH_FW_array_dist\figs_3Seq';
png_dir = 'G:\PaperPrep\PostNatNeuroRevision\SPK_TH_FW_array_dist\pngs_3Seq';
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end
if ~exist(png_dir,'dir')
    mkdir(png_dir);
end

for iFile = 1:num_files
    curr_tag = file_tags{iFile};
    load(sprintf('%s\\ocean2021-%s-%s\\Neuron.mat',spike_dir,curr_tag(1:2),curr_tag(end-1:end)),'Spk_channel')
    load(sprintf('%s\\%s-%s_weight_enc250.mat',idx_dir,curr_tag(end-3:end-2),curr_tag(end-1:end)),'ch_info_all');
    chIn256 = double(cell2mat(ch_info_all(:,2)));
    
    [~, uniqueIdx] = unique(Spk_channel);
    dupeIdx = ismember( Spk_channel, Spk_channel( setdiff( 1:numel(Spk_channel), uniqueIdx ) ) );
    dupes = unique(Spk_channel(dupeIdx)); dupe_locs = find(dupeIdx);

    Spk_chan_trim = setdiff(chIn256,dupes);
    [~,uF_sel] = intersect(F_idx,Spk_chan_trim);

    uF_idx_full = chIn256;
    [~,uF_sel2] = setdiff(uF_idx_full,dupes); % subselect Spk FW
    
    color_hold = nan(num_chans,5); % dim2: entry/mem1/TH-entry/mem2
    
    % Mem1 Post-Delay
    load(sprintf('%s\\tdr_frontal_spk_L3\\tdr_mem_delay\\%s-%s_weight_mem_delay.mat',loaddir,curr_tag(end-3:end-2),curr_tag(end-1:end)),'r1_weight') % SPK FW
    SPK_S1_FW = r1_weight(uF_sel2);
    color_hold(uF_sel,1) = SPK_S1_FW;
    
    % Entry FW SPK
    load(sprintf('%s\\tdr_frontal_spk_L3\\tdr_sens\\%s-%s_weight_enc250.mat',loaddir,curr_tag(end-3:end-2),curr_tag(end-1:end)),'enc_weight') % SPK FW
    SPK_Ent_FW = enc_weight(uF_sel2);
    color_hold(uF_sel,2) = SPK_Ent_FW;
    
    % Mem2 Post-Delay
    load(sprintf('%s\\tdr_frontal_spk_L3\\tdr_mem_delay\\%s-%s_weight_mem_delay.mat',loaddir,curr_tag(end-3:end-2),curr_tag(end-1:end)),'r2_weight') % SPK FW
    SPK_S2_FW = r2_weight(uF_sel2);
    color_hold(uF_sel,3) = SPK_S2_FW;
    
    % Entry FW TH
    load(sprintf('%s\\tdr_frontal_LFP_L3\\tdr_sens\\%s-%s_weight_enc250_4_8.mat',loaddir,curr_tag(end-3:end-2),curr_tag(end-1:end)),'enc_weight') % TH FW
    TH_FW = enc_weight(uF_sel);
    color_hold(uF_sel,4) = TH_FW;
    
    % Mem3 Post-Delay
    load(sprintf('%s\\tdr_frontal_spk_L3\\tdr_mem_delay\\%s-%s_weight_mem_delay.mat',loaddir,curr_tag(end-3:end-2),curr_tag(end-1:end)),'r3_weight') % SPk FW
    SPK_S3_FW = r3_weight(uF_sel2);
    color_hold(uF_sel,5) = SPK_S3_FW;
    color_hold_set(:,:,iFile) = color_hold;
    
end

color_hold_ave = squeeze(nanmean(color_hold_set,3));
for iP = 1:5
    h = figure; set(h,'Position',fig_pos)
        panel_str = panel_str_set{iP};
        scatter(-RAS_F(:,2),RAS_F(:,1),outer_size,[0.5,0.5,0.5]); ylim(y_bounds);
        hold on; scatter(-RAS_F(:,2),RAS_F(:,1),inner_size,color_hold_ave(:,iP),'filled','MarkerEdgeColor',[0.5,0.5,0.5]); 
        colorbar;
        title(panel_str);

    if iP == 1
        colormap(SPK_S1_map);
    elseif iP == 2
        colormap(SPK_Ent_map);
    elseif iP == 3
        colormap(SPK_S2_map);
    elseif iP == 4
        colormap(TH_Ent_map);
    elseif iP == 5
        colormap(SPK_S3_map);
    end
    caxis([0.075,0.2]);
    set(gca,'FontSize',10);
    saveas(h,sprintf('%s\\Ocean13_3Seq_SPK_TH_FW_array_dist_ZC_%s.fig',fig_dir,panel_str));
    export_fig(sprintf('%s\\Ocean13_3Seq_SPK_TH_FW_array_dist_ZC_%s.pdf',fig_dir,panel_str),h);
    close(h);
end
    
save('G:\PaperPrep\PostNatNeuroRevision\SPK_TH_FW_array_dist\Ocean13_SPK_TH_FW_array_dist_hold_3Seq.mat',...
    'F_idx','color_hold_set','RAS_F','file_tags');
    
%% mapping SPK FW (Ocean), Seq3, rand sel peak shuffle + marginal dist
load('G:\PaperPrep\PostNatNeuroRevision\SPK_TH_FW_array_dist\Ocean13_SPK_TH_FW_array_dist_hold_3Seq.mat',...
    'F_idx','color_hold_set','RAS_F','file_tags');

num_shuf = 1000;
num_chans = length(F_idx);
scope_size = 40;
peak_size = 4;

pos_hold_SPK_Mem1 = zeros(num_shuf,2);
pos_hold_SPK_Ent = zeros(num_shuf,2);
pos_hold_SPK_Mem2 = zeros(num_shuf,2);
pos_hold_TH_Ent = zeros(num_shuf,2);
pos_hold_SPK_Mem3 = zeros(num_shuf,2);

% color_hold_ave = squeeze(nanmean(color_hold_set,3));

for iS = 1:num_shuf
    scope_idx = randperm(num_chans,scope_size);
    color_hold_ave = squeeze(nanmean(color_hold_set(scope_idx,:,:),3));
    FW_in_scope = color_hold_ave;
%     FW_in_scope = color_hold_ave(scope_idx,:);
    RAS_in_scope = RAS_F(scope_idx,:);
    
    [~,idx] = sort(FW_in_scope(:,1),'descend','MissingPlacement','last');
    pos_hold_SPK_Mem1(iS,:) = mean(RAS_in_scope(idx(1:peak_size),1:2),1);
    
    [~,idx] = sort(FW_in_scope(:,2),'descend','MissingPlacement','last');
    pos_hold_SPK_Ent(iS,:) = mean(RAS_in_scope(idx(1:peak_size),1:2),1);
    
    [~,idx] = sort(FW_in_scope(:,3),'descend','MissingPlacement','last');
    pos_hold_SPK_Mem2(iS,:) = mean(RAS_in_scope(idx(1:peak_size),1:2),1);
    
    [~,idx] = sort(FW_in_scope(:,4),'descend','MissingPlacement','last');
    pos_hold_TH_Ent(iS,:) = mean(RAS_in_scope(idx(1:peak_size),1:2),1);
    
    [~,idx] = sort(FW_in_scope(:,5),'descend','MissingPlacement','last');
    pos_hold_SPK_Mem3(iS,:) = mean(RAS_in_scope(idx(1:peak_size),1:2),1);
end

SPK_Ent_pos = [-pos_hold_SPK_Ent(:,2),pos_hold_SPK_Ent(:,1)];
SPK_S1_pos = [-pos_hold_SPK_Mem1(:,2),pos_hold_SPK_Mem1(:,1)];
TH_Ent_pos = [-pos_hold_TH_Ent(:,2),pos_hold_TH_Ent(:,1)];
SPK_S2_pos = [-pos_hold_SPK_Mem2(:,2),pos_hold_SPK_Mem2(:,1)];
SPK_S3_pos = [-pos_hold_SPK_Mem3(:,2),pos_hold_SPK_Mem3(:,1)];
save('G:\PaperPrep\PostNatNeuroRevision\SPK_TH_FW_array_dist\Ocean13_scope40_peak4_3Seq_Ent_Mem_pos.mat','*_pos');

%% mapping SPK FW (Groot), Seq3
load('E:\ZC_data\groot_data\groot_pos_labels.mat','RAS_sorted');
RAS_full = [RAS_sorted(1:143,:);[NaN,NaN,NaN];RAS_sorted(144:157,:);[NaN,NaN,NaN];[NaN,NaN,NaN];RAS_sorted(158:end,:)];
load('G:\PaperPrep\SPK_TH_FW_array_dist\new_color_maps.mat');

loaddir = 'E:\ZC_data\groot_data\sharing';
spike_dir = 'C:\Users\DELL\Documents\ZC_data\groot_data\sorting';
idx_dir = 'E:\ZC_data\groot_data\sharing\13days_fr_L3_win100step50_frontal';

load('C:\Users\DELL\Documents\ZC_data\groot_data\new_area_codes_Mar2024.mat','good_chans','area_code');
F_idx = intersect(good_chans,find(area_code==2 | area_code==3));
RAS_F = RAS_full(F_idx,:);

file_tags = {'0311','0314','0315','0316','0319','0321','0322','0323','0325','0326','0327'};
num_files = length(file_tags);
num_chans = length(F_idx);

% plot properties
fig_pos = [622,585,369,290];
y_bounds = [-22,-2];
outer_size = 38;
inner_size = 38;
panel_str_set = {'SPK Mem1','SPK Ent','SPK Mem2','TH Ent','SPK Mem3'};
color_hold_set = nan(num_chans,5,num_files);
fig_dir = 'G:\PaperPrep\PostNatNeuroRevision\SPK_TH_FW_array_dist\figs_3Seq';
png_dir = 'G:\PaperPrep\PostNatNeuroRevision\SPK_TH_FW_array_dist\pngs_3Seq';
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end
if ~exist(png_dir,'dir')
    mkdir(png_dir);
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
    
    color_hold = nan(num_chans,5); % dim2: mem1/mem2/mem3
    
    % Mem1 Post-Delay
    load(sprintf('%s\\tdr_mem_delay1_LongD\\%s-%s_weight_mem_delay.mat',loaddir,curr_tag(end-3:end-2),curr_tag(end-1:end)),'r1_weight');
    SPK_S1_FW = r1_weight(uF_sel2);
    color_hold(uF_sel,1) = SPK_S1_FW;
    
    % Entry FW SPK
    load(sprintf('%s\\grootL3_weight_enc250\\%s-%s_weight_enc250.mat',loaddir,curr_tag(end-3:end-2),curr_tag(end-1:end)),'enc_weight') % SPK FW
    SPK_Ent_FW = enc_weight(uF_sel2);
    color_hold(uF_sel,2) = SPK_Ent_FW;
    
    % Mem2 Post-Delay
    load(sprintf('%s\\tdr_mem_delay1_LongD\\%s-%s_weight_mem_delay.mat',loaddir,curr_tag(end-3:end-2),curr_tag(end-1:end)),'r2_weight') % SPK FW
    SPK_S2_FW = r2_weight(uF_sel2);
    color_hold(uF_sel,3) = SPK_S2_FW;
    
    % Entry FW TH
    load(sprintf('%s\\grootL3_weight_enc250_4_8\\%s-%s_weight_enc250_4_8.mat',loaddir,curr_tag(end-3:end-2),curr_tag(end-1:end)),'enc_weight') % TH FW
    TH_FW = enc_weight(uF_sel);
    color_hold(uF_sel,4) = TH_FW;
    
    % Mem3 Post-Delay
    load(sprintf('%s\\tdr_mem_delay1_LongD\\%s-%s_weight_mem_delay.mat',loaddir,curr_tag(end-3:end-2),curr_tag(end-1:end)),'r3_weight') % SPk FW
    SPK_S3_FW = r3_weight(uF_sel2);
    color_hold(uF_sel,5) = SPK_S3_FW;
    color_hold_set(:,:,iFile) = color_hold;
    
end

color_hold_ave = squeeze(nanmean(color_hold_set,3));
for iP = 1:5
    h = figure; set(h,'Position',fig_pos)
        panel_str = panel_str_set{iP};
        scatter(-RAS_F(:,2),RAS_F(:,1),outer_size,[0.5,0.5,0.5]); ylim(y_bounds);
        hold on; scatter(-RAS_F(:,2),RAS_F(:,1),inner_size,color_hold_ave(:,iP),'filled','MarkerEdgeColor',[0.5,0.5,0.5]); 
        colorbar;
        title(panel_str);

    if iP == 1
        colormap(SPK_S1_map);
    elseif iP == 2
        colormap(SPK_Ent_map);
    elseif iP == 3
        colormap(SPK_S2_map);
    elseif iP == 4
        colormap(TH_Ent_map);
    elseif iP == 5
        colormap(SPK_S3_map);
    end
    caxis([0.075,0.2]);
    set(gca,'FontSize',10);
    saveas(h,sprintf('%s\\Groot11_3Seq_SPK_TH_FW_array_dist_ZC_%s.fig',fig_dir,panel_str));
    export_fig(sprintf('%s\\Groot11_3Seq_SPK_TH_FW_array_dist_ZC_%s.pdf',fig_dir,panel_str),h);
    close(h);
end
    
save('G:\PaperPrep\PostNatNeuroRevision\SPK_TH_FW_array_dist\Groot11_SPK_TH_FW_array_dist_hold_3Seq.mat',...
    'F_idx','color_hold_set','RAS_F','file_tags');
    
%% mapping SPK FW (Groot), Seq3, rand sel peak shuffle + marginal dist
load('G:\PaperPrep\PostNatNeuroRevision\SPK_TH_FW_array_dist\Groot11_SPK_TH_FW_array_dist_hold_3Seq.mat',...
    'F_idx','color_hold_set','RAS_F','file_tags');
num_shuf = 1000;
num_chans = length(F_idx);
scope_size = 40;
peak_size = 4;

pos_hold_SPK_Mem1 = zeros(num_shuf,2);
pos_hold_SPK_Ent = zeros(num_shuf,2);
pos_hold_SPK_Mem2 = zeros(num_shuf,2);
pos_hold_TH_Ent = zeros(num_shuf,2);
pos_hold_SPK_Mem3 = zeros(num_shuf,2);

% color_hold_ave = squeeze(nanmean(color_hold_set,3));

for iS = 1:num_shuf
    scope_idx = randperm(num_chans,scope_size);
    color_hold_ave = squeeze(nanmean(color_hold_set(scope_idx,:,:),3));
    FW_in_scope = color_hold_ave;
%     FW_in_scope = color_hold_ave(scope_idx,:);
    RAS_in_scope = RAS_F(scope_idx,:);
    
    [~,idx] = sort(FW_in_scope(:,1),'descend','MissingPlacement','last');
    pos_hold_SPK_Mem1(iS,:) = mean(RAS_in_scope(idx(1:peak_size),1:2),1);
    
    [~,idx] = sort(FW_in_scope(:,2),'descend','MissingPlacement','last');
    pos_hold_SPK_Ent(iS,:) = mean(RAS_in_scope(idx(1:peak_size),1:2),1);
    
    [~,idx] = sort(FW_in_scope(:,3),'descend','MissingPlacement','last');
    pos_hold_SPK_Mem2(iS,:) = mean(RAS_in_scope(idx(1:peak_size),1:2),1);
    
    [~,idx] = sort(FW_in_scope(:,4),'descend','MissingPlacement','last');
    pos_hold_TH_Ent(iS,:) = mean(RAS_in_scope(idx(1:peak_size),1:2),1);
    
    [~,idx] = sort(FW_in_scope(:,5),'descend','MissingPlacement','last');
    pos_hold_SPK_Mem3(iS,:) = mean(RAS_in_scope(idx(1:peak_size),1:2),1);
end

SPK_S1_pos = [-pos_hold_SPK_Mem1(:,2),pos_hold_SPK_Mem1(:,1)];
SPK_Ent_pos = [-pos_hold_SPK_Ent(:,2),pos_hold_SPK_Ent(:,1)];
SPK_S2_pos = [-pos_hold_SPK_Mem2(:,2),pos_hold_SPK_Mem2(:,1)];
TH_Ent_pos = [-pos_hold_TH_Ent(:,2),pos_hold_TH_Ent(:,1)];
SPK_S3_pos = [-pos_hold_SPK_Mem3(:,2),pos_hold_SPK_Mem3(:,1)];
save('G:\PaperPrep\PostNatNeuroRevision\SPK_TH_FW_array_dist\Groot11_scope40_peak4_3Seq_Ent_Mem_pos.mat','*_pos');
