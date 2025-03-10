%% individual + concat SFC slide maps, sort by pT1/pT2
savedir = 'G:\PaperPrep\CurrBiolRevision';
vline_pos = 16.5;
x_bounds = [6,32];
fig_pos = [548,50,824,946];

Ocean_dir = 'G:\FW_data\NewJuly\ocean_data\preprocessed\prelim_analyses\SFC_CrossChan_Slide';
Ocean_tags = {'0709','0710','0713','0714','0715','0720','0721','0722','0723','0726','0727','0729','0730'};

Groot_dir = 'E:\ZC_data\groot_data_switch\prelim_analyses\SFC_CrossChan_Slide';
Groot_tags = {'0415','0507','0509','0513','0514','0515','0516','0520','0521','0522'};

monkey_set = {'Ocean','Groot'};
num_M = length(monkey_set);

for iM = 1:num_M
    curr_M = monkey_set{iM};
    
    switch curr_M
        case {'Ocean'}
            load('G:\FW_data\NewJuly\ocean_data\preprocessed\new_area_codes_July2021.mat','good_chans','area_code')
            [F_chans,F_idx] = intersect(good_chans,setdiff(find(area_code == 2 | area_code == 3),12));
            load('G:\PaperPrep\PostNeuronRevision\Ocean13Days_ROI_idx_concat.mat','ch256_idx','day_idx');
%             fig_dir = sprintf('%s\\figs\\Ocean13Days_SFC_Slide_pT1_pT2_sort',savedir);
%             png_dir = sprintf('%s\\pngs\\Ocean13Days_SFC_Slide_pT1_pT2_sort',savedir);
            fig_dir = sprintf('%s\\figs\\Ocean13Days_SFC_Slide_pT1_pT2_NOsort',savedir);
            png_dir = sprintf('%s\\pngs\\Ocean13Days_SFC_Slide_pT1_pT2_NOsort',savedir);
        case {'Groot'}
            load('E:\ZC_data\groot_data_switch\new_area_codes_Mar2024.mat','good_chans','area_code');
            [F_chans,F_idx] = intersect(good_chans,find(area_code == 2 | area_code == 3));
            load('G:\PaperPrep\PostNeuronRevision\Groot10Days_ROI_idx_concat.mat','ch256_idx','day_idx');
%             fig_dir = sprintf('%s\\figs\\Groot10Days_SFC_Slide_pT1_pT2_sort',savedir);
%             png_dir = sprintf('%s\\pngs\\Groot10Days_SFC_Slide_pT1_pT2_sort',savedir);
            fig_dir = sprintf('%s\\figs\\Groot10Days_SFC_Slide_pT1_pT2_NOsort',savedir);
            png_dir = sprintf('%s\\pngs\\Groot10Days_SFC_Slide_pT1_pT2_NOsort',savedir);
    end
    
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir)
    end
    if ~exist(png_dir,'dir')
        mkdir(png_dir)
    end
            
    eval(sprintf('curr_dir = %s_dir;',curr_M));
    eval(sprintf('num_files = length(%s_tags);',curr_M));
    C_cat = [];
    for iFile = 1:num_files
        eval(sprintf('curr_tag = %s_tags{iFile};',curr_M));
        load(sprintf('%s\\%s%s_T1_Trig_Theta_SFC_trialAve.mat',curr_dir,curr_M,curr_tag),'C_set');
        C_trim = C_set(F_idx,F_idx,:);
        C_mat = squeeze(nansum(C_trim,1));
        [~,F_hits] = intersect(F_chans,ch256_idx(day_idx == iFile));
        C_mat_trim = C_mat(F_hits,:);
        C_cat = [C_cat;C_mat_trim];
        
        h1 = figure; set(h1,'Position',fig_pos);
        
%         [~,T1_idx] = mat_peaksort(C_mat_trim(:,6:16),2);
%         [~,T2_idx] = mat_peaksort(C_mat_trim(:,17:27),2);

%         subplot(1,2,1); imagesc(C_mat_trim(T1_idx,:)); 
        subplot(1,2,1); imagesc(C_mat_trim); 
        colormap hot; colorbar;
        xlim(x_bounds); vline(vline_pos,'w--'); ylabel('Channels');
        set(gca,'XTick',6:5:32);
        set(gca,'XTickLabel',0:250:1300);
        xlabel('Time (ms)')
        title('sort by post-S1-onset 0-500ms peak val times')
        
%         subplot(1,2,2); imagesc(C_mat_trim(T2_idx,:)); 
        subplot(1,2,2); imagesc(C_mat_trim); 
        colorbar; xlim(x_bounds); vline(vline_pos,'w--');
        ylabel('Channels');
        set(gca,'XTick',6:5:32);
        set(gca,'XTickLabel',0:250:1300);
        xlabel('Time (ms)')
        title('sort by post-S2-onset 0-500ms peak val times')
        sgtitle([curr_M,curr_tag,' Theta SFC strength, +-250 ms bins, 50 ms steps, 0-1300 ms from S1 onset'])
        
        saveas(h1,sprintf('%s\\%s%s_SFC_Slide_pT1_pT2_sort.fig',fig_dir,curr_M,curr_tag));
        saveas(h1,sprintf('%s\\%s%s_SFC_Slide_pT1_pT2_sort.png',png_dir,curr_M,curr_tag));
        close(h1);
    end
    
    h2 = figure; set(h2,'Position',fig_pos);
    
%     [~,T1_idx] = mat_peaksort(C_cat(:,6:16),2);
%     [~,T2_idx] = mat_peaksort(C_cat(:,17:27),2);

%     subplot(1,2,1); imagesc(C_cat(T1_idx,:));
    subplot(1,2,1); imagesc(C_cat);
    colormap hot; colorbar;
    xlim(x_bounds); vline(vline_pos,'w--'); ylabel('Channels');
    set(gca,'XTick',6:5:32);
    set(gca,'XTickLabel',0:250:1300);
    xlabel('Time (ms)')
    title('sort by post-S1-onset 0-500ms peak val times')
    
%     subplot(1,2,2); imagesc(C_cat(T2_idx,:));
    subplot(1,2,2); imagesc(C_cat);
    colorbar; xlim(x_bounds); vline(vline_pos,'w--');
    ylabel('Channels');
    set(gca,'XTick',6:5:32);
    set(gca,'XTickLabel',0:250:1300);
    xlabel('Time (ms)')
    title('sort by post-S2-onset 0-500ms peak val times')
    sgtitle([curr_M,' Theta SFC strength, +-250 ms bins, 50 ms steps, 0-1300 ms from S1 onset'])
    
    saveas(h2,sprintf('%s\\%sALLDays_SFC_Slide_pT1_pT2_sort.fig',fig_dir,curr_M));
    saveas(h2,sprintf('%s\\%sALLDays_SFC_Slide_pT1_pT2_sort.png',png_dir,curr_M));
    close(h2);
end

%% Frontal SFC x SPK Mem FW shuffles, 2-Seq, Ocean Switch, ANOVA-AllSig
bin_width_S1 = 0.015; x_range_S1 = [-0.5,0.5];
bin_width_S2 = 0.015; x_range_S2 = [-0.5,0.5];
fig_pos = [505,397,896,565];

loaddir = 'G:\PaperPrep\PostNatNeuroRevision';
load(sprintf('%s\\Ocean13Days_SFCvec_concat_BaseComp_NoSubtract_ZC_AddEnt_TrueIdxVerify_Frontal_NoMask.mat',loaddir))
load(sprintf('%s\\Ocean13Days_FW_concat_srm_ZC_Frontal.mat',loaddir));
load(sprintf('%s\\Ocean13Days_ANOVAvecs_concat_Frontal.mat',loaddir));

ANOVA_sel = ~isnan(ANOVA_vec_SPK_D1);
SFC_T1c_CombZ = SFC_T1c_CombZ(ANOVA_sel);
SFC_T2c_CombZ = SFC_T2c_CombZ(ANOVA_sel);
SPK_S1_D_srm_Comb_ansc = SPK_S1_D_srm_Comb_ansc(ANOVA_sel);
SPK_S2_D_srm_Comb_ansc = SPK_S2_D_srm_Comb_ansc(ANOVA_sel);

orig_mdl_S1 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S1_D_srm_Comb_ansc);
S1_R1_acc_beta = orig_mdl_S1.Coefficients.Estimate(2);
S1_R1_acc_p = orig_mdl_S1.Coefficients.pValue(2);
S1_R2_acc_beta = orig_mdl_S1.Coefficients.Estimate(3);
S1_R2_acc_p = orig_mdl_S1.Coefficients.pValue(3);

orig_mdl_S2 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S2_D_srm_Comb_ansc);
S2_R1_acc_beta = orig_mdl_S2.Coefficients.Estimate(2);
S2_R1_acc_p = orig_mdl_S2.Coefficients.pValue(2);
S2_R2_acc_beta = orig_mdl_S2.Coefficients.Estimate(3);
S2_R2_acc_p = orig_mdl_S2.Coefficients.pValue(3);

num_shuf = 1000;
beta_hold_S1 = zeros(2,num_shuf); % dim1: mod/acc
beta_hold_S2 = zeros(2,num_shuf);

num_chans = length(SPK_S1_D_srm_Comb_ansc);

for iShuf = 1:num_shuf
    rand_idx = randperm(num_chans);
    rand_mdl_S1 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S1_D_srm_Comb_ansc(rand_idx));
    rand_mdl_S2 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S2_D_srm_Comb_ansc(rand_idx));
    
    beta_hold_S1(:,iShuf) = rand_mdl_S1.Coefficients.Estimate(2:3);
    beta_hold_S2(:,iShuf) = rand_mdl_S2.Coefficients.Estimate(2:3);
end

save(sprintf('%s\\SFC_x_FW_stats\\Ocean13Days_ANOVA-AllSig_SPK_Mem_FW_x_SFCacc_perm_stats_2Seq.mat',loaddir),...
    'beta_hold_*','S1_*','S2_*','orig_*')

h = figure; set(h,'Position',fig_pos);
subplot(2,2,1);
histogram(beta_hold_S1(1,:),'BinWidth',bin_width_S1);
vline(S1_R1_acc_beta,'r');
p_perm_R1_S1_acc = sum(S1_R1_acc_beta < beta_hold_S1(1,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ SFC-S1-acc, pP = %.3f, pO = %.3f',p_perm_R1_S1_acc,S1_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,2);
histogram(beta_hold_S1(2,:),'BinWidth',bin_width_S1);
vline(S1_R2_acc_beta,'r');
p_perm_R2_S1_acc = sum(S1_R2_acc_beta < beta_hold_S1(2,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ SFC-S2-acc, pP = %.3f, pO = %.3f',p_perm_R2_S1_acc,S1_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,3);
histogram(beta_hold_S2(1,:),'BinWidth',bin_width_S2);
vline(S2_R1_acc_beta,'r');
p_perm_R1_S2_acc = sum(S2_R1_acc_beta < beta_hold_S2(1,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ SFC-S1-acc, pP = %.3f, pO = %.3f',p_perm_R1_S2_acc,S2_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

subplot(2,2,4);
histogram(beta_hold_S2(2,:),'BinWidth',bin_width_S2);
vline(S2_R2_acc_beta,'r');
p_perm_R2_S2_acc = sum(S2_R2_acc_beta < beta_hold_S2(2,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ SFC-S2-acc, pP = %.3f, pO = %.3f',p_perm_R2_S2_acc,S2_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

saveas(h,sprintf('%s\\figures\\figs\\Ocean13Days_ANOVA-AllSig_SPK_Mem_FW_x_SFCacc_perm_stats_2Seq.fig',loaddir))
saveas(h,sprintf('%s\\figures\\pngs\\Ocean13Days_ANOVA-AllSig_SPK_Mem_FW_x_SFCacc_perm_stats_2Seq.png',loaddir))

%% Frontal SFC x SPK Mem FW shuffles, 2-Seq, Groot Switch, ANOVA-AllSig
bin_width_S1 = 0.015; x_range_S1 = [-0.5,0.5];
bin_width_S2 = 0.015; x_range_S2 = [-0.5,0.5];
fig_pos = [505,397,896,565];

loaddir = 'G:\PaperPrep\PostNatNeuroRevision';
load(sprintf('%s\\GrootSwitch_SFCvec_Frontal_NoMask_11days.mat',loaddir))
load(sprintf('%s\\GrootSwitch_FW_concat_Frontal.mat',loaddir));
load(sprintf('%s\\GrootSwitch_Encoding_SPK_ANOVA_vec.mat',loaddir));

day_sel = 1:824; % 10 days
SFC_T1c_CombZ = SFC_T1c_CombZ(day_sel);
SFC_T2c_CombZ = SFC_T2c_CombZ(day_sel);
SPK_S1_D_Comb_ansc = SPK_S1_D_Comb_ansc(day_sel);
SPK_S2_D_Comb_ansc = SPK_S2_D_Comb_ansc(day_sel);
ANOVA_vec = ANOVA_vec(day_sel);

ANOVA_sel = ~isnan(ANOVA_vec);
SFC_T1c_CombZ = SFC_T1c_CombZ(ANOVA_sel);
SFC_T2c_CombZ = SFC_T2c_CombZ(ANOVA_sel);
SPK_S1_D_Comb_ansc = SPK_S1_D_Comb_ansc(ANOVA_sel);
SPK_S2_D_Comb_ansc = SPK_S2_D_Comb_ansc(ANOVA_sel);

orig_mdl_S1 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S1_D_Comb_ansc);
S1_R1_acc_beta = orig_mdl_S1.Coefficients.Estimate(2);
S1_R1_acc_p = orig_mdl_S1.Coefficients.pValue(2);
S1_R2_acc_beta = orig_mdl_S1.Coefficients.Estimate(3);
S1_R2_acc_p = orig_mdl_S1.Coefficients.pValue(3);

orig_mdl_S2 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S2_D_Comb_ansc);
S2_R1_acc_beta = orig_mdl_S2.Coefficients.Estimate(2);
S2_R1_acc_p = orig_mdl_S2.Coefficients.pValue(2);
S2_R2_acc_beta = orig_mdl_S2.Coefficients.Estimate(3);
S2_R2_acc_p = orig_mdl_S2.Coefficients.pValue(3);

num_shuf = 1000;
beta_hold_S1 = zeros(2,num_shuf); % dim1: mod/acc
beta_hold_S2 = zeros(2,num_shuf);

num_chans = length(SPK_S1_D_Comb_ansc);

for iShuf = 1:num_shuf
    rand_idx = randperm(num_chans);
    rand_mdl_S1 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S1_D_Comb_ansc(rand_idx));
    rand_mdl_S2 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S2_D_Comb_ansc(rand_idx));
    
    beta_hold_S1(:,iShuf) = rand_mdl_S1.Coefficients.Estimate(2:3);
    beta_hold_S2(:,iShuf) = rand_mdl_S2.Coefficients.Estimate(2:3);
end

save(sprintf('%s\\SFC_x_FW_stats\\GrootSwitch10Days_ANOVA-AllSig_SPK_Mem_FW_x_SFCacc_perm_stats_2Seq.mat',loaddir),...
    'beta_hold_*','S1_*','S2_*','orig_*')

h = figure; set(h,'Position',fig_pos);
subplot(2,2,1);
histogram(beta_hold_S1(1,:),'BinWidth',bin_width_S1);
vline(S1_R1_acc_beta,'r');
p_perm_R1_S1_acc = sum(S1_R1_acc_beta < beta_hold_S1(1,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ SFC-S1-acc, pP = %.3f, pO = %.3f',p_perm_R1_S1_acc,S1_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,2);
histogram(beta_hold_S1(2,:),'BinWidth',bin_width_S1);
vline(S1_R2_acc_beta,'r');
p_perm_R2_S1_acc = sum(S1_R2_acc_beta < beta_hold_S1(2,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ SFC-S2-acc, pP = %.3f, pO = %.3f',p_perm_R2_S1_acc,S1_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,3);
histogram(beta_hold_S2(1,:),'BinWidth',bin_width_S2);
vline(S2_R1_acc_beta,'r');
p_perm_R1_S2_acc = sum(S2_R1_acc_beta < beta_hold_S2(1,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ SFC-S1-acc, pP = %.3f, pO = %.3f',p_perm_R1_S2_acc,S2_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

subplot(2,2,4);
histogram(beta_hold_S2(2,:),'BinWidth',bin_width_S2);
vline(S2_R2_acc_beta,'r');
p_perm_R2_S2_acc = sum(S2_R2_acc_beta < beta_hold_S2(2,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ SFC-S2-acc, pP = %.3f, pO = %.3f',p_perm_R2_S2_acc,S2_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

saveas(h,sprintf('%s\\figures\\figs\\GrootSwitch10Days_ANOVA-AllSig_SPK_Mem_FW_x_SFCacc_perm_stats_2Seq.fig',loaddir))
saveas(h,sprintf('%s\\figures\\pngs\\GrootSwitch10Days_ANOVA-AllSig_SPK_Mem_FW_x_SFCacc_perm_stats_2Seq.png',loaddir))

%% Frontal SFC x SPK Mem FW shuffles, 2-Seq, Ocean Switch Err
loaddir = 'G:\PaperPrep\PostNatNeuroRevision';
load(sprintf('%s\\Ocean13Days_SFCvec_concat_BaseComp_NoSubtract_AddEnt_TrueIdxVerify_Frontal_AllErr.mat',loaddir))
SFC_T1c_CombZ_Err = SFC_T1c_CombZ; SFC_T2c_CombZ_Err = SFC_T2c_CombZ;
load(sprintf('%s\\Ocean13Days_SFCvec_concat_BaseComp_NoSubtract_ZC_AddEnt_TrueIdxVerify_Frontal_NoMask.mat',loaddir))
load(sprintf('%s\\Ocean13Days_FW_concat_srm_ZC_Frontal.mat',loaddir));

orig_mdl_S1 = fitlm([SFC_T1c_CombZ_Err,SFC_T2c_CombZ_Err,SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S1_D_srm_Comb_ansc);
S1_R1_acc_Err_beta = orig_mdl_S1.Coefficients.Estimate(2);
S1_R1_acc_Err_p = orig_mdl_S1.Coefficients.pValue(2);
S1_R2_acc_Err_beta = orig_mdl_S1.Coefficients.Estimate(3);
S1_R2_acc_Err_p = orig_mdl_S1.Coefficients.pValue(3);
S1_R1_acc_beta = orig_mdl_S1.Coefficients.Estimate(4);
S1_R1_acc_p = orig_mdl_S1.Coefficients.pValue(4);
S1_R2_acc_beta = orig_mdl_S1.Coefficients.Estimate(5);
S1_R2_acc_p = orig_mdl_S1.Coefficients.pValue(5);

orig_mdl_S2 = fitlm([SFC_T1c_CombZ_Err,SFC_T2c_CombZ_Err,SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S2_D_srm_Comb_ansc);
S2_R1_acc_Err_beta = orig_mdl_S2.Coefficients.Estimate(2);
S2_R1_acc_Err_p = orig_mdl_S2.Coefficients.pValue(2);
S2_R2_acc_Err_beta = orig_mdl_S2.Coefficients.Estimate(3);
S2_R2_acc_Err_p = orig_mdl_S2.Coefficients.pValue(3);
S2_R1_acc_beta = orig_mdl_S2.Coefficients.Estimate(4);
S2_R1_acc_p = orig_mdl_S2.Coefficients.pValue(4);
S2_R2_acc_beta = orig_mdl_S2.Coefficients.Estimate(5);
S2_R2_acc_p = orig_mdl_S2.Coefficients.pValue(5);

num_shuf = 1000;
beta_hold_S1 = zeros(4,num_shuf); % dim1: T1-acc Err, T2-acc Err, T1-acc Corr, T2-acc Corr
beta_hold_S2 = zeros(4,num_shuf);

num_chans = length(SPK_S1_D_srm_Comb_ansc);

for iShuf = 1:num_shuf
    rand_idx = randperm(num_chans);
    rand_mdl_S1 = fitlm([SFC_T1c_CombZ_Err,SFC_T2c_CombZ_Err,SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S1_D_srm_Comb_ansc(rand_idx));
    rand_mdl_S2 = fitlm([SFC_T1c_CombZ_Err,SFC_T2c_CombZ_Err,SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S2_D_srm_Comb_ansc(rand_idx));
    
    beta_hold_S1(:,iShuf) = rand_mdl_S1.Coefficients.Estimate(2:5);
    beta_hold_S2(:,iShuf) = rand_mdl_S2.Coefficients.Estimate(2:5);
end

save(sprintf('%s\\SFC_x_FW_stats\\Ocean13Days_SPK_Mem_FW_x_SFCacc_ErrorAndCorrect_perm_stats_2Seq.mat',loaddir),...
    'beta_hold_*','S1_*','S2_*','orig_*')

%% Frontal SFC x SPK Mem FW shuffles, 2-Seq, Groot Switch Err
loaddir = 'G:\PaperPrep\PostNatNeuroRevision';
day_sel = 1:824;
load(sprintf('%s\\GrootSwitch_SFCvec_Frontal_NoMask_11days_AllErr.mat',loaddir))
SFC_T1c_CombZ_Err = SFC_T1c_CombZ; SFC_T2c_CombZ_Err = SFC_T2c_CombZ;
load(sprintf('%s\\GrootSwitch_SFCvec_Frontal_NoMask_11days.mat',loaddir))
load(sprintf('%s\\GrootSwitch_FW_concat_Frontal.mat',loaddir));

SFC_T1c_CombZ_Err = SFC_T1c_CombZ_Err(day_sel);
SFC_T2c_CombZ_Err = SFC_T2c_CombZ_Err(day_sel);
SFC_T1c_CombZ = SFC_T1c_CombZ(day_sel);
SFC_T2c_CombZ = SFC_T2c_CombZ(day_sel);
SPK_S1_D_Comb_ansc = SPK_S1_D_Comb_ansc(day_sel);
SPK_S2_D_Comb_ansc = SPK_S2_D_Comb_ansc(day_sel);

orig_mdl_S1 = fitlm([SFC_T1c_CombZ_Err,SFC_T2c_CombZ_Err,SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S1_D_Comb_ansc);
S1_R1_acc_Err_beta = orig_mdl_S1.Coefficients.Estimate(2);
S1_R1_acc_Err_p = orig_mdl_S1.Coefficients.pValue(2);
S1_R2_acc_Err_beta = orig_mdl_S1.Coefficients.Estimate(3);
S1_R2_acc_Err_p = orig_mdl_S1.Coefficients.pValue(3);
S1_R1_acc_beta = orig_mdl_S1.Coefficients.Estimate(4);
S1_R1_acc_p = orig_mdl_S1.Coefficients.pValue(4);
S1_R2_acc_beta = orig_mdl_S1.Coefficients.Estimate(5);
S1_R2_acc_p = orig_mdl_S1.Coefficients.pValue(5);

orig_mdl_S2 = fitlm([SFC_T1c_CombZ_Err,SFC_T2c_CombZ_Err,SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S2_D_Comb_ansc);
S2_R1_acc_Err_beta = orig_mdl_S2.Coefficients.Estimate(2);
S2_R1_acc_Err_p = orig_mdl_S2.Coefficients.pValue(2);
S2_R2_acc_Err_beta = orig_mdl_S2.Coefficients.Estimate(3);
S2_R2_acc_Err_p = orig_mdl_S2.Coefficients.pValue(3);
S2_R1_acc_beta = orig_mdl_S2.Coefficients.Estimate(4);
S2_R1_acc_p = orig_mdl_S2.Coefficients.pValue(4);
S2_R2_acc_beta = orig_mdl_S2.Coefficients.Estimate(5);
S2_R2_acc_p = orig_mdl_S2.Coefficients.pValue(5);

num_shuf = 1000;
beta_hold_S1 = zeros(4,num_shuf); % dim1: T1-acc Err, T2-acc Err, T1-acc Corr, T2-acc Corr
beta_hold_S2 = zeros(4,num_shuf);

num_chans = length(SPK_S1_D_Comb_ansc);

for iShuf = 1:num_shuf
    rand_idx = randperm(num_chans);
    rand_mdl_S1 = fitlm([SFC_T1c_CombZ_Err,SFC_T2c_CombZ_Err,SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S1_D_Comb_ansc(rand_idx));
    rand_mdl_S2 = fitlm([SFC_T1c_CombZ_Err,SFC_T2c_CombZ_Err,SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S2_D_Comb_ansc(rand_idx));
    
    beta_hold_S1(:,iShuf) = rand_mdl_S1.Coefficients.Estimate(2:5);
    beta_hold_S2(:,iShuf) = rand_mdl_S2.Coefficients.Estimate(2:5);
end

save(sprintf('%s\\SFC_x_FW_stats\\GrootSwitch10Days_SPK_Mem_FW_x_SFCacc_ErrorAndCorrect_perm_stats_2Seq.mat',loaddir),...
    'beta_hold_*','S1_*','S2_*','orig_*')

%% Frontal SFC x SPK Mem FW shuffles, 2-Seq, Ocean Switch, FR-balanced
bin_width_S1 = 0.015; x_range_S1 = [-0.5,0.5];
bin_width_S2 = 0.015; x_range_S2 = [-0.5,0.5];
fig_pos = [505,397,896,565];

loaddir = 'G:\PaperPrep\PostNatNeuroRevision';
loaddir2 = 'G:\PaperPrep\PostNeuronRevision';
load(sprintf('%s\\Ocean13Days_SFCvec_concat_BaseComp_NoSubtract_ZC_AddEnt_TrueIdxVerify_Frontal_NoMask.mat',loaddir))
load(sprintf('%s\\Ocean13Days_FW_concat_srm_ZC_Frontal.mat',loaddir));
load(sprintf('%s\\Ocean13Days_ANOVAvecs_concat_Frontal.mat',loaddir));
load(sprintf('%s\\TwoMonkeys_SFCacc_FR_diff_vecs.mat',loaddir2),'Ocean_vec');

ANOVA_sel = intersect(find(Ocean_vec==0),find(~isnan(ANOVA_vec_SPK_D1)));
SFC_T1c_CombZ = SFC_T1c_CombZ(ANOVA_sel);
SFC_T2c_CombZ = SFC_T2c_CombZ(ANOVA_sel);
SPK_S1_D_srm_Comb_ansc = SPK_S1_D_srm_Comb_ansc(ANOVA_sel);
SPK_S2_D_srm_Comb_ansc = SPK_S2_D_srm_Comb_ansc(ANOVA_sel);

orig_mdl_S1 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S1_D_srm_Comb_ansc);
S1_R1_acc_beta = orig_mdl_S1.Coefficients.Estimate(2);
S1_R1_acc_p = orig_mdl_S1.Coefficients.pValue(2);
S1_R2_acc_beta = orig_mdl_S1.Coefficients.Estimate(3);
S1_R2_acc_p = orig_mdl_S1.Coefficients.pValue(3);

orig_mdl_S2 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S2_D_srm_Comb_ansc);
S2_R1_acc_beta = orig_mdl_S2.Coefficients.Estimate(2);
S2_R1_acc_p = orig_mdl_S2.Coefficients.pValue(2);
S2_R2_acc_beta = orig_mdl_S2.Coefficients.Estimate(3);
S2_R2_acc_p = orig_mdl_S2.Coefficients.pValue(3);

num_shuf = 1000;
beta_hold_S1 = zeros(2,num_shuf); % dim1: mod/acc
beta_hold_S2 = zeros(2,num_shuf);

num_chans = length(SPK_S1_D_srm_Comb_ansc);

for iShuf = 1:num_shuf
    rand_idx = randperm(num_chans);
    rand_mdl_S1 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S1_D_srm_Comb_ansc(rand_idx));
    rand_mdl_S2 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S2_D_srm_Comb_ansc(rand_idx));
    
    beta_hold_S1(:,iShuf) = rand_mdl_S1.Coefficients.Estimate(2:3);
    beta_hold_S2(:,iShuf) = rand_mdl_S2.Coefficients.Estimate(2:3);
end

save(sprintf('%s\\SFC_x_FW_stats\\Ocean13Days_S1vsS2_FR_NoDiff_SPK_Mem_FW_x_SFCacc_perm_stats_2Seq.mat',loaddir2),...
    'beta_hold_*','S1_*','S2_*','orig_*')

h = figure; set(h,'Position',fig_pos);
subplot(2,2,1);
histogram(beta_hold_S1(1,:),'BinWidth',bin_width_S1);
vline(S1_R1_acc_beta,'r');
p_perm_R1_S1_acc = sum(S1_R1_acc_beta < beta_hold_S1(1,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ SFC-S1-acc, pP = %.3f, pO = %.3f',p_perm_R1_S1_acc,S1_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,2);
histogram(beta_hold_S1(2,:),'BinWidth',bin_width_S1);
vline(S1_R2_acc_beta,'r');
p_perm_R2_S1_acc = sum(S1_R2_acc_beta < beta_hold_S1(2,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ SFC-S2-acc, pP = %.3f, pO = %.3f',p_perm_R2_S1_acc,S1_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,3);
histogram(beta_hold_S2(1,:),'BinWidth',bin_width_S2);
vline(S2_R1_acc_beta,'r');
p_perm_R1_S2_acc = sum(S2_R1_acc_beta < beta_hold_S2(1,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ SFC-S1-acc, pP = %.3f, pO = %.3f',p_perm_R1_S2_acc,S2_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

subplot(2,2,4);
histogram(beta_hold_S2(2,:),'BinWidth',bin_width_S2);
vline(S2_R2_acc_beta,'r');
p_perm_R2_S2_acc = sum(S2_R2_acc_beta < beta_hold_S2(2,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ SFC-S2-acc, pP = %.3f, pO = %.3f',p_perm_R2_S2_acc,S2_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

saveas(h,sprintf('%s\\figures\\figs\\Ocean13Days_S1vsS2_FR_NoDiff_SPK_Mem_FW_x_SFCacc_perm_stats_2Seq.fig',loaddir2))
saveas(h,sprintf('%s\\figures\\pngs\\Ocean13Days_S1vsS2_FR_NoDiff_SPK_Mem_FW_x_SFCacc_perm_stats_2Seq.png',loaddir2))

%% Frontal SFC x SPK Mem FW shuffles, 2-Seq, Groot Switch, FR-balanced
bin_width_S1 = 0.015; x_range_S1 = [-0.5,0.5];
bin_width_S2 = 0.015; x_range_S2 = [-0.5,0.5];
fig_pos = [505,397,896,565];

loaddir = 'G:\PaperPrep\PostNatNeuroRevision';
loaddir2 = 'G:\PaperPrep\PostNeuronRevision';
load(sprintf('%s\\GrootSwitch_SFCvec_Frontal_NoMask_11days.mat',loaddir))
load(sprintf('%s\\GrootSwitch_FW_concat_Frontal.mat',loaddir));
load(sprintf('%s\\GrootSwitch_Encoding_SPK_ANOVA_vec.mat',loaddir));
load(sprintf('%s\\TwoMonkeys_SFCacc_FR_diff_vecs.mat',loaddir2),'Groot_vec');

day_sel = 1:824; % 10 days
SFC_T1c_CombZ = SFC_T1c_CombZ(day_sel);
SFC_T2c_CombZ = SFC_T2c_CombZ(day_sel);
SPK_S1_D_Comb_ansc = SPK_S1_D_Comb_ansc(day_sel);
SPK_S2_D_Comb_ansc = SPK_S2_D_Comb_ansc(day_sel);
ANOVA_vec = ANOVA_vec(day_sel);

ANOVA_sel = intersect(find(Groot_vec == 0),find(~isnan(ANOVA_vec)));
SFC_T1c_CombZ = SFC_T1c_CombZ(ANOVA_sel);
SFC_T2c_CombZ = SFC_T2c_CombZ(ANOVA_sel);
SPK_S1_D_Comb_ansc = SPK_S1_D_Comb_ansc(ANOVA_sel);
SPK_S2_D_Comb_ansc = SPK_S2_D_Comb_ansc(ANOVA_sel);

orig_mdl_S1 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S1_D_Comb_ansc);
S1_R1_acc_beta = orig_mdl_S1.Coefficients.Estimate(2);
S1_R1_acc_p = orig_mdl_S1.Coefficients.pValue(2);
S1_R2_acc_beta = orig_mdl_S1.Coefficients.Estimate(3);
S1_R2_acc_p = orig_mdl_S1.Coefficients.pValue(3);

orig_mdl_S2 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S2_D_Comb_ansc);
S2_R1_acc_beta = orig_mdl_S2.Coefficients.Estimate(2);
S2_R1_acc_p = orig_mdl_S2.Coefficients.pValue(2);
S2_R2_acc_beta = orig_mdl_S2.Coefficients.Estimate(3);
S2_R2_acc_p = orig_mdl_S2.Coefficients.pValue(3);

num_shuf = 1000;
beta_hold_S1 = zeros(2,num_shuf); % dim1: mod/acc
beta_hold_S2 = zeros(2,num_shuf);

num_chans = length(SPK_S1_D_Comb_ansc);

for iShuf = 1:num_shuf
    rand_idx = randperm(num_chans);
    rand_mdl_S1 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S1_D_Comb_ansc(rand_idx));
    rand_mdl_S2 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S2_D_Comb_ansc(rand_idx));
    
    beta_hold_S1(:,iShuf) = rand_mdl_S1.Coefficients.Estimate(2:3);
    beta_hold_S2(:,iShuf) = rand_mdl_S2.Coefficients.Estimate(2:3);
end

save(sprintf('%s\\SFC_x_FW_stats\\GrootSwitch10Days_S1vsS2_FR_NoDiff_SPK_Mem_FW_x_SFCacc_perm_stats_2Seq.mat',loaddir2),...
    'beta_hold_*','S1_*','S2_*','orig_*')

h = figure; set(h,'Position',fig_pos);
subplot(2,2,1);
histogram(beta_hold_S1(1,:),'BinWidth',bin_width_S1);
vline(S1_R1_acc_beta,'r');
p_perm_R1_S1_acc = sum(S1_R1_acc_beta < beta_hold_S1(1,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ SFC-S1-acc, pP = %.3f, pO = %.3f',p_perm_R1_S1_acc,S1_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,2);
histogram(beta_hold_S1(2,:),'BinWidth',bin_width_S1);
vline(S1_R2_acc_beta,'r');
p_perm_R2_S1_acc = sum(S1_R2_acc_beta < beta_hold_S1(2,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ SFC-S2-acc, pP = %.3f, pO = %.3f',p_perm_R2_S1_acc,S1_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,3);
histogram(beta_hold_S2(1,:),'BinWidth',bin_width_S2);
vline(S2_R1_acc_beta,'r');
p_perm_R1_S2_acc = sum(S2_R1_acc_beta < beta_hold_S2(1,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ SFC-S1-acc, pP = %.3f, pO = %.3f',p_perm_R1_S2_acc,S2_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

subplot(2,2,4);
histogram(beta_hold_S2(2,:),'BinWidth',bin_width_S2);
vline(S2_R2_acc_beta,'r');
p_perm_R2_S2_acc = sum(S2_R2_acc_beta < beta_hold_S2(2,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ SFC-S2-acc, pP = %.3f, pO = %.3f',p_perm_R2_S2_acc,S2_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

saveas(h,sprintf('%s\\figures\\figs\\GrootSwitch10Days_S1vsS2_FR_NoDiff_SPK_Mem_FW_x_SFCacc_perm_stats_2Seq.fig',loaddir2))
saveas(h,sprintf('%s\\figures\\pngs\\GrootSwitch10Days_S1vsS2_FR_NoDiff_SPK_Mem_FW_x_SFCacc_perm_stats_2Seq.png',loaddir2))

%% Frontal FR x SPK Mem FW shuffles, 2-Seq, Ocean Switch, ANOVA-AllSig
bin_width_S1 = 0.04; x_range_S1 = [-1.1,1.1];
bin_width_S2 = 0.04; x_range_S2 = [-1.1,1.1];
fig_pos = [505,397,896,565];

loaddir = 'G:\PaperPrep\PostNatNeuroRevision';
load(sprintf('%s\\SpikeRateCtrl\\Ocean13Days_TrialEvent_Spike_Rates_vecs.mat',loaddir))
load(sprintf('%s\\Ocean13Days_FW_concat_srm_ZC_Frontal.mat',loaddir));
load(sprintf('%s\\Ocean13Days_ANOVAvecs_concat_Frontal.mat',loaddir));

ANOVA_sel = ~isnan(ANOVA_vec_SPK_D1);
rate_AZ_vec_pT1_SPK = rate_AZ_vec_pT1_SPK(ANOVA_sel);
rate_AZ_vec_pT2_SPK = rate_AZ_vec_pT2_SPK(ANOVA_sel);
SPK_S1_D_srm_Comb_ansc = SPK_S1_D_srm_Comb_ansc(ANOVA_sel);
SPK_S2_D_srm_Comb_ansc = SPK_S2_D_srm_Comb_ansc(ANOVA_sel);

orig_mdl_S1 = fitlm([rate_AZ_vec_pT1_SPK,rate_AZ_vec_pT2_SPK],SPK_S1_D_srm_Comb_ansc);
S1_R1_acc_beta = orig_mdl_S1.Coefficients.Estimate(2);
S1_R1_acc_p = orig_mdl_S1.Coefficients.pValue(2);
S1_R2_acc_beta = orig_mdl_S1.Coefficients.Estimate(3);
S1_R2_acc_p = orig_mdl_S1.Coefficients.pValue(3);

orig_mdl_S2 = fitlm([rate_AZ_vec_pT1_SPK,rate_AZ_vec_pT2_SPK],SPK_S2_D_srm_Comb_ansc);
S2_R1_acc_beta = orig_mdl_S2.Coefficients.Estimate(2);
S2_R1_acc_p = orig_mdl_S2.Coefficients.pValue(2);
S2_R2_acc_beta = orig_mdl_S2.Coefficients.Estimate(3);
S2_R2_acc_p = orig_mdl_S2.Coefficients.pValue(3);

num_shuf = 1000;
beta_hold_S1 = zeros(2,num_shuf); % dim1: mod/acc
beta_hold_S2 = zeros(2,num_shuf);

num_chans = length(SPK_S1_D_srm_Comb_ansc);

for iShuf = 1:num_shuf
    rand_idx = randperm(num_chans);
    rand_mdl_S1 = fitlm([rate_AZ_vec_pT1_SPK,rate_AZ_vec_pT2_SPK],SPK_S1_D_srm_Comb_ansc(rand_idx));
    rand_mdl_S2 = fitlm([rate_AZ_vec_pT1_SPK,rate_AZ_vec_pT2_SPK],SPK_S2_D_srm_Comb_ansc(rand_idx));
    
    beta_hold_S1(:,iShuf) = rand_mdl_S1.Coefficients.Estimate(2:3);
    beta_hold_S2(:,iShuf) = rand_mdl_S2.Coefficients.Estimate(2:3);
end

save(sprintf('%s\\SFC_x_FW_stats\\Ocean13Days_ANOVA-AllSig_SPK_Mem_FW_x_FRctrl_perm_stats_2Seq.mat',loaddir),...
    'beta_hold_*','S1_*','S2_*','orig_*')

h = figure; set(h,'Position',fig_pos);
subplot(2,2,1);
histogram(beta_hold_S1(1,:),'BinWidth',bin_width_S1);
vline(S1_R1_acc_beta,'r');
p_perm_R1_S1_acc = sum(S1_R1_acc_beta < beta_hold_S1(1,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ FR-S1, pP = %.3f, pO = %.3f',p_perm_R1_S1_acc,S1_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,2);
histogram(beta_hold_S1(2,:),'BinWidth',bin_width_S1);
vline(S1_R2_acc_beta,'r');
p_perm_R2_S1_acc = sum(S1_R2_acc_beta < beta_hold_S1(2,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ FR-S2, pP = %.3f, pO = %.3f',p_perm_R2_S1_acc,S1_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,3);
histogram(beta_hold_S2(1,:),'BinWidth',bin_width_S2);
vline(S2_R1_acc_beta,'r');
p_perm_R1_S2_acc = sum(S2_R1_acc_beta < beta_hold_S2(1,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ FR-S1, pP = %.3f, pO = %.3f',p_perm_R1_S2_acc,S2_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

subplot(2,2,4);
histogram(beta_hold_S2(2,:),'BinWidth',bin_width_S2);
vline(S2_R2_acc_beta,'r');
p_perm_R2_S2_acc = sum(S2_R2_acc_beta < beta_hold_S2(2,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ FR-S2, pP = %.3f, pO = %.3f',p_perm_R2_S2_acc,S2_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

saveas(h,sprintf('%s\\figures\\figs\\Ocean13Days_ANOVA-AllSig_SPK_Mem_FW_x_FRctrl_perm_stats_2Seq.fig',loaddir))
saveas(h,sprintf('%s\\figures\\pngs\\Ocean13Days_ANOVA-AllSig_SPK_Mem_FW_x_FRctrl_perm_stats_2Seq.png',loaddir))

%% Frontal FR x SPK Mem FW shuffles, 2-Seq, Groot Switch, ANOVA-AllSig
bin_width_S1 = 0.025; x_range_S1 = [-0.9,0.9];
bin_width_S2 = 0.025; x_range_S2 = [-0.9,0.9];
fig_pos = [505,397,896,565];

loaddir = 'G:\PaperPrep\PostNatNeuroRevision';
load(sprintf('%s\\SpikeRateCtrl\\GrootSwitch11Days_TrialEvent_Spike_Rates_vecs.mat',loaddir))
load(sprintf('%s\\GrootSwitch_FW_concat_Frontal.mat',loaddir));
load(sprintf('%s\\GrootSwitch_Encoding_SPK_ANOVA_vec.mat',loaddir));

day_sel = 1:824; % 10 days
rate_AZ_vec_pT1_SPK = rate_AZ_vec_pT1_SPK(day_sel);
rate_AZ_vec_pT2_SPK = rate_AZ_vec_pT2_SPK(day_sel);
SPK_S1_D_Comb_ansc = SPK_S1_D_Comb_ansc(day_sel);
SPK_S2_D_Comb_ansc = SPK_S2_D_Comb_ansc(day_sel);
ANOVA_vec = ANOVA_vec(day_sel);

ANOVA_sel = ~isnan(ANOVA_vec);
rate_AZ_vec_pT1_SPK = rate_AZ_vec_pT1_SPK(ANOVA_sel);
rate_AZ_vec_pT2_SPK = rate_AZ_vec_pT2_SPK(ANOVA_sel);
SPK_S1_D_Comb_ansc = SPK_S1_D_Comb_ansc(ANOVA_sel);
SPK_S2_D_Comb_ansc = SPK_S2_D_Comb_ansc(ANOVA_sel);

orig_mdl_S1 = fitlm([rate_AZ_vec_pT1_SPK,rate_AZ_vec_pT2_SPK],SPK_S1_D_Comb_ansc);
S1_R1_acc_beta = orig_mdl_S1.Coefficients.Estimate(2);
S1_R1_acc_p = orig_mdl_S1.Coefficients.pValue(2);
S1_R2_acc_beta = orig_mdl_S1.Coefficients.Estimate(3);
S1_R2_acc_p = orig_mdl_S1.Coefficients.pValue(3);

orig_mdl_S2 = fitlm([rate_AZ_vec_pT1_SPK,rate_AZ_vec_pT2_SPK],SPK_S2_D_Comb_ansc);
S2_R1_acc_beta = orig_mdl_S2.Coefficients.Estimate(2);
S2_R1_acc_p = orig_mdl_S2.Coefficients.pValue(2);
S2_R2_acc_beta = orig_mdl_S2.Coefficients.Estimate(3);
S2_R2_acc_p = orig_mdl_S2.Coefficients.pValue(3);

num_shuf = 1000;
beta_hold_S1 = zeros(2,num_shuf); % dim1: mod/acc
beta_hold_S2 = zeros(2,num_shuf);

num_chans = length(SPK_S1_D_Comb_ansc);

for iShuf = 1:num_shuf
    rand_idx = randperm(num_chans);
    rand_mdl_S1 = fitlm([rate_AZ_vec_pT1_SPK,rate_AZ_vec_pT2_SPK],SPK_S1_D_Comb_ansc(rand_idx));
    rand_mdl_S2 = fitlm([rate_AZ_vec_pT1_SPK,rate_AZ_vec_pT2_SPK],SPK_S2_D_Comb_ansc(rand_idx));
    
    beta_hold_S1(:,iShuf) = rand_mdl_S1.Coefficients.Estimate(2:3);
    beta_hold_S2(:,iShuf) = rand_mdl_S2.Coefficients.Estimate(2:3);
end

save(sprintf('%s\\SFC_x_FW_stats\\GrootSwitch10Days_ANOVA-AllSig_SPK_Mem_FW_x_FRctrl_perm_stats_2Seq.mat',loaddir),...
    'beta_hold_*','S1_*','S2_*','orig_*')

h = figure; set(h,'Position',fig_pos);
subplot(2,2,1);
histogram(beta_hold_S1(1,:),'BinWidth',bin_width_S1);
vline(S1_R1_acc_beta,'r');
p_perm_R1_S1_acc = sum(S1_R1_acc_beta < beta_hold_S1(1,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ FR-S1, pP = %.3f, pO = %.3f',p_perm_R1_S1_acc,S1_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,2);
histogram(beta_hold_S1(2,:),'BinWidth',bin_width_S1);
vline(S1_R2_acc_beta,'r');
p_perm_R2_S1_acc = sum(S1_R2_acc_beta < beta_hold_S1(2,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ FR-S2, pP = %.3f, pO = %.3f',p_perm_R2_S1_acc,S1_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,3);
histogram(beta_hold_S2(1,:),'BinWidth',bin_width_S2);
vline(S2_R1_acc_beta,'r');
p_perm_R1_S2_acc = sum(S2_R1_acc_beta < beta_hold_S2(1,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ FR-S1, pP = %.3f, pO = %.3f',p_perm_R1_S2_acc,S2_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

subplot(2,2,4);
histogram(beta_hold_S2(2,:),'BinWidth',bin_width_S2);
vline(S2_R2_acc_beta,'r');
p_perm_R2_S2_acc = sum(S2_R2_acc_beta < beta_hold_S2(2,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ FR-S2, pP = %.3f, pO = %.3f',p_perm_R2_S2_acc,S2_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

saveas(h,sprintf('%s\\figures\\figs\\GrootSwitch10Days_ANOVA-AllSig_SPK_Mem_FW_x_FRctrl_perm_stats_2Seq.fig',loaddir))
saveas(h,sprintf('%s\\figures\\pngs\\GrootSwitch10Days_ANOVA-AllSig_SPK_Mem_FW_x_FRctrl_perm_stats_2Seq.png',loaddir))

%% Frontal TH-ERSP x SPK Mem FW shuffles, 2-Seq, Ocean Switch, ANOVA-AllSig
bin_width_S1 = 0.015; x_range_S1 = [-0.5,0.5];
bin_width_S2 = 0.015; x_range_S2 = [-0.5,0.5];
fig_pos = [505,397,896,565];

loaddir = 'G:\PaperPrep\PostNatNeuroRevision';
load(sprintf('%s\\Ocean13Days_FW_concat_srm_ZC_Frontal.mat',loaddir));
load(sprintf('%s\\Ocean13Days_ANOVAvecs_concat_Frontal.mat',loaddir));

loaddir2 = 'G:\PaperPrep\PostNeuronRevision';
load(sprintf('%s\\Ocean13Days_TH_ERSP_daynorm_vecs.mat',loaddir2))

ANOVA_sel = ~isnan(ANOVA_vec_SPK_D1);
TH_T1_CombZ = TH_pT1_O_norm(ANOVA_sel);
TH_T2_CombZ = TH_pT2_O_norm(ANOVA_sel);
SPK_S1_D_srm_Comb_ansc = SPK_S1_D_srm_Comb_ansc(ANOVA_sel);
SPK_S2_D_srm_Comb_ansc = SPK_S2_D_srm_Comb_ansc(ANOVA_sel);

orig_mdl_S1 = fitlm([TH_T1_CombZ,TH_T2_CombZ],SPK_S1_D_srm_Comb_ansc);
S1_R1_acc_beta = orig_mdl_S1.Coefficients.Estimate(2);
S1_R1_acc_p = orig_mdl_S1.Coefficients.pValue(2);
S1_R2_acc_beta = orig_mdl_S1.Coefficients.Estimate(3);
S1_R2_acc_p = orig_mdl_S1.Coefficients.pValue(3);

orig_mdl_S2 = fitlm([TH_T1_CombZ,TH_T2_CombZ],SPK_S2_D_srm_Comb_ansc);
S2_R1_acc_beta = orig_mdl_S2.Coefficients.Estimate(2);
S2_R1_acc_p = orig_mdl_S2.Coefficients.pValue(2);
S2_R2_acc_beta = orig_mdl_S2.Coefficients.Estimate(3);
S2_R2_acc_p = orig_mdl_S2.Coefficients.pValue(3);

num_shuf = 1000;
beta_hold_S1 = zeros(2,num_shuf); % dim1: mod/acc
beta_hold_S2 = zeros(2,num_shuf);

num_chans = length(SPK_S1_D_srm_Comb_ansc);

for iShuf = 1:num_shuf
    rand_idx = randperm(num_chans);
    rand_mdl_S1 = fitlm([TH_T1_CombZ,TH_T2_CombZ],SPK_S1_D_srm_Comb_ansc(rand_idx));
    rand_mdl_S2 = fitlm([TH_T1_CombZ,TH_T2_CombZ],SPK_S2_D_srm_Comb_ansc(rand_idx));
    
    beta_hold_S1(:,iShuf) = rand_mdl_S1.Coefficients.Estimate(2:3);
    beta_hold_S2(:,iShuf) = rand_mdl_S2.Coefficients.Estimate(2:3);
end

save(sprintf('%s\\TH_x_FW_stats\\Ocean13Days_ANOVA-AllSig_SPK_Mem_FW_x_TH-ERSP_perm_stats_2Seq.mat',loaddir2),...
    'beta_hold_*','S1_*','S2_*','orig_*')

h = figure; set(h,'Position',fig_pos);
subplot(2,2,1);
histogram(beta_hold_S1(1,:),'BinWidth',bin_width_S1);
vline(S1_R1_acc_beta,'r');
p_perm_R1_S1_acc = sum(S1_R1_acc_beta < beta_hold_S1(1,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ TH-S1-ERSP, pP = %.3f, pO = %.3f',p_perm_R1_S1_acc,S1_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,2);
histogram(beta_hold_S1(2,:),'BinWidth',bin_width_S1);
vline(S1_R2_acc_beta,'r');
p_perm_R2_S1_acc = sum(S1_R2_acc_beta < beta_hold_S1(2,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ TH-S2-ERSP, pP = %.3f, pO = %.3f',p_perm_R2_S1_acc,S1_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,3);
histogram(beta_hold_S2(1,:),'BinWidth',bin_width_S2);
vline(S2_R1_acc_beta,'r');
p_perm_R1_S2_acc = sum(S2_R1_acc_beta < beta_hold_S2(1,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ TH-S1-ERSP, pP = %.3f, pO = %.3f',p_perm_R1_S2_acc,S2_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

subplot(2,2,4);
histogram(beta_hold_S2(2,:),'BinWidth',bin_width_S2);
vline(S2_R2_acc_beta,'r');
p_perm_R2_S2_acc = sum(S2_R2_acc_beta < beta_hold_S2(2,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ TH-S2-ERSP, pP = %.3f, pO = %.3f',p_perm_R2_S2_acc,S2_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

saveas(h,sprintf('%s\\figures\\figs\\Ocean13Days_ANOVA-AllSig_SPK_Mem_FW_x_TH-ERSP_perm_stats_2Seq.fig',loaddir2))
saveas(h,sprintf('%s\\figures\\pngs\\Ocean13Days_ANOVA-AllSig_SPK_Mem_FW_x_TH-ERSP_perm_stats_2Seq.png',loaddir2))

%% Frontal TH-ERSP x SPK Mem FW shuffles, 2-Seq, Groot Switch, ANOVA-AllSig
bin_width_S1 = 0.015; x_range_S1 = [-0.5,0.5];
bin_width_S2 = 0.015; x_range_S2 = [-0.5,0.5];
fig_pos = [505,397,896,565];

loaddir = 'G:\PaperPrep\PostNatNeuroRevision';
load(sprintf('%s\\GrootSwitch_FW_concat_Frontal.mat',loaddir));
load(sprintf('%s\\GrootSwitch_Encoding_SPK_ANOVA_vec.mat',loaddir));

loaddir2 = 'G:\PaperPrep\PostNeuronRevision';
load(sprintf('%s\\Groot10Days_TH_ERSP_daynorm_vecs.mat',loaddir2))

day_sel = 1:824; % 10 days
SPK_S1_D_Comb_ansc = SPK_S1_D_Comb_ansc(day_sel);
SPK_S2_D_Comb_ansc = SPK_S2_D_Comb_ansc(day_sel);
ANOVA_vec = ANOVA_vec(day_sel);

ANOVA_sel = ~isnan(ANOVA_vec);
TH_T1_CombZ = TH_pT1_G_norm(ANOVA_sel);
TH_T2_CombZ = TH_pT2_G_norm(ANOVA_sel);
SPK_S1_D_Comb_ansc = SPK_S1_D_Comb_ansc(ANOVA_sel);
SPK_S2_D_Comb_ansc = SPK_S2_D_Comb_ansc(ANOVA_sel);

orig_mdl_S1 = fitlm([TH_T1_CombZ,TH_T2_CombZ],SPK_S1_D_Comb_ansc);
S1_R1_acc_beta = orig_mdl_S1.Coefficients.Estimate(2);
S1_R1_acc_p = orig_mdl_S1.Coefficients.pValue(2);
S1_R2_acc_beta = orig_mdl_S1.Coefficients.Estimate(3);
S1_R2_acc_p = orig_mdl_S1.Coefficients.pValue(3);

orig_mdl_S2 = fitlm([TH_T1_CombZ,TH_T2_CombZ],SPK_S2_D_Comb_ansc);
S2_R1_acc_beta = orig_mdl_S2.Coefficients.Estimate(2);
S2_R1_acc_p = orig_mdl_S2.Coefficients.pValue(2);
S2_R2_acc_beta = orig_mdl_S2.Coefficients.Estimate(3);
S2_R2_acc_p = orig_mdl_S2.Coefficients.pValue(3);

num_shuf = 1000;
beta_hold_S1 = zeros(2,num_shuf); % dim1: mod/acc
beta_hold_S2 = zeros(2,num_shuf);

num_chans = length(SPK_S1_D_Comb_ansc);

for iShuf = 1:num_shuf
    rand_idx = randperm(num_chans);
    rand_mdl_S1 = fitlm([TH_T1_CombZ,TH_T2_CombZ],SPK_S1_D_Comb_ansc(rand_idx));
    rand_mdl_S2 = fitlm([TH_T1_CombZ,TH_T2_CombZ],SPK_S2_D_Comb_ansc(rand_idx));
    
    beta_hold_S1(:,iShuf) = rand_mdl_S1.Coefficients.Estimate(2:3);
    beta_hold_S2(:,iShuf) = rand_mdl_S2.Coefficients.Estimate(2:3);
end

save(sprintf('%s\\TH_x_FW_stats\\GrootSwitch10Days_ANOVA-AllSig_SPK_Mem_FW_x_TH-ERSP_perm_stats_2Seq.mat',loaddir2),...
    'beta_hold_*','S1_*','S2_*','orig_*')

h = figure; set(h,'Position',fig_pos);
subplot(2,2,1);
histogram(beta_hold_S1(1,:),'BinWidth',bin_width_S1);
vline(S1_R1_acc_beta,'r');
p_perm_R1_S1_acc = sum(S1_R1_acc_beta < beta_hold_S1(1,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ TH-S1-ERSP, pP = %.3f, pO = %.3f',p_perm_R1_S1_acc,S1_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,2);
histogram(beta_hold_S1(2,:),'BinWidth',bin_width_S1);
vline(S1_R2_acc_beta,'r');
p_perm_R2_S1_acc = sum(S1_R2_acc_beta < beta_hold_S1(2,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ TH-S2-ERSP, pP = %.3f, pO = %.3f',p_perm_R2_S1_acc,S1_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,3);
histogram(beta_hold_S2(1,:),'BinWidth',bin_width_S2);
vline(S2_R1_acc_beta,'r');
p_perm_R1_S2_acc = sum(S2_R1_acc_beta < beta_hold_S2(1,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ TH-S1-ERSP, pP = %.3f, pO = %.3f',p_perm_R1_S2_acc,S2_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

subplot(2,2,4);
histogram(beta_hold_S2(2,:),'BinWidth',bin_width_S2);
vline(S2_R2_acc_beta,'r');
p_perm_R2_S2_acc = sum(S2_R2_acc_beta < beta_hold_S2(2,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ TH-S2-ERSP, pP = %.3f, pO = %.3f',p_perm_R2_S2_acc,S2_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

saveas(h,sprintf('%s\\figures\\figs\\GrootSwitch10Days_ANOVA-AllSig_SPK_Mem_FW_x_TH-ERSP_perm_stats_2Seq.fig',loaddir2))
saveas(h,sprintf('%s\\figures\\pngs\\GrootSwitch10Days_ANOVA-AllSig_SPK_Mem_FW_x_TH-ERSP_perm_stats_2Seq.png',loaddir2))

%% Frontal Beta-SFC x SPK Mem FW shuffles, 2-Seq, Ocean Switch, ANOVA-AllSig
bin_width_S1 = 0.015; x_range_S1 = [-0.5,0.5];
bin_width_S2 = 0.015; x_range_S2 = [-0.5,0.5];
fig_pos = [505,397,896,565];

loaddir = 'G:\PaperPrep\CurrBiolRevision';
load(sprintf('%s\\Ocean13Days_FW_concat_srm_ZC_Frontal.mat',loaddir));
load(sprintf('%s\\Ocean13Days_ANOVAvec_concat_Frontal.mat',loaddir));

loaddir2 = 'G:\PaperPrep\CurrBiolRevision';
load(sprintf('%s\\Ocean13Days_BetaSFCvec_concat.mat',loaddir2))

ANOVA_sel = ~isnan(ANOVA_vec_SPK_D1);
BT_T1_CombZ = SFC_T1c_CombZ(ANOVA_sel);
BT_T2_CombZ = SFC_T2c_CombZ(ANOVA_sel);
SPK_S1_D_srm_Comb_ansc = SPK_S1_D_srm_Comb_ansc(ANOVA_sel);
SPK_S2_D_srm_Comb_ansc = SPK_S2_D_srm_Comb_ansc(ANOVA_sel);

orig_mdl_S1 = fitlm([BT_T1_CombZ,BT_T2_CombZ],SPK_S1_D_srm_Comb_ansc);
S1_R1_acc_beta = orig_mdl_S1.Coefficients.Estimate(2);
S1_R1_acc_p = orig_mdl_S1.Coefficients.pValue(2);
S1_R2_acc_beta = orig_mdl_S1.Coefficients.Estimate(3);
S1_R2_acc_p = orig_mdl_S1.Coefficients.pValue(3);

orig_mdl_S2 = fitlm([BT_T1_CombZ,BT_T2_CombZ],SPK_S2_D_srm_Comb_ansc);
S2_R1_acc_beta = orig_mdl_S2.Coefficients.Estimate(2);
S2_R1_acc_p = orig_mdl_S2.Coefficients.pValue(2);
S2_R2_acc_beta = orig_mdl_S2.Coefficients.Estimate(3);
S2_R2_acc_p = orig_mdl_S2.Coefficients.pValue(3);

num_shuf = 1000;
beta_hold_S1 = zeros(2,num_shuf); % dim1: mod/acc
beta_hold_S2 = zeros(2,num_shuf);

num_chans = length(SPK_S1_D_srm_Comb_ansc);

for iShuf = 1:num_shuf
    rand_idx = randperm(num_chans);
    rand_mdl_S1 = fitlm([BT_T1_CombZ,BT_T2_CombZ],SPK_S1_D_srm_Comb_ansc(rand_idx));
    rand_mdl_S2 = fitlm([BT_T1_CombZ,BT_T2_CombZ],SPK_S2_D_srm_Comb_ansc(rand_idx));
    
    beta_hold_S1(:,iShuf) = rand_mdl_S1.Coefficients.Estimate(2:3);
    beta_hold_S2(:,iShuf) = rand_mdl_S2.Coefficients.Estimate(2:3);
end

save(sprintf('%s\\BT_SFC_x_FW_stats\\Ocean13Days_ANOVA-AllSig_SPK_Mem_FW_x_BT-SFC_perm_stats_2Seq.mat',loaddir2),...
    'beta_hold_*','S1_*','S2_*','orig_*')

h = figure; set(h,'Position',fig_pos);
subplot(2,2,1);
histogram(beta_hold_S1(1,:),'BinWidth',bin_width_S1);
vline(S1_R1_acc_beta,'r');
p_perm_R1_S1_acc = sum(S1_R1_acc_beta < beta_hold_S1(1,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ BT-S1-SFC, pP = %.3f, pO = %.3f',p_perm_R1_S1_acc,S1_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,2);
histogram(beta_hold_S1(2,:),'BinWidth',bin_width_S1);
vline(S1_R2_acc_beta,'r');
p_perm_R2_S1_acc = sum(S1_R2_acc_beta < beta_hold_S1(2,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ BT-S2-SFC, pP = %.3f, pO = %.3f',p_perm_R2_S1_acc,S1_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,3);
histogram(beta_hold_S2(1,:),'BinWidth',bin_width_S2);
vline(S2_R1_acc_beta,'r');
p_perm_R1_S2_acc = sum(S2_R1_acc_beta < beta_hold_S2(1,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ BT-S1-SFC, pP = %.3f, pO = %.3f',p_perm_R1_S2_acc,S2_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

subplot(2,2,4);
histogram(beta_hold_S2(2,:),'BinWidth',bin_width_S2);
vline(S2_R2_acc_beta,'r');
p_perm_R2_S2_acc = sum(S2_R2_acc_beta < beta_hold_S2(2,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ BT-S2-SFC, pP = %.3f, pO = %.3f',p_perm_R2_S2_acc,S2_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

saveas(h,sprintf('%s\\figs\\Ocean13Days_ANOVA-AllSig_SPK_Mem_FW_x_BT-SFC_perm_stats_2Seq.fig',loaddir2))
saveas(h,sprintf('%s\\pngs\\Ocean13Days_ANOVA-AllSig_SPK_Mem_FW_x_BT-SFC_perm_stats_2Seq.png',loaddir2))

%% Frontal Beta-SFC x SPK Mem FW shuffles, 2-Seq, Groot Switch, ANOVA-AllSig
bin_width_S1 = 0.015; x_range_S1 = [-0.5,0.5];
bin_width_S2 = 0.015; x_range_S2 = [-0.5,0.5];
fig_pos = [505,397,896,565];

loaddir = 'G:\PaperPrep\CurrBiolRevision';
load(sprintf('%s\\GrootSwitch_FW_concat_Frontal.mat',loaddir));
load(sprintf('%s\\Groot10Days_ANOVAvec_concat_Frontal.mat',loaddir));

loaddir2 = 'G:\PaperPrep\CurrBiolRevision';
load(sprintf('%s\\Groot10Days_BetaSFCvec_concat.mat',loaddir2))

day_sel = 1:824; % 10 days
SPK_S1_D_Comb_ansc = SPK_S1_D_Comb_ansc(day_sel);
SPK_S2_D_Comb_ansc = SPK_S2_D_Comb_ansc(day_sel);
ANOVA_vec = ANOVA_vec_SPK_D1(day_sel);

ANOVA_sel = ~isnan(ANOVA_vec);
BT_T1_CombZ = SFC_T1c_CombZ(ANOVA_sel);
BT_T2_CombZ = SFC_T2c_CombZ(ANOVA_sel);
SPK_S1_D_Comb_ansc = SPK_S1_D_Comb_ansc(ANOVA_sel);
SPK_S2_D_Comb_ansc = SPK_S2_D_Comb_ansc(ANOVA_sel);

orig_mdl_S1 = fitlm([BT_T1_CombZ,BT_T2_CombZ],SPK_S1_D_Comb_ansc);
S1_R1_acc_beta = orig_mdl_S1.Coefficients.Estimate(2);
S1_R1_acc_p = orig_mdl_S1.Coefficients.pValue(2);
S1_R2_acc_beta = orig_mdl_S1.Coefficients.Estimate(3);
S1_R2_acc_p = orig_mdl_S1.Coefficients.pValue(3);

orig_mdl_S2 = fitlm([BT_T1_CombZ,BT_T2_CombZ],SPK_S2_D_Comb_ansc);
S2_R1_acc_beta = orig_mdl_S2.Coefficients.Estimate(2);
S2_R1_acc_p = orig_mdl_S2.Coefficients.pValue(2);
S2_R2_acc_beta = orig_mdl_S2.Coefficients.Estimate(3);
S2_R2_acc_p = orig_mdl_S2.Coefficients.pValue(3);

num_shuf = 1000;
beta_hold_S1 = zeros(2,num_shuf); % dim1: mod/acc
beta_hold_S2 = zeros(2,num_shuf);

num_chans = length(SPK_S1_D_Comb_ansc);

for iShuf = 1:num_shuf
    rand_idx = randperm(num_chans);
    rand_mdl_S1 = fitlm([BT_T1_CombZ,BT_T2_CombZ],SPK_S1_D_Comb_ansc(rand_idx));
    rand_mdl_S2 = fitlm([BT_T1_CombZ,BT_T2_CombZ],SPK_S2_D_Comb_ansc(rand_idx));
    
    beta_hold_S1(:,iShuf) = rand_mdl_S1.Coefficients.Estimate(2:3);
    beta_hold_S2(:,iShuf) = rand_mdl_S2.Coefficients.Estimate(2:3);
end

save(sprintf('%s\\BT_SFC_x_FW_stats\\GrootSwitch10Days_ANOVA-AllSig_SPK_Mem_FW_x_BT-SFC_perm_stats_2Seq.mat',loaddir2),...
    'beta_hold_*','S1_*','S2_*','orig_*')

h = figure; set(h,'Position',fig_pos);
subplot(2,2,1);
histogram(beta_hold_S1(1,:),'BinWidth',bin_width_S1);
vline(S1_R1_acc_beta,'r');
p_perm_R1_S1_acc = sum(S1_R1_acc_beta < beta_hold_S1(1,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ BT-S1-SFC, pP = %.3f, pO = %.3f',p_perm_R1_S1_acc,S1_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,2);
histogram(beta_hold_S1(2,:),'BinWidth',bin_width_S1);
vline(S1_R2_acc_beta,'r');
p_perm_R2_S1_acc = sum(S1_R2_acc_beta < beta_hold_S1(2,:))/num_shuf;
title(sprintf('SPK S1 Mem FW ~ BT-S2-SFC, pP = %.3f, pO = %.3f',p_perm_R2_S1_acc,S1_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S1);

subplot(2,2,3);
histogram(beta_hold_S2(1,:),'BinWidth',bin_width_S2);
vline(S2_R1_acc_beta,'r');
p_perm_R1_S2_acc = sum(S2_R1_acc_beta < beta_hold_S2(1,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ BT-S1-SFC, pP = %.3f, pO = %.3f',p_perm_R1_S2_acc,S2_R1_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

subplot(2,2,4);
histogram(beta_hold_S2(2,:),'BinWidth',bin_width_S2);
vline(S2_R2_acc_beta,'r');
p_perm_R2_S2_acc = sum(S2_R2_acc_beta < beta_hold_S2(2,:))/num_shuf;
title(sprintf('SPK S2 Mem FW ~ BT-S2-SFC, pP = %.3f, pO = %.3f',p_perm_R2_S2_acc,S2_R2_acc_p))
xlabel('LM Beta'); ylabel('# Perms');
xlim(x_range_S2);

saveas(h,sprintf('%s\\figs\\GrootSwitch10Days_ANOVA-AllSig_SPK_Mem_FW_x_BT-SFC_perm_stats_2Seq.fig',loaddir2))
saveas(h,sprintf('%s\\pngs\\GrootSwitch10Days_ANOVA-AllSig_SPK_Mem_FW_x_BT-SFC_perm_stats_2Seq.png',loaddir2))

%% SFC phase distrib. plots, Mem1 vs. Mem2 in pT1 vs. pT2
EDGES = 0:(2*pi/30):2*pi;
h = figure;
subplot(1,2,1);
[N1,~] = histcounts(phase_T1_S1_O,'BinEdges',EDGES);
[N2,~] = histcounts(phase_T2_S2_O,'BinEdges',EDGES);
bar(N1./sum(N1),0.9,'FaceAlpha',0.5,'LineStyle','none');
hold on; bar(N2./sum(N2),0.9,'FaceAlpha',0.5,'LineStyle','none');
ylim([0.015,0.055]);
title('Ocean, Mem1 vs. Mem2')

subplot(1,2,2);
[N1,~] = histcounts(phase_T1_S1_G,'BinEdges',EDGES);
[N2,~] = histcounts(phase_T2_S2_G,'BinEdges',EDGES);
bar(N1./sum(N1),0.9,'FaceAlpha',0.5,'LineStyle','none');
hold on; bar(N2./sum(N2),0.9,'FaceAlpha',0.5,'LineStyle','none');
ylim([0.015,0.055]);
title('Groot, Mem1 vs. Mem2')
