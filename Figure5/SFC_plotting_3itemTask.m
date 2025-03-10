%% Frontal SFC x SPK Mem FW shuffles, 3-Seq, Ocean, 3x predictors stepLM, ANOVA-AllSig
loaddir = 'G:\PaperPrep\PostNatNeuroRevision';

load(sprintf('%s\\Ocean3Seq_FW_concat_NewFW2024_srm_ZC_Frontal_321_adj.mat',loaddir));
load(sprintf('%s\\Ocean3Seq_SFCvec_concat_BaseComp_NoSubtract_AddEnt_TrueIdxVerify_Frontal_NoMask.mat',loaddir));
load(sprintf('%s\\Ocean3Seq_SPK_ANOVA_concat_Frontal.mat',loaddir));
ANOVA_sel = find(~isnan(ANOVA_vec_R3));
SFC_T1c_CombZ = SFC_T1c_CombZ(ANOVA_sel);
SFC_T2c_CombZ = SFC_T2c_CombZ(ANOVA_sel);
SFC_T3c_CombZ = SFC_T3c_CombZ(ANOVA_sel);
SPK_S1_D_Comb_ansc = SPK_S1_D_Comb_ansc(ANOVA_sel);
SPK_S2_D_Comb_ansc = SPK_S2_D_Comb_ansc(ANOVA_sel);
SPK_S3_D_Comb_ansc = SPK_S3_D_Comb_ansc(ANOVA_sel);

tbl = table(SFC_T1c_CombZ,SFC_T2c_CombZ,SFC_T3c_CombZ,...
    SPK_S1_D_Comb_ansc,SPK_S2_D_Comb_ansc,SPK_S3_D_Comb_ansc,...
    'VariableNames',{'SFC_T1acc','SFC_T2acc','SFC_T3acc',...
    'SPK_S1_D','SPK_S2_D','SPK_S3_D'});

% S1 model: SPK_S1_D ~ SFC_T1acc + SFC_T3acc
orig_mdl_S1 = stepwiselm(tbl,'ResponseVar','SPK_S1_D',...
    'PredictorVars',{'SFC_T1acc','SFC_T2acc','SFC_T3acc'},'Upper','linear');
S1_R1_acc_beta = orig_mdl_S1.Coefficients.Estimate(2);
S1_R1_acc_p = orig_mdl_S1.Coefficients.pValue(2);
S1_R3_acc_beta = orig_mdl_S1.Coefficients.Estimate(3);
S1_R3_acc_p = orig_mdl_S1.Coefficients.pValue(3);

% S2 model: SPK_S2_D ~ SFC_T1acc + SFC_T2acc
orig_mdl_S2 = stepwiselm(tbl,'ResponseVar','SPK_S2_D',...
    'PredictorVars',{'SFC_T1acc','SFC_T2acc','SFC_T3acc'},'Upper','linear');
S2_R1_acc_beta = orig_mdl_S2.Coefficients.Estimate(2);
S2_R1_acc_p = orig_mdl_S2.Coefficients.pValue(2);
S2_R2_acc_beta = orig_mdl_S2.Coefficients.Estimate(3);
S2_R2_acc_p = orig_mdl_S2.Coefficients.pValue(3);

% S3 model: SPK_S3_D ~ SFC_T3acc
orig_mdl_S3 = stepwiselm(tbl,'ResponseVar','SPK_S3_D',...
    'PredictorVars',{'SFC_T1acc','SFC_T2acc','SFC_T3acc'},'Upper','linear');
S3_R3_acc_beta = orig_mdl_S3.Coefficients.Estimate(2);
S3_R3_acc_p = orig_mdl_S3.Coefficients.pValue(2);

num_shuf = 1000;
beta_hold_S1 = zeros(2,num_shuf); % dim1: R1/R3 acc
beta_hold_S2 = zeros(2,num_shuf); % dim2: R1/R2 acc
beta_hold_S3 = zeros(1,num_shuf); % dim3: R3 acc

num_chans = length(SPK_S1_D_Comb_ansc);

for iShuf = 1:num_shuf
    rand_idx = randperm(num_chans);

    rand_mdl_S1 = fitlm([SFC_T1c_CombZ,SFC_T3c_CombZ],SPK_S1_D_Comb_ansc(rand_idx));
    rand_mdl_S2 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S2_D_Comb_ansc(rand_idx));
    rand_mdl_S3 = fitlm(SFC_T3c_CombZ,SPK_S3_D_Comb_ansc(rand_idx));
    
    beta_hold_S1(:,iShuf) = rand_mdl_S1.Coefficients.Estimate(2:3);
    beta_hold_S2(:,iShuf) = rand_mdl_S2.Coefficients.Estimate(2:3);
    beta_hold_S3(:,iShuf) = rand_mdl_S3.Coefficients.Estimate(2);
end

save(sprintf('%s\\SFC_x_FW_stats\\Ocean_SPK_Mem_FW_x_SFCacc_perm_stats_3Seq_StepLM.mat',loaddir),...
    'beta_hold_*','S1_*','S2_*','S3_*','orig_*')

%% Frontal SFC x SPK Mem FW shuffles, 3-Seq, Groot, 3x predictors stepLM, ANOVA-AllSig
loaddir = 'G:\PaperPrep\PostNatNeuroRevision';

load(sprintf('%s\\Groot3Seq_FW_concat_NewFW_srm_ZC_Frontal_new3_11days_LongD.mat',loaddir));
load(sprintf('%s\\Groot3Seq_SFCvec_concat_BaseComp_Frontal_NoMask_11days_debug.mat',loaddir));
load(sprintf('%s\\Groot3Seq_Encoding_SPK_ANOVA_vec_3Seq.mat',loaddir));
ANOVA_sel = find(~isnan(ANOVA_vec_R3));
SFC_T1c_CombZ = SFC_T1c_CombZ(ANOVA_sel);
SFC_T2c_CombZ = SFC_T2c_CombZ(ANOVA_sel);
SFC_T3c_CombZ = SFC_T3c_CombZ(ANOVA_sel);
SPK_S1_D_Comb_ansc = SPK_S1_D_Comb_ansc(ANOVA_sel);
SPK_S2_D_Comb_ansc = SPK_S2_D_Comb_ansc(ANOVA_sel);
SPK_S3_D_Comb_ansc = SPK_S3_D_Comb_ansc(ANOVA_sel);

tbl = table(SFC_T1c_CombZ,SFC_T2c_CombZ,SFC_T3c_CombZ,...
    SPK_S1_D_Comb_ansc,SPK_S2_D_Comb_ansc,SPK_S3_D_Comb_ansc,...
    'VariableNames',{'SFC_T1acc','SFC_T2acc','SFC_T3acc',...
    'SPK_S1_D','SPK_S2_D','SPK_S3_D'});

% S1 model: SPK_S1_D ~ SFC_T1acc + SFC_T3acc
orig_mdl_S1 = stepwiselm(tbl,'ResponseVar','SPK_S1_D',...
    'PredictorVars',{'SFC_T1acc','SFC_T2acc','SFC_T3acc'},'Upper','linear');
S1_R1_acc_beta = orig_mdl_S1.Coefficients.Estimate(2);
S1_R1_acc_p = orig_mdl_S1.Coefficients.pValue(2);
S1_R3_acc_beta = orig_mdl_S1.Coefficients.Estimate(3);
S1_R3_acc_p = orig_mdl_S1.Coefficients.pValue(3);

% S2 model: SPK_S2_D ~ SFC_T1acc + SFC_T2acc
orig_mdl_S2 = stepwiselm(tbl,'ResponseVar','SPK_S2_D',...
    'PredictorVars',{'SFC_T1acc','SFC_T2acc','SFC_T3acc'},'Upper','linear');
S2_R1_acc_beta = orig_mdl_S2.Coefficients.Estimate(2);
S2_R1_acc_p = orig_mdl_S2.Coefficients.pValue(2);
S2_R2_acc_beta = orig_mdl_S2.Coefficients.Estimate(3);
S2_R2_acc_p = orig_mdl_S2.Coefficients.pValue(3);

% S3 model: SPK_S3_D ~ SFC_T3acc
orig_mdl_S3 = stepwiselm(tbl,'ResponseVar','SPK_S3_D',...
    'PredictorVars',{'SFC_T1acc','SFC_T2acc','SFC_T3acc'},'Upper','linear');
S3_R3_acc_beta = orig_mdl_S3.Coefficients.Estimate(2);
S3_R3_acc_p = orig_mdl_S3.Coefficients.pValue(2);

num_shuf = 1000;
beta_hold_S1 = zeros(2,num_shuf); % dim1: R1/R3 acc
beta_hold_S2 = zeros(2,num_shuf); % dim2: R1/R2 acc
beta_hold_S3 = zeros(1,num_shuf); % dim3: R3 acc

num_chans = length(SPK_S1_D_Comb_ansc);

for iShuf = 1:num_shuf
    rand_idx = randperm(num_chans);

    rand_mdl_S1 = fitlm([SFC_T1c_CombZ,SFC_T3c_CombZ],SPK_S1_D_Comb_ansc(rand_idx));
    rand_mdl_S2 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S2_D_Comb_ansc(rand_idx));
    rand_mdl_S3 = fitlm(SFC_T3c_CombZ,SPK_S3_D_Comb_ansc(rand_idx));
    
    beta_hold_S1(:,iShuf) = rand_mdl_S1.Coefficients.Estimate(2:3);
    beta_hold_S2(:,iShuf) = rand_mdl_S2.Coefficients.Estimate(2:3);
    beta_hold_S3(:,iShuf) = rand_mdl_S3.Coefficients.Estimate(2);
end

save(sprintf('%s\\SFC_x_FW_stats\\Groot_SPK_Mem_FW_x_SFCacc_perm_stats_3Seq_StepLM.mat',loaddir),...
    'beta_hold_*','S1_*','S2_*','S3_*','orig_*')

%% SFC x SPK Mem FW shuffles, 3-Seq Err, Ocean, 3x predictors stepLM, ANOVA-AllSig
loaddir = 'G:\PaperPrep\CurrBiolRevision';

load(sprintf('%s\\Ocean3Seq_FW_concat_NewFW2024_srm_ZC_Frontal_321_adj.mat',loaddir));
load(sprintf('%s\\Ocean3Seq_SFCvec_concat_Frontal_AllErr.mat',loaddir));
load(sprintf('%s\\Ocean3Seq_SPK_ANOVA_concat_Frontal.mat',loaddir));
ANOVA_sel = find(~isnan(ANOVA_vec_R3));
SFC_T1c_CombZ = SFC_T1c_CombZ(ANOVA_sel);
SFC_T2c_CombZ = SFC_T2c_CombZ(ANOVA_sel);
SFC_T3c_CombZ = SFC_T3c_CombZ(ANOVA_sel);
SPK_S1_D_Comb_ansc = SPK_S1_D_Comb_ansc(ANOVA_sel);
SPK_S2_D_Comb_ansc = SPK_S2_D_Comb_ansc(ANOVA_sel);
SPK_S3_D_Comb_ansc = SPK_S3_D_Comb_ansc(ANOVA_sel);

tbl = table(SFC_T1c_CombZ,SFC_T2c_CombZ,SFC_T3c_CombZ,...
    SPK_S1_D_Comb_ansc,SPK_S2_D_Comb_ansc,SPK_S3_D_Comb_ansc,...
    'VariableNames',{'SFC_T1acc','SFC_T2acc','SFC_T3acc',...
    'SPK_S1_D','SPK_S2_D','SPK_S3_D'});

% S1 model: SPK_S1_D ~ SFC_T1acc + SFC_T3acc
orig_mdl_S1 = stepwiselm(tbl,'ResponseVar','SPK_S1_D',...
    'PredictorVars',{'SFC_T1acc','SFC_T2acc','SFC_T3acc'},'Upper','linear');
S1_R1_acc_beta = orig_mdl_S1.Coefficients.Estimate(2);
S1_R1_acc_p = orig_mdl_S1.Coefficients.pValue(2);
S1_R3_acc_beta = orig_mdl_S1.Coefficients.Estimate(3);
S1_R3_acc_p = orig_mdl_S1.Coefficients.pValue(3);

% S2 model: SPK_S2_D ~ SFC_T1acc + SFC_T3acc
orig_mdl_S2 = stepwiselm(tbl,'ResponseVar','SPK_S2_D',...
    'PredictorVars',{'SFC_T1acc','SFC_T2acc','SFC_T3acc'},'Upper','linear');
S2_R1_acc_beta = orig_mdl_S2.Coefficients.Estimate(2);
S2_R1_acc_p = orig_mdl_S2.Coefficients.pValue(2);
S2_R3_acc_beta = orig_mdl_S2.Coefficients.Estimate(3);
S2_R3_acc_p = orig_mdl_S2.Coefficients.pValue(3);

% S3 model: SPK_S3_D ~ SFC_T1acc + SFC_T2acc
orig_mdl_S3 = stepwiselm(tbl,'ResponseVar','SPK_S3_D',...
    'PredictorVars',{'SFC_T1acc','SFC_T2acc','SFC_T3acc'},'Upper','linear');
S3_R1_acc_beta = orig_mdl_S3.Coefficients.Estimate(2);
S3_R1_acc_p = orig_mdl_S3.Coefficients.pValue(2);
S3_R2_acc_beta = orig_mdl_S3.Coefficients.Estimate(2);
S3_R2_acc_p = orig_mdl_S3.Coefficients.pValue(2);

num_shuf = 1000;
beta_hold_S1 = zeros(2,num_shuf); % dim1: R1/R3 acc
beta_hold_S2 = zeros(2,num_shuf); % dim2: R1/R3 acc
beta_hold_S3 = zeros(2,num_shuf); % dim3: R1/R2 acc

num_chans = length(SPK_S1_D_Comb_ansc);

for iShuf = 1:num_shuf
    rand_idx = randperm(num_chans);

    rand_mdl_S1 = fitlm([SFC_T1c_CombZ,SFC_T3c_CombZ],SPK_S1_D_Comb_ansc(rand_idx));
    rand_mdl_S2 = fitlm([SFC_T1c_CombZ,SFC_T3c_CombZ],SPK_S2_D_Comb_ansc(rand_idx));
    rand_mdl_S3 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S3_D_Comb_ansc(rand_idx));
    
    beta_hold_S1(:,iShuf) = rand_mdl_S1.Coefficients.Estimate(2:3);
    beta_hold_S2(:,iShuf) = rand_mdl_S2.Coefficients.Estimate(2:3);
    beta_hold_S3(:,iShuf) = rand_mdl_S3.Coefficients.Estimate(2:3);
end

save(sprintf('%s\\SFC_x_FW_stats\\Ocean_SPK_Mem_FW_x_SFCacc_perm_stats_3Seq_StepLM_Err.mat',loaddir),...
    'beta_hold_*','S1_*','S2_*','S3_*','orig_*')

%% SFC x SPK Mem FW shuffles, 3-Seq Err, Groot, 3x predictors stepLM, ANOVA-AllSig
loaddir = 'G:\PaperPrep\CurrBiolRevision';

load(sprintf('%s\\Groot3Seq_FW_concat_NewFW_srm_ZC_Frontal_new3_11days_LongD.mat',loaddir));
load(sprintf('%s\\Groot3Seq_SFCvec_concat_Frontal_AllErr.mat',loaddir));
load(sprintf('%s\\Groot3Seq_Encoding_SPK_ANOVA_vec_3Seq.mat',loaddir));
ANOVA_sel = find(~isnan(ANOVA_vec_R3));
SFC_T1c_CombZ = SFC_T1c_CombZ(ANOVA_sel);
SFC_T2c_CombZ = SFC_T2c_CombZ(ANOVA_sel);
SFC_T3c_CombZ = SFC_T3c_CombZ(ANOVA_sel);
SPK_S1_D_Comb_ansc = SPK_S1_D_Comb_ansc(ANOVA_sel);
SPK_S2_D_Comb_ansc = SPK_S2_D_Comb_ansc(ANOVA_sel);
SPK_S3_D_Comb_ansc = SPK_S3_D_Comb_ansc(ANOVA_sel);

tbl = table(SFC_T1c_CombZ,SFC_T2c_CombZ,SFC_T3c_CombZ,...
    SPK_S1_D_Comb_ansc,SPK_S2_D_Comb_ansc,SPK_S3_D_Comb_ansc,...
    'VariableNames',{'SFC_T1acc','SFC_T2acc','SFC_T3acc',...
    'SPK_S1_D','SPK_S2_D','SPK_S3_D'});

% S1 model: SPK_S1_D ~ SFC_T1acc
orig_mdl_S1 = stepwiselm(tbl,'ResponseVar','SPK_S1_D',...
    'PredictorVars',{'SFC_T1acc','SFC_T2acc','SFC_T3acc'},'Upper','linear');
S1_R1_acc_beta = orig_mdl_S1.Coefficients.Estimate(2);
S1_R1_acc_p = orig_mdl_S1.Coefficients.pValue(2);

% S2 model: SPK_S2_D ~ SFC_T1acc + SFC_T2acc
orig_mdl_S2 = stepwiselm(tbl,'ResponseVar','SPK_S2_D',...
    'PredictorVars',{'SFC_T1acc','SFC_T2acc'},'Upper','linear');
S2_R1_acc_beta = orig_mdl_S2.Coefficients.Estimate(2);
S2_R1_acc_p = orig_mdl_S2.Coefficients.pValue(2);
S2_R2_acc_beta = orig_mdl_S2.Coefficients.Estimate(3);
S2_R2_acc_p = orig_mdl_S2.Coefficients.pValue(3);

% S3 model: SPK_S3_D ~ SFC_T3acc
orig_mdl_S3 = stepwiselm(tbl,'ResponseVar','SPK_S3_D',...
    'PredictorVars',{'SFC_T1acc','SFC_T2acc','SFC_T3acc'},'Upper','linear');
S3_R3_acc_beta = orig_mdl_S3.Coefficients.Estimate(2);
S3_R3_acc_p = orig_mdl_S3.Coefficients.pValue(2);

num_shuf = 1000;
beta_hold_S1 = zeros(1,num_shuf); % dim1: R1/R3 acc
beta_hold_S2 = zeros(2,num_shuf); % dim2: R1/R2 acc
beta_hold_S3 = zeros(1,num_shuf); % dim3: R3 acc

num_chans = length(SPK_S1_D_Comb_ansc);

for iShuf = 1:num_shuf
    rand_idx = randperm(num_chans);

    rand_mdl_S1 = fitlm(SFC_T1c_CombZ,SPK_S1_D_Comb_ansc(rand_idx));
    rand_mdl_S2 = fitlm([SFC_T1c_CombZ,SFC_T2c_CombZ],SPK_S2_D_Comb_ansc(rand_idx));
    rand_mdl_S3 = fitlm(SFC_T3c_CombZ,SPK_S3_D_Comb_ansc(rand_idx));
    
    beta_hold_S1(:,iShuf) = rand_mdl_S1.Coefficients.Estimate(2);
    beta_hold_S2(:,iShuf) = rand_mdl_S2.Coefficients.Estimate(2:3);
    beta_hold_S3(:,iShuf) = rand_mdl_S3.Coefficients.Estimate(2);
end

save(sprintf('%s\\SFC_x_FW_stats\\Groot_SPK_Mem_FW_x_SFCacc_perm_stats_3Seq_StepLM_Err.mat',loaddir),...
    'beta_hold_*','S1_*','S2_*','S3_*','orig_*')