%% SFC, cross chan, Ocean
addpath(genpath('D:\Scripts\chronux_2_12.v03\chronux_2_12\chronux_2_12'));

file_names = {'ocean2021-08-27','ocean2021-08-30','ocean2021-08-31','ocean2021-09-01','ocean2021-09-02',...
    'ocean2021-09-06','ocean2021-09-07','ocean2021-09-08','ocean2021-09-09','ocean2021-09-10'};
file_tags = {'Ocean0827','Ocean0830','Ocean0831','Ocean0901','Ocean0902',...
    'Ocean0906','Ocean0907','Ocean0908','Ocean0909','Ocean0910'};
spike_dir = 'D:\ocean_data\sorting';
loaddir = 'D:\ocean_data\preprocessed';
savedir = 'D:\ocean_data\preprocessed\prelim_analyses\SFC_CrossChan';
if ~exist(savedir,'dir')
    mkdir(savedir);
end

load(sprintf('%s\\new_area_codes_July2021.mat',loaddir),'good_chans')

% estim_start = 250;
% estim_end = 1300;
% num_bins = ((estim_start+estim_end)/step_size)+1;
% bin_centers = 0:50:(estim_start+estim_end);

% eval_labels = {'Cue_Trig','T1_Trig','T2_Trig'};
% img_bounds_set = {[400,1000,1150],[250,1000,1300],[750,1000,800]};
% event_col = [1,2,4]; % Cue/T1/T2 trig for Encoding

% eval_labels = {'T1_Trig','T2_Trig'};
% img_bounds_set = {[250,1000,1300],[750,1000,800]};
% event_col = [2,4]; % T1/T2 trig for Encoding
% bin_centers = 0:50:1550;


eval_labels = {'Base','T1','T2','T3','GC'};
img_size = 500;
event_col = [1,2,4,6,10]; % pre-CueFrame for Baseline
img_shift_set = [-500,250,250,250,250];
num_ts = length(eval_labels);

eval_size = 250; %+-X ms
    fscorr = 1; % correction for finite number of spikes per trial
    sfreq = 1000;
    trial_LB = 1500;
    num_chans = length(good_chans); 
    
freq_str = 'Theta';
% freq_band = [4,8];

params.Fs = sfreq;
params.tapers = [3,5];
params.trialave = 1;
params.fpass = [2 10];
% num_freq = 4; % changes with tapers
f_select = 2; % 4-8 Hz, assume 4Hz bandwidth;

num_eval = length(eval_labels);
num_files = length(file_tags);
for iFile = 1:num_files
% for iFile = 3:num_files

    curr_file = file_names{iFile};
    curr_tag = file_tags{iFile};
    disp(['Processing ',curr_tag,'...']);
    
    load(sprintf('%s\\trial_cuts\\%s_CueFrameTrialData_LB1500_RB5500_denoise.mat',loaddir,curr_file),'LFP_trial_data','trial_events');
    load(sprintf('%s\\trial_cuts\\%s_TrialInfo_split.mat',loaddir,curr_file),'TrueCorrect','Target');
    load(sprintf('%s\\%s\\Neuron.mat',spike_dir,curr_file),'Spk','Spk_channel');
    
    Correct_idx = find(TrueCorrect==1 & ~isnan(Target(:,3)));
    trial_events = trial_events(Correct_idx,:);
    num_trials = length(Correct_idx);
    LFP_trial_data = LFP_trial_data(good_chans,:,Correct_idx);
    
% for Chronux (coherencycpt.m)
%       data1        (continuous data in time x trials form) -- required
%       data2        (structure array of spike times with dimension trials; 
%                     also accepts 1d array of spike times) -- required
        if params.trialave == 0
            if exist(sprintf('%s\\%s_%s_TimeSeg_SFC_trialFull_temp.mat',savedir,curr_tag,freq_str),'file')
                load(sprintf('%s\\%s_%s_TimeSeg_SFC_trialFull_temp.mat',savedir,curr_tag,freq_str),'C_set','iEval');
                start_bin = iEval+1;
            else
                start_bin = 1;
                C_set = nan(num_chans,num_chans,num_trials,num_ts);
            end
        else
            if exist(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve_temp.mat',savedir,curr_tag,freq_str),'file')
                load(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve_temp.mat',savedir,curr_tag,freq_str),'C_set','iEval')
                start_bin = iEval+1;
            else
                start_bin = 1;
                C_set = nan(num_chans,num_chans,num_ts);
            end
        end

        for iEval = start_bin:num_ts
            tic;
            for iCL = 1:num_chans
                event_var = round((trial_events(:,event_col(iEval)) - trial_events(:,1)).*sfreq)+trial_LB+img_shift_set(iEval); 
                event_ts = trial_events(:,event_col(iEval)) + img_shift_set(iEval)/sfreq;
            
                LFP_recut = zeros(eval_size*2+1,num_trials);
                for i = 1:num_trials
                    LFP_recut(:,i) = squeeze(LFP_trial_data(iCL,event_var(i)-eval_size:event_var(i)+eval_size,i));
                end
                for iChan = 1:num_chans
                    curr_chan = good_chans(iChan); 
                    spike_loc = find(Spk_channel == curr_chan);
                    if isempty(spike_loc) || length(spike_loc) > 1
                        continue
                    end
                    spike_idx = Spk{spike_loc};
                    N = num_trials;
        
                    b = repmat(struct('x',1),N,1);

                    for i = 1:num_trials
                        spike_recut_temp = spike_idx(spike_idx > event_ts(i)-(eval_size/sfreq) & spike_idx < event_ts(i)+(eval_size/sfreq));
                        spike_recut_temp = spike_recut_temp - (event_ts(i)-(eval_size/sfreq)); % realign timestamps
                        b(i).x = spike_recut_temp;
                    end
                    [C,~,~,~,~,f] = coherencycpt(LFP_recut,b,params,fscorr);
                    
                    if params.trialave == 0
                        C_set(iCL,iChan,:,iEval) = C(f_select,:); % 4-8 Hz, assume 4Hz bandwidth;
                    else
                        C_set(iCL,iChan,iEval) = C(f_select);
                    end
                end
            end
            if params.trialave == 0
                save(sprintf('%s\\%s_%s_TimeSeg_SFC_trialFull_temp.mat',savedir,curr_tag,freq_str),'C_set','iEval','-v7.3')
            else
                save(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve_temp.mat',savedir,curr_tag,freq_str),'C_set','iEval','-v7.3')
            end
            toc;
        end
        
        disp('Saving...')   
        if params.trialave == 0
            save(sprintf('%s\\%s_%s_TimeSeg_SFC_trialFull.mat',savedir,curr_tag,freq_str),'C_set','-v7.3')
        else
            save(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve.mat',savedir,curr_tag,freq_str),'C_set','-v7.3')
        end
        
        disp('Done.')
end

%% SFC, cross chan, error trials, Ocean
addpath(genpath('D:\Scripts\chronux_2_12.v03\chronux_2_12\chronux_2_12'));

file_names = {'ocean2021-08-27','ocean2021-08-30','ocean2021-08-31','ocean2021-09-01','ocean2021-09-02',...
    'ocean2021-09-06','ocean2021-09-07','ocean2021-09-08','ocean2021-09-09','ocean2021-09-10'};
file_tags = {'Ocean0827','Ocean0830','Ocean0831','Ocean0901','Ocean0902',...
    'Ocean0906','Ocean0907','Ocean0908','Ocean0909','Ocean0910'};

spike_dir = 'D:\ocean_data\sorting';
loaddir = 'D:\ocean_data\preprocessed';
savedir = 'D:\ocean_data\preprocessed\prelim_analyses\SFC_CrossChan\Error';
if ~exist(savedir,'dir')
    mkdir(savedir);
end

load(sprintf('%s\\new_area_codes_July2021.mat',loaddir),'good_chans')

% estim_start = 250;
% estim_end = 1300;
% num_bins = ((estim_start+estim_end)/step_size)+1;
% bin_centers = 0:50:(estim_start+estim_end);

% eval_labels = {'Cue_Trig','T1_Trig','T2_Trig'};
% img_bounds_set = {[400,1000,1150],[250,1000,1300],[750,1000,800]};
% event_col = [1,2,4]; % Cue/T1/T2 trig for Encoding

% eval_labels = {'T1_Trig','T2_Trig'};
% img_bounds_set = {[250,1000,1300],[750,1000,800]};
% event_col = [2,4]; % T1/T2 trig for Encoding
% bin_centers = 0:50:1550;


eval_labels = {'Base','T1','T2','T3','GC'};
img_size = 500;
event_col = [1,2,4,6,10]; % pre-CueFrame for Baseline
img_shift_set = [-500,250,250,250,250];
num_ts = length(eval_labels);

eval_size = 250; %+-X ms
    fscorr = 1; % correction for finite number of spikes per trial
    sfreq = 1000;
    trial_LB = 1500;
    num_chans = length(good_chans); 
    
freq_str = 'Theta';
% freq_band = [4,8];

params.Fs = sfreq;
params.tapers = [3,5];
params.trialave = 1;
params.fpass = [2 10];
% num_freq = 4; % changes with tapers
f_select = 2; % 4-8 Hz, assume 4Hz bandwidth;

error_labels = {'AllErr','OrderErr'};

num_eval = length(eval_labels);
num_files = length(file_tags);
for iFile = 1:num_files

    curr_file = file_names{iFile};
    curr_tag = file_tags{iFile};
    disp(['Processing ',curr_tag,'...']);
    
    load(sprintf('%s\\trial_cuts\\%s_CueFrameTrialData_LB1500_RB5500_denoise.mat',loaddir,curr_file),'LFP_trial_data','trial_events');
    load(sprintf('%s\\trial_cuts\\%s_TrialInfo_split.mat',loaddir,curr_file),'TrueCorrect','Rule','Response','Target');
    load(sprintf('%s\\%s\\Neuron.mat',spike_dir,curr_file),'Spk','Spk_channel');
    
    Error_idx_all = find(TrueCorrect==0 & ~isnan(Rule) & ~isnan(Response(:,1)) & ~isnan(Target(:,3)));
    Error_idx_order = find(TrueCorrect==0 & ~isnan(Rule) & ~isnan(Response(:,1)) & ~isnan(Target(:,3)) & Target(:,3) == Response(:,2));
    trial_events_orig = trial_events;
    LFP_trial_data_orig = LFP_trial_data;
    
    Error_idx_set = cell(1,2);
    Error_idx_set{1} = Error_idx_all;
    Error_idx_set{2} = Error_idx_order;
    
% for Chronux (coherencycpt.m)
%       data1        (continuous data in time x trials form) -- required
%       data2        (structure array of spike times with dimension trials; 
%                     also accepts 1d array of spike times) -- required
    for iTT = 1:2
        Error_idx = Error_idx_set{iTT};
        err_str = error_labels{iTT};
        trial_events = trial_events_orig(Error_idx,:);
        num_trials = length(Error_idx);
        LFP_trial_data = LFP_trial_data_orig(good_chans,:,Error_idx);
        
        if params.trialave == 0
            if exist(sprintf('%s\\%s_%s_TimeSeg_SFC_trialFull_Err_%s_temp.mat',savedir,curr_tag,freq_str,err_str),'file')
                load(sprintf('%s\\%s_%s_TimeSeg_SFC_trialFull_Err_%s_temp.mat',savedir,curr_tag,freq_str,err_str),'C_set','iEval');
                start_bin = iEval+1;
            else
                start_bin = 1;
                C_set = nan(num_chans,num_chans,num_trials,num_ts);
            end
        else
            if exist(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve_Err_%s_temp.mat',savedir,curr_tag,freq_str,err_str),'file')
                load(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve_Err_%s_temp.mat',savedir,curr_tag,freq_str,err_str),'C_set','iEval')
                start_bin = iEval+1;
            else
                start_bin = 1;
                C_set = nan(num_chans,num_chans,num_ts);
            end
        end

        for iEval = start_bin:num_ts
            tic;
            for iCL = 1:num_chans
                event_var = round((trial_events(:,event_col(iEval)) - trial_events(:,1)).*sfreq)+trial_LB+img_shift_set(iEval); 
                event_ts = trial_events(:,event_col(iEval)) + img_shift_set(iEval)/sfreq;
            
                LFP_recut = zeros(eval_size*2+1,num_trials);
                for i = 1:num_trials
                    LFP_recut(:,i) = squeeze(LFP_trial_data(iCL,event_var(i)-eval_size:event_var(i)+eval_size,i));
                end
                for iChan = 1:num_chans
                    curr_chan = good_chans(iChan); 
                    spike_loc = find(Spk_channel == curr_chan);
                    if isempty(spike_loc) || length(spike_loc) > 1
                        continue
                    end
                    spike_idx = Spk{spike_loc};
                    N = num_trials;
        
                    b = repmat(struct('x',1),N,1);

                    for i = 1:num_trials
                        spike_recut_temp = spike_idx(spike_idx > event_ts(i)-(eval_size/sfreq) & spike_idx < event_ts(i)+(eval_size/sfreq));
                        spike_recut_temp = spike_recut_temp - (event_ts(i)-(eval_size/sfreq)); % realign timestamps
                        b(i).x = spike_recut_temp;
                    end
                    [C,~,~,~,~,f] = coherencycpt(LFP_recut,b,params,fscorr);
                    
                    if params.trialave == 0
                        C_set(iCL,iChan,:,iEval) = C(f_select,:); % 4-8 Hz, assume 4Hz bandwidth;
                    else
                        C_set(iCL,iChan,iEval) = C(f_select);
                    end
                end
            end
            if params.trialave == 0
                save(sprintf('%s\\%s_%s_TimeSeg_SFC_trialFull_Err_%s_temp.mat',savedir,curr_tag,freq_str,err_str),'C_set','iEval','-v7.3')
            else
                save(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve_Err_%s_temp.mat',savedir,curr_tag,freq_str,err_str),'C_set','iEval','-v7.3')
            end
            toc;
        end
        
        disp('Saving...')   
        if params.trialave == 0
            save(sprintf('%s\\%s_%s_TimeSeg_SFC_trialFull_Err_%s.mat',savedir,curr_tag,freq_str,err_str),'C_set','-v7.3')
        else
            save(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve_Err_%s.mat',savedir,curr_tag,freq_str,err_str),'C_set','-v7.3')
        end
        
        disp('Done.')
    end
end

%% SFC, cross chan, Groot
addpath(genpath('D:\Scripts\chronux_2_12.v03\chronux_2_12\chronux_2_12'));

file_names = {'groot2024-03-11','groot2024-03-14','groot2024-03-15','groot2024-03-16','groot2024-03-18',...
    'groot2024-03-19','groot2024-03-20','groot2024-03-21','groot2024-03-22','groot2024-03-23',...
    'groot2024-03-25','groot2024-03-26','groot2024-03-27'};
file_tags = {'Groot0311','Groot0314','Groot0315','Groot0316','Groot0318',...
    'Groot0319','Groot0320','Groot0321','Groot0322','Groot0323',...
    'Groot0325','Groot0326','Groot0327'};
spike_dir = 'E:\Groot_3Seq\sorting';
loaddir = 'E:\Groot_3Seq';
savedir = 'E:\Groot_3Seq\preprocessed\prelim_analyses\SFC_CrossChan';
if ~exist(savedir,'dir')
    mkdir(savedir);
end

% estim_start = 250;
% estim_end = 1300;
% num_bins = ((estim_start+estim_end)/step_size)+1;
% bin_centers = 0:50:(estim_start+estim_end);

% eval_labels = {'Cue_Trig','T1_Trig','T2_Trig'};
% img_bounds_set = {[400,1000,1150],[250,1000,1300],[750,1000,800]};
% event_col = [1,2,4]; % Cue/T1/T2 trig for Encoding

% eval_labels = {'T1_Trig','T2_Trig'};
% img_bounds_set = {[250,1000,1300],[750,1000,800]};
% event_col = [2,4]; % T1/T2 trig for Encoding
% bin_centers = 0:50:1550;


eval_labels = {'Base','T1','T2','T3','GC'};
img_size = 500;
event_col = [1,2,4,6,10]; % pre-CueFrame for Baseline
img_shift_set = [-500,250,250,250,250];
num_ts = length(eval_labels);

eval_size = 250; %+-X ms
    fscorr = 1; % correction for finite number of spikes per trial
    sfreq = 1000;
    trial_LB = 1500;
%     num_chans = length(good_chans); 
    num_chans = 256;
    bad_chans = [15,110,144,159,160];
    
freq_str = 'Theta';
% freq_band = [4,8];

params.Fs = sfreq;
params.tapers = [3,5];
params.trialave = 1;
params.fpass = [2 10];
% num_freq = 4; % changes with tapers
f_select = 2; % 4-8 Hz, assume 4Hz bandwidth;

num_eval = length(eval_labels);
num_files = length(file_tags);
for iFile = 1:num_files

    curr_file = file_names{iFile};
    curr_tag = file_tags{iFile};
    disp(['Processing ',curr_tag,'...']);
    
    load(sprintf('%s\\trial_cuts\\%s_CueFrameTrialData_LB1500_RB5500_denoise.mat',loaddir,curr_file),'LFP_trial_data','trial_events','bad_trial_idx');
    load(sprintf('%s\\trial_cuts\\%s_TrialInfo_split.mat',loaddir,curr_file),'TrueCorrect','Target');
    load(sprintf('%s\\%s\\Neuron.mat',spike_dir,curr_file),'Spk','Spk_channel');
    Spk_channel = Spk_channel(:,1);
    
    bad_trial_vec = logical(sum(bad_trial_idx,1))';
    
    Correct_idx = find(TrueCorrect==1 & ~isnan(Target(:,3)) & bad_trial_vec == 0);
    trial_events = trial_events(Correct_idx,:);
    num_trials = length(Correct_idx);
    LFP_trial_data = LFP_trial_data(:,:,Correct_idx);
    
% for Chronux (coherencycpt.m)
%       data1        (continuous data in time x trials form) -- required
%       data2        (structure array of spike times with dimension trials; 
%                     also accepts 1d array of spike times) -- required
        if params.trialave == 0
            if exist(sprintf('%s\\%s_%s_TimeSeg_SFC_trialFull_temp.mat',savedir,curr_tag,freq_str),'file')
                load(sprintf('%s\\%s_%s_TimeSeg_SFC_trialFull_temp.mat',savedir,curr_tag,freq_str),'C_set','iEval');
                start_bin = iEval+1;
            else
                start_bin = 1;
                C_set = nan(num_chans,num_chans,num_trials,num_ts);
            end
        else
            if exist(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve_temp.mat',savedir,curr_tag,freq_str),'file')
                load(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve_temp.mat',savedir,curr_tag,freq_str),'C_set','iEval')
                start_bin = iEval+1;
            else
                start_bin = 1;
                C_set = nan(num_chans,num_chans,num_ts);
            end
        end

        for iEval = start_bin:num_ts
            tic;
            for iCL = 1:num_chans
                if ismember(iCL,bad_chans)
                    continue
                end
                event_var = round((trial_events(:,event_col(iEval)) - trial_events(:,1)).*sfreq)+trial_LB+img_shift_set(iEval); 
                event_ts = trial_events(:,event_col(iEval)) + img_shift_set(iEval)/sfreq;
            
                LFP_recut = zeros(eval_size*2+1,num_trials);
                for i = 1:num_trials
                    LFP_recut(:,i) = squeeze(LFP_trial_data(iCL,event_var(i)-eval_size:event_var(i)+eval_size,i));
                end
                for iChan = 1:num_chans
                    curr_chan = iChan; 
                    spike_loc = find(Spk_channel == curr_chan);
                    if isempty(spike_loc) || length(spike_loc) > 1
                        continue
                    end
                    spike_idx = Spk{spike_loc};
                    N = num_trials;
        
                    b = repmat(struct('x',1),N,1);

                    for i = 1:num_trials
                        spike_recut_temp = spike_idx(spike_idx > event_ts(i)-(eval_size/sfreq) & spike_idx < event_ts(i)+(eval_size/sfreq));
                        spike_recut_temp = spike_recut_temp - (event_ts(i)-(eval_size/sfreq)); % realign timestamps
                        b(i).x = spike_recut_temp;
                    end
                    [C,~,~,~,~,f] = coherencycpt(LFP_recut,b,params,fscorr);
                    
                    if params.trialave == 0
                        C_set(iCL,iChan,:,iEval) = C(f_select,:); % 4-8 Hz, assume 4Hz bandwidth;
                    else
                        C_set(iCL,iChan,iEval) = C(f_select);
                    end
                end
            end
            if params.trialave == 0
                save(sprintf('%s\\%s_%s_TimeSeg_SFC_trialFull_temp.mat',savedir,curr_tag,freq_str),'C_set','iEval','-v7.3')
            else
                save(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve_temp.mat',savedir,curr_tag,freq_str),'C_set','iEval','-v7.3')
            end
            toc;
        end
        
        disp('Saving...')   
        if params.trialave == 0
            save(sprintf('%s\\%s_%s_TimeSeg_SFC_trialFull.mat',savedir,curr_tag,freq_str),'C_set','-v7.3')
        else
            save(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve.mat',savedir,curr_tag,freq_str),'C_set','-v7.3')
        end
        
        disp('Done.')
end

%% SFC, cross chan, error trials, Groot
addpath(genpath('D:\Scripts\chronux_2_12.v03\chronux_2_12\chronux_2_12'));

file_names = {'groot2024-03-11','groot2024-03-14','groot2024-03-15','groot2024-03-16','groot2024-03-18',...
    'groot2024-03-19','groot2024-03-20','groot2024-03-21','groot2024-03-22','groot2024-03-23',...
    'groot2024-03-25','groot2024-03-26','groot2024-03-27'};
file_tags = {'Groot0311','Groot0314','Groot0315','Groot0316','Groot0318',...
    'Groot0319','Groot0320','Groot0321','Groot0322','Groot0323',...
    'Groot0325','Groot0326','Groot0327'};

spike_dir = 'E:\Groot_3Seq\sorting';
loaddir = 'E:\Groot_3Seq';
savedir = 'E:\Groot_3Seq\preprocessed\prelim_analyses\SFC_CrossChan\Error';
if ~exist(savedir,'dir')
    mkdir(savedir);
end

% estim_start = 250;
% estim_end = 1300;
% num_bins = ((estim_start+estim_end)/step_size)+1;
% bin_centers = 0:50:(estim_start+estim_end);

% eval_labels = {'Cue_Trig','T1_Trig','T2_Trig'};
% img_bounds_set = {[400,1000,1150],[250,1000,1300],[750,1000,800]};
% event_col = [1,2,4]; % Cue/T1/T2 trig for Encoding

% eval_labels = {'T1_Trig','T2_Trig'};
% img_bounds_set = {[250,1000,1300],[750,1000,800]};
% event_col = [2,4]; % T1/T2 trig for Encoding
% bin_centers = 0:50:1550;

eval_labels = {'Base','T1','T2','T3','GC'};
img_size = 500;
event_col = [1,2,4,6,10]; % pre-CueFrame for Baseline
img_shift_set = [-500,250,250,250,250];
num_ts = length(eval_labels);

eval_size = 250; %+-X ms
    fscorr = 1; % correction for finite number of spikes per trial
    sfreq = 1000;
    trial_LB = 1500;
    num_chans = 256; 
    
freq_str = 'Theta';
% freq_band = [4,8];

params.Fs = sfreq;
params.tapers = [3,5];
params.trialave = 1;
params.fpass = [2 10];
% num_freq = 4; % changes with tapers
f_select = 2; % 4-8 Hz, assume 4Hz bandwidth;

error_labels = {'AllErr','OrderErr'};

num_eval = length(eval_labels);
num_files = length(file_tags);
for iFile = 4:num_files

    curr_file = file_names{iFile};
    curr_tag = file_tags{iFile};
    disp(['Processing ',curr_tag,'...']);
    
    load(sprintf('%s\\trial_cuts\\%s_CueFrameTrialData_LB1500_RB5500_denoise.mat',loaddir,curr_file),'LFP_trial_data','trial_events');
    load(sprintf('%s\\trial_cuts\\%s_TrialInfo_split.mat',loaddir,curr_file),'TrueCorrect','Response','Target','CueFrameOnTime');
    load(sprintf('%s\\%s\\Neuron.mat',spike_dir,curr_file),'Spk','Spk_channel');
    Spk_channel = Spk_channel(:,1);
    
    Error_idx_all = find(TrueCorrect==0 & ~isnan(Target(:,3)) & ~isnan(Response(:,1)) & ~isnan(CueFrameOnTime));
    Error_idx_order = find(TrueCorrect==0 & ~isnan(Target(:,3)) & ~isnan(Response(:,1)) & ~isnan(CueFrameOnTime) & (Target(:,3) == Response(:,2) | Target(:,2) == Response(:,1) | Target(:,3) == Response(:,1)));
    trial_events_orig = trial_events;
    LFP_trial_data_orig = LFP_trial_data;
    
    Error_idx_set = cell(1,2);
    Error_idx_set{1} = Error_idx_all;
    Error_idx_set{2} = Error_idx_order;
    
% for Chronux (coherencycpt.m)
%       data1        (continuous data in time x trials form) -- required
%       data2        (structure array of spike times with dimension trials; 
%                     also accepts 1d array of spike times) -- required
    for iTT = 1:2
        Error_idx = Error_idx_set{iTT};
        err_str = error_labels{iTT};
        trial_events = trial_events_orig(Error_idx,:);
        num_trials = length(Error_idx);
        LFP_trial_data = LFP_trial_data_orig(:,:,Error_idx);
        
        if params.trialave == 0
            if exist(sprintf('%s\\%s_%s_TimeSeg_SFC_trialFull_Err_%s_temp.mat',savedir,curr_tag,freq_str,err_str),'file')
                load(sprintf('%s\\%s_%s_TimeSeg_SFC_trialFull_Err_%s_temp.mat',savedir,curr_tag,freq_str,err_str),'C_set','iEval');
                start_bin = iEval+1;
            else
                start_bin = 1;
                C_set = nan(num_chans,num_chans,num_trials,num_ts);
            end
        else
            if exist(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve_Err_%s_temp.mat',savedir,curr_tag,freq_str,err_str),'file')
                load(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve_Err_%s_temp.mat',savedir,curr_tag,freq_str,err_str),'C_set','iEval')
                start_bin = iEval+1;
            else
                start_bin = 1;
                C_set = nan(num_chans,num_chans,num_ts);
            end
        end

        for iEval = start_bin:num_ts
            tic;
            for iCL = 1:num_chans
                event_var = round((trial_events(:,event_col(iEval)) - trial_events(:,1)).*sfreq)+trial_LB+img_shift_set(iEval); 
                event_ts = trial_events(:,event_col(iEval)) + img_shift_set(iEval)/sfreq;
            
                LFP_recut = zeros(eval_size*2+1,num_trials);
                for i = 1:num_trials
                    if isnan(event_var(i))
                        continue
                    end
                    LFP_recut(:,i) = squeeze(LFP_trial_data(iCL,event_var(i)-eval_size:event_var(i)+eval_size,i));
                end
                LFP_recut(:,sum(LFP_recut,1)==0) = [];
                if isempty(LFP_recut)
                    continue
                end
                num_trials = size(LFP_recut,2);
                for iChan = 1:num_chans
                    curr_chan = iChan; 
                    spike_loc = find(Spk_channel == curr_chan);
                    if isempty(spike_loc) || length(spike_loc) > 1
                        continue
                    end
                    spike_idx = Spk{spike_loc};
                    N = num_trials;
        
                    b = repmat(struct('x',1),N,1);

                    for i = 1:num_trials
                        spike_recut_temp = spike_idx(spike_idx > event_ts(i)-(eval_size/sfreq) & spike_idx < event_ts(i)+(eval_size/sfreq));
                        spike_recut_temp = spike_recut_temp - (event_ts(i)-(eval_size/sfreq)); % realign timestamps
                        b(i).x = spike_recut_temp;
                    end
                    [C,~,~,~,~,f] = coherencycpt(LFP_recut,b,params,fscorr);
                    
                    if params.trialave == 0
                        C_set(iCL,iChan,:,iEval) = C(f_select,:); % 4-8 Hz, assume 4Hz bandwidth;
                    else
                        C_set(iCL,iChan,iEval) = C(f_select);
                    end
                end
            end
            if params.trialave == 0
                save(sprintf('%s\\%s_%s_TimeSeg_SFC_trialFull_Err_%s_temp.mat',savedir,curr_tag,freq_str,err_str),'C_set','iEval','-v7.3')
            else
                save(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve_Err_%s_temp.mat',savedir,curr_tag,freq_str,err_str),'C_set','iEval','-v7.3')
            end
            toc;
        end
        
        disp('Saving...')   
        if params.trialave == 0
            save(sprintf('%s\\%s_%s_TimeSeg_SFC_trialFull_Err_%s.mat',savedir,curr_tag,freq_str,err_str),'C_set','-v7.3')
        else
            save(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve_Err_%s.mat',savedir,curr_tag,freq_str,err_str),'C_set','-v7.3')
        end
        
        disp('Done.')
    end
end

%% SFC x SPK Mem FW shuffles, 3-Seq Err, Ocean, get SFC vectors first
SFC_dir = 'G:\FW_data\Seq3\SFC_CrossChan_3Seq\SFC_CrossChan\Error';

file_tags = {'0827','0830','0831','0901','0902','0906','0907','0908','0909','0910'};
savedir = 'G:\PaperPrep\CurrBiolRevision';
if ~exist(savedir,'dir')
    mkdir(savedir)
end
spike_dir = 'G:\FW_data\NewJuly\ocean_data\sorting';

load('C:\Users\DELL\Documents\FW_data\ocean_data\preprocessed\new_area_codes_July2021.mat','good_chans','area_code');

idx_dir = 'C:\Users\DELL\Documents\FW_data\sharing\Unit_PEV\PostSept\tdr_frontal_spk_L3\tdr_sens';

SFC_T1r_Comb = []; SFC_T2r_Comb = []; SFC_T3r_Comb = []; 
SFC_T1c_Comb = []; SFC_T2c_Comb = []; SFC_T3c_Comb = [];
SFC_T1r_CombZ = []; SFC_T2r_CombZ = []; SFC_T3r_CombZ = []; 
SFC_T1c_CombZ = []; SFC_T2c_CombZ = []; SFC_T3c_CombZ = [];
SFC_EntR_Comb = []; SFC_EntC_Comb = []; 
SFC_EntR_CombZ = []; SFC_EntC_CombZ = [];
SFC_T1d_Comb = []; SFC_T2d_Comb = []; SFC_T3d_Comb = []; SFC_EntD_Comb = [];
SFC_T1d_CombZ = []; SFC_T2d_CombZ = []; SFC_T3d_CombZ = []; SFC_EntD_CombZ = [];
num_files = length(file_tags);
for iFile = 1:num_files
    curr_tag = file_tags{iFile};
    
%     load(sprintf('%s\\Ocean%s_Theta_TimeSeg_SFC_trialAve_Err_AllErr.mat',SFC_dir,curr_tag),'C_set')
    load(sprintf('%s\\Ocean%s_Theta_TimeSeg_SFC_trialAve_Err_OrderErr.mat',SFC_dir,curr_tag),'C_set')

    load(sprintf('%s\\ocean2021-%s-%s\\Neuron.mat',spike_dir,curr_tag(end-3:end-2),curr_tag(end-1:end)),'Spk_channel')
    
    load(sprintf('%s\\%s-%s_weight_enc250.mat',idx_dir,curr_tag(end-3:end-2),curr_tag(end-1:end)),'ch_info_all');
    chIn256 = double(cell2mat(ch_info_all(:,2)));
    
    [~, uniqueIdx] = unique(Spk_channel);
    dupeIdx = ismember( Spk_channel, Spk_channel( setdiff( 1:numel(Spk_channel), uniqueIdx ) ) );
    dupes = unique(Spk_channel(dupeIdx)); dupe_locs = find(dupeIdx);
    Spk_chan_trim = setdiff(chIn256,dupes);
    [~,uPFC_sel] = intersect(good_chans,Spk_chan_trim);
    
    mat_B_mean = squeeze(C_set(uPFC_sel,uPFC_sel,1));
    mat_T1_mean = squeeze(C_set(uPFC_sel,uPFC_sel,2));
    mat_T2_mean = squeeze(C_set(uPFC_sel,uPFC_sel,3));
    mat_T3_mean = squeeze(C_set(uPFC_sel,uPFC_sel,4));

    T1_mat_trim = mat_T1_mean;
    T2_mat_trim = mat_T2_mean;
    T3_mat_trim = mat_T3_mean;
    Ent_mat_trim = mat_T1_mean.*mat_T2_mean.*mat_T3_mean;
    T1_diag = diag(T1_mat_trim); T2_diag = diag(T2_mat_trim); T3_diag = diag(T3_mat_trim);
    Ent_diag = diag(Ent_mat_trim);
    
    row_vec_T1 = sum(T1_mat_trim,2); row_vec_T2 = sum(T2_mat_trim,2); row_vec_T3 = sum(T3_mat_trim,2);
    col_vec_T1 = sum(T1_mat_trim,1); col_vec_T2 = sum(T2_mat_trim,1); col_vec_T3 = sum(T3_mat_trim,1);
    ent_vec_r = sum(Ent_mat_trim,2); ent_vec_c = sum(Ent_mat_trim,1);
    
    SFC_T1r_Comb = [SFC_T1r_Comb;row_vec_T1]; SFC_T2r_Comb = [SFC_T2r_Comb;row_vec_T2]; 
    SFC_T3r_Comb = [SFC_T3r_Comb;row_vec_T3]; 
    SFC_T1c_Comb = [SFC_T1c_Comb;col_vec_T1']; SFC_T2c_Comb = [SFC_T2c_Comb;col_vec_T2'];
    SFC_T3c_Comb = [SFC_T3c_Comb;col_vec_T3'];
    SFC_EntR_Comb = [SFC_EntR_Comb;ent_vec_r]; SFC_EntC_Comb = [SFC_EntC_Comb;ent_vec_c'];
    SFC_T1d_Comb = [SFC_T1d_Comb;T1_diag]; SFC_T2d_Comb = [SFC_T2d_Comb;T2_diag]; 
    SFC_T3d_Comb = [SFC_T3d_Comb;T3_diag]; 
    SFC_EntD_Comb = [SFC_EntD_Comb;Ent_diag];

    SFC_T1r_CombZ = [SFC_T1r_CombZ;anscombe(row_vec_T1)]; SFC_T2r_CombZ = [SFC_T2r_CombZ;anscombe(row_vec_T2)]; 
    SFC_T3r_CombZ = [SFC_T3r_CombZ;anscombe(row_vec_T3)]; 
    SFC_T1c_CombZ = [SFC_T1c_CombZ;anscombe(col_vec_T1')]; SFC_T2c_CombZ = [SFC_T2c_CombZ;anscombe(col_vec_T2')];
    SFC_T3c_CombZ = [SFC_T3c_CombZ;anscombe(col_vec_T3')];
    SFC_EntR_CombZ = [SFC_EntR_CombZ;anscombe(ent_vec_r)]; SFC_EntC_CombZ = [SFC_EntC_CombZ;anscombe(ent_vec_c')];
    SFC_T1d_CombZ = [SFC_T1d_CombZ;anscombe(T1_diag)]; SFC_T2d_CombZ = [SFC_T2d_CombZ;anscombe(T2_diag)]; 
    SFC_T3d_CombZ = [SFC_T3d_CombZ;anscombe(T3_diag)]; 
    SFC_EntD_CombZ = [SFC_EntD_CombZ;anscombe(Ent_diag)];
end

% save(sprintf('%s\\Ocean3Seq_SFCvec_concat_Frontal_AllErr.mat',savedir),'*_Comb','*_CombZ');
save(sprintf('%s\\Ocean3Seq_SFCvec_concat_Frontal_OrderErr.mat',savedir),'*_Comb','*_CombZ');

%% SFC x SPK Mem FW shuffles, 3-Seq Err, Groot, get SFC vectors first
SFC_dir = 'E:\ZC_data\groot_data\preprocessed\prelim_analyses\SFC_CrossChan_3Seq_debug\Error';
file_tags = {'0311','0314','0315','0316','0319','0321','0322','0323','0325','0326','0327'};

savedir = 'G:\PaperPrep\CurrBiolRevision';
if ~exist(savedir,'dir')
    mkdir(savedir)
end
spike_dir = 'E:\ZC_data\groot_data\sorting';

load('C:\Users\DELL\Documents\ZC_data\groot_data\new_area_codes_Mar2024.mat','good_chans','area_code');

idx_dir = 'E:\ZC_data\groot_data\sharing\13days_fr_L3_win100step50_frontal';

SFC_T1r_Comb = []; SFC_T2r_Comb = []; SFC_T3r_Comb = []; 
SFC_T1c_Comb = []; SFC_T2c_Comb = []; SFC_T3c_Comb = [];
SFC_T1r_CombZ = []; SFC_T2r_CombZ = []; SFC_T3r_CombZ = []; 
SFC_T1c_CombZ = []; SFC_T2c_CombZ = []; SFC_T3c_CombZ = [];
SFC_EntR_Comb = []; SFC_EntC_Comb = []; 
SFC_EntR_CombZ = []; SFC_EntC_CombZ = [];
SFC_T1d_Comb = []; SFC_T2d_Comb = []; SFC_T3d_Comb = []; SFC_EntD_Comb = [];
SFC_T1d_CombZ = []; SFC_T2d_CombZ = []; SFC_T3d_CombZ = []; SFC_EntD_CombZ = [];
num_files = length(file_tags);
for iFile = 1:num_files
    curr_tag = file_tags{iFile};
    
%     load(sprintf('%s\\Groot%s_Theta_TimeSeg_SFC_trialAve_Err_AllErr.mat',SFC_dir,curr_tag),'C_set')
    load(sprintf('%s\\Groot%s_Theta_TimeSeg_SFC_trialAve_Err_OrderErr.mat',SFC_dir,curr_tag),'C_set')

    C_set = C_set(good_chans,good_chans,:);

    load(sprintf('%s\\groot2024-%s-%s\\Neuron.mat',spike_dir,curr_tag(end-3:end-2),curr_tag(end-1:end)),'Spk_channel')
    Spk_channel = Spk_channel(:,1);
    
    load(sprintf('%s\\%s_ch_idx.mat',idx_dir,curr_tag),'chIn256');
    
    [~, uniqueIdx] = unique(Spk_channel);
    dupeIdx = ismember( Spk_channel, Spk_channel( setdiff( 1:numel(Spk_channel), uniqueIdx ) ) );
    dupes = unique(Spk_channel(dupeIdx)); dupe_locs = find(dupeIdx);
    Spk_chan_trim = setdiff(chIn256,dupes);
    [~,uPFC_sel] = intersect(good_chans,Spk_chan_trim);
    
    mat_B_mean = squeeze(C_set(uPFC_sel,uPFC_sel,1));
    mat_T1_mean = squeeze(C_set(uPFC_sel,uPFC_sel,2));
    mat_T2_mean = squeeze(C_set(uPFC_sel,uPFC_sel,3));
    mat_T3_mean = squeeze(C_set(uPFC_sel,uPFC_sel,4));

    T1_mat_trim = mat_T1_mean;
    T2_mat_trim = mat_T2_mean;
    T3_mat_trim = mat_T3_mean;
    Ent_mat_trim = mat_T1_mean.*mat_T2_mean.*mat_T3_mean;
    T1_diag = diag(T1_mat_trim); T2_diag = diag(T2_mat_trim); T3_diag = diag(T3_mat_trim);
    Ent_diag = diag(Ent_mat_trim);
    
    row_vec_T1 = nansum(T1_mat_trim,2); row_vec_T2 = nansum(T2_mat_trim,2); row_vec_T3 = nansum(T3_mat_trim,2);
    col_vec_T1 = nansum(T1_mat_trim,1); col_vec_T2 = nansum(T2_mat_trim,1); col_vec_T3 = nansum(T3_mat_trim,1);
    ent_vec_r = nansum(Ent_mat_trim,2); ent_vec_c = nansum(Ent_mat_trim,1);
    
    SFC_T1r_Comb = [SFC_T1r_Comb;row_vec_T1]; SFC_T2r_Comb = [SFC_T2r_Comb;row_vec_T2]; 
    SFC_T3r_Comb = [SFC_T3r_Comb;row_vec_T3]; 
    SFC_T1c_Comb = [SFC_T1c_Comb;col_vec_T1']; SFC_T2c_Comb = [SFC_T2c_Comb;col_vec_T2'];
    SFC_T3c_Comb = [SFC_T3c_Comb;col_vec_T3'];
    SFC_EntR_Comb = [SFC_EntR_Comb;ent_vec_r]; SFC_EntC_Comb = [SFC_EntC_Comb;ent_vec_c'];
    SFC_T1d_Comb = [SFC_T1d_Comb;T1_diag]; SFC_T2d_Comb = [SFC_T2d_Comb;T2_diag]; 
    SFC_T3d_Comb = [SFC_T3d_Comb;T3_diag]; 
    SFC_EntD_Comb = [SFC_EntD_Comb;Ent_diag];

    SFC_T1r_CombZ = [SFC_T1r_CombZ;anscombe(row_vec_T1)]; SFC_T2r_CombZ = [SFC_T2r_CombZ;anscombe(row_vec_T2)]; 
    SFC_T3r_CombZ = [SFC_T3r_CombZ;anscombe(row_vec_T3)]; 
    SFC_T1c_CombZ = [SFC_T1c_CombZ;anscombe(col_vec_T1')]; SFC_T2c_CombZ = [SFC_T2c_CombZ;anscombe(col_vec_T2')];
    SFC_T3c_CombZ = [SFC_T3c_CombZ;anscombe(col_vec_T3')];
    SFC_EntR_CombZ = [SFC_EntR_CombZ;anscombe(ent_vec_r)]; SFC_EntC_CombZ = [SFC_EntC_CombZ;anscombe(ent_vec_c')];
    SFC_T1d_CombZ = [SFC_T1d_CombZ;anscombe(T1_diag)]; SFC_T2d_CombZ = [SFC_T2d_CombZ;anscombe(T2_diag)]; 
    SFC_T3d_CombZ = [SFC_T3d_CombZ;anscombe(T3_diag)]; 
    SFC_EntD_CombZ = [SFC_EntD_CombZ;anscombe(Ent_diag)];
end

% save(sprintf('%s\\Groot3Seq_SFCvec_concat_Frontal_AllErr.mat',savedir),'*_Comb','*_CombZ');
save(sprintf('%s\\Groot3Seq_SFCvec_concat_Frontal_OrderErr.mat',savedir),'*_Comb','*_CombZ');
