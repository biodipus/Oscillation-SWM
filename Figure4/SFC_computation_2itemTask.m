%% SFC, cross chan, Ocean 2-item
addpath(genpath('C:\Users\DELL\Documents\scripts\chronux_2_12.v03\chronux_2_12\chronux_2_12'));

% file_names = {'ocean2021-07-09','ocean2021-07-10','ocean2021-07-13','ocean2021-07-14',...
%     'ocean2021-07-15','ocean2021-07-20','ocean2021-07-21','ocean2021-07-22',...
%     'ocean2021-07-23','ocean2021-07-26','ocean2021-07-27','ocean2021-07-29','ocean2021-07-30'};
% file_tags = {'Ocean0709','Ocean0710','Ocean0713','Ocean0714','Ocean0715',...
%     'Ocean0720','Ocean0721','Ocean0722','Ocean0723','Ocean0726',...
%     'Ocean0727','Ocean0729','Ocean0730'};
file_names = {'ocean2021-07-09','ocean2021-07-13','ocean2021-07-14',...
    'ocean2021-07-20','ocean2021-07-21','ocean2021-07-22',...
    'ocean2021-07-26','ocean2021-07-27','ocean2021-07-29','ocean2021-07-30'};
file_tags = {'Ocean0709','Ocean0713','Ocean0714',...
    'Ocean0720','Ocean0721','Ocean0722','Ocean0726',...
    'Ocean0727','Ocean0729','Ocean0730'};
spike_dir = 'G:\FW_data\NewJuly\ocean_data\sorting';
loaddir = 'G:\FW_data\NewJuly\ocean_data\preprocessed';

savedir = 'G:\FW_data\NewJuly\ocean_data\preprocessed\prelim_analyses\SFC_CrossChan'; % Theta
% savedir = 'G:\FW_data\NewJuly\ocean_data\preprocessed\prelim_analyses\SFC_CrossChan_Alpha';
% savedir = 'G:\FW_data\NewJuly\ocean_data\preprocessed\prelim_analyses\SFC_CrossChan_Beta';

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


eval_labels = {'Base','T1','T2','RC','GC'};
img_size = 500;
event_col = [1,2,4,8,10]; % pre-CueFrame for Baseline
img_shift_set = [-500,250,250,250,250];
num_ts = length(eval_labels);

eval_size = 250; %+-X ms
    fscorr = 1; % correction for finite number of spikes per trial
    sfreq = 1000;
    trial_LB = 1500;
    num_chans = length(good_chans); 
    
freq_str = 'Theta';
% freq_str = 'Alpha';
% freq_str = 'Beta';
% freq_band = [4,8];

params.Fs = sfreq;
params.tapers = [3,5];
params.trialave = 1;

% for Theta
params.fpass = [2 10];
% num_freq = 4; % changes with tapers
f_select = 2; % 4-8 Hz, assume 4Hz bandwidth;

% for Alpha
% params.fpass = [4 16];
% f_select = 3; % 8-12 Hz, assume 4Hz bandwidth;

% for (low) Beta
% params.fpass = [8 24];
% f_select = 3:5; % 12-20 Hz, assume 4Hz bandwidth;

num_eval = length(eval_labels);
num_files = length(file_tags);
for iFile = 1:num_files
% for iFile = 3:num_files

    curr_file = file_names{iFile};
    curr_tag = file_tags{iFile};
    disp(['Processing ',curr_tag,'...']);
    
    load(sprintf('%s\\trial_cuts\\%s_CueFrameTrialData_LB1500_RB5500_denoise.mat',loaddir,curr_file),'LFP_trial_data','trial_events');
    load(sprintf('%s\\trial_cuts\\%s_TrialInfo_split.mat',loaddir,curr_file),'TrueCorrect','Rule');
    load(sprintf('%s\\%s\\Neuron.mat',spike_dir,curr_file),'Spk','Spk_channel');
    
    Correct_idx = find(TrueCorrect==1 & ~isnan(Rule));
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
                        if length(f_select) == 1
                            C_set(iCL,iChan,:,iEval) = C(f_select,:); % 4-8 Hz, assume 4Hz bandwidth;
                        else
                            C_set(iCL,iChan,:,iEval) = mean(C(f_select,:),1);
                        end
                    else
                        if length(f_select) == 1
                            C_set(iCL,iChan,iEval) = C(f_select);
                        else
                            C_set(iCL,iChan,iEval) = mean(C(f_select));
                        end
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

%% SFC, cross chan, error trials, Ocean 2-item
addpath(genpath('C:\Users\DELL\Documents\scripts\chronux_2_12.v03\chronux_2_12\chronux_2_12'));

file_names = {'ocean2021-07-09','ocean2021-07-10','ocean2021-07-13','ocean2021-07-14',...
    'ocean2021-07-15','ocean2021-07-20','ocean2021-07-21','ocean2021-07-22',...
    'ocean2021-07-23','ocean2021-07-26','ocean2021-07-27','ocean2021-07-29','ocean2021-07-30'};
file_tags = {'Ocean0709','Ocean0710','Ocean0713','Ocean0714','Ocean0715',...
    'Ocean0720','Ocean0721','Ocean0722','Ocean0723','Ocean0726',...
    'Ocean0727','Ocean0729','Ocean0730'};
spike_dir = 'G:\FW_data\NewJuly\ocean_data\sorting';
loaddir = 'G:\FW_data\NewJuly\ocean_data\preprocessed';
savedir = 'G:\FW_data\NewJuly\ocean_data\preprocessed\prelim_analyses\SFC_CrossChan\Error';
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


eval_labels = {'Base','T1','T2','RC','GC'};
img_size = 500;
event_col = [1,2,4,8,10]; % pre-CueFrame for Baseline
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

error_labels = {'AllErr','OrderErr','ItemErr'};

num_eval = length(eval_labels);
num_files = length(file_tags);
for iFile = 1:num_files
% for iFile = 3:num_files

    curr_file = file_names{iFile};
    curr_tag = file_tags{iFile};
    disp(['Processing ',curr_tag,'...']);
    
    load(sprintf('%s\\trial_cuts\\%s_CueFrameTrialData_LB1500_RB5500_denoise.mat',loaddir,curr_file),'LFP_trial_data','trial_events');
    load(sprintf('%s\\trial_cuts\\%s_TrialInfo_split.mat',loaddir,curr_file),'TrueCorrect','Rule','Response','Target');
    load(sprintf('%s\\%s\\Neuron.mat',spike_dir,curr_file),'Spk','Spk_channel');
    
    Error_idx_all = find(TrueCorrect==0 & ~isnan(Rule) & ~isnan(Response(:,1)));
    Error_idx_order = find(TrueCorrect==0 & ~isnan(Rule) & ~isnan(Response(:,1)) & (Target(:,1) == Response(:,2) | Target(:,2) == Response(:,1)));
    Error_idx_item = setdiff(Error_idx_all,Error_idx_order);
    trial_events_orig = trial_events;
    LFP_trial_data_orig = LFP_trial_data;
    
    Error_idx_set = cell(1,3);
    Error_idx_set{1} = Error_idx_all;
    Error_idx_set{2} = Error_idx_order;
    Error_idx_set{3} = Error_idx_item;
    
% for Chronux (coherencycpt.m)
%       data1        (continuous data in time x trials form) -- required
%       data2        (structure array of spike times with dimension trials; 
%                     also accepts 1d array of spike times) -- required
    for iTT = 1:3
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

%% SFC, cross chan, Groot 2-item
addpath(genpath('C:\Users\DELL\Documents\scripts\chronux_2_12.v03\chronux_2_12\chronux_2_12'));

file_names = {'groot2024-04-15','groot2024-05-07','groot2024-05-09',...
    'groot2024-05-13','groot2024-05-14','groot2024-05-15',...
    'groot2024-05-16','groot2024-05-20','groot2024-05-21',...
    'groot2024-05-22','groot2024-05-23'};
file_tags = {'Groot0415','Groot0507','Groot0509',...
    'Groot0513','Groot0514','Groot0515',...
    'Groot0516','Groot0520','Groot0521',...
    'Groot0522','Groot0523'};

spike_dir = 'E:\ZC_data\groot_data_switch\Groot_switch_sorting\sorting';
loaddir = 'E:\ZC_data\groot_data_switch';

savedir = 'E:\ZC_data\groot_data_switch\prelim_analyses\SFC_CrossChan';
if ~exist(savedir,'dir')
    mkdir(savedir);
end

load(sprintf('%s\\new_area_codes_Mar2024.mat',loaddir),'good_chans')

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


eval_labels = {'Base','T1','T2'};
img_size = 500;
event_col = [1,2,4]; % pre-CueFrame for Baseline
img_shift_set = [-500,250,250,250,250];
num_ts = length(eval_labels);

eval_size = 250; %+-X ms
    fscorr = 1; % correction for finite number of spikes per trial
    sfreq = 1000;
    trial_LB = 1500;
    num_chans = length(good_chans); 
    
freq_str = 'Theta';
% freq_str = 'Alpha';
% freq_str = 'Beta';
% freq_band = [4,8];

params.Fs = sfreq;
params.tapers = [3,5];
params.trialave = 1;

% for Theta
params.fpass = [2 10];
% num_freq = 4; % changes with tapers
f_select = 2; % 4-8 Hz, assume 4Hz bandwidth;

% for Alpha
% params.fpass = [4 16];
% f_select = 3; % 8-12 Hz, assume 4Hz bandwidth;

% for (low) Beta
% params.fpass = [8 24];
% f_select = 3:5; % 12-20 Hz, assume 4Hz bandwidth;

num_eval = length(eval_labels);
num_files = length(file_tags);
for iFile = 1:num_files
% for iFile = 3:num_files

    curr_file = file_names{iFile};
    curr_tag = file_tags{iFile};
    disp(['Processing ',curr_tag,'...']);
    
    load(sprintf('%s\\trial_cuts\\%s_CueFrameTrialData_LB1500_RB5500_denoise.mat',loaddir,curr_file),'LFP_trial_data','trial_events');
    load(sprintf('%s\\trial_cuts\\%s_TrialInfo_split.mat',loaddir,curr_file),'TrueCorrect','Rule','CueFrameOnTime');
    load(sprintf('%s\\%s\\Neuron.mat',spike_dir,curr_file),'Spk','Spk_channel');
    Spk_channel = Spk_channel(:,1);
    
    Correct_idx = find(TrueCorrect==1 & ~isnan(Rule) & ~isnan(CueFrameOnTime));
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
                        if length(f_select) == 1
                            C_set(iCL,iChan,:,iEval) = C(f_select,:); % 4-8 Hz, assume 4Hz bandwidth;
                        else
                            C_set(iCL,iChan,:,iEval) = mean(C(f_select,:),1);
                        end
                    else
                        if length(f_select) == 1
                            C_set(iCL,iChan,iEval) = C(f_select);
                        else
                            C_set(iCL,iChan,iEval) = mean(C(f_select));
                        end
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

%% SFC, cross chan, error trials, Groot 2-item
addpath(genpath('C:\Users\DELL\Documents\scripts\chronux_2_12.v03\chronux_2_12\chronux_2_12'));

file_names = {'groot2024-04-15','groot2024-05-07','groot2024-05-09',...
    'groot2024-05-13','groot2024-05-14','groot2024-05-15',...
    'groot2024-05-16','groot2024-05-20','groot2024-05-21',...
    'groot2024-05-22','groot2024-05-23'};
file_tags = {'Groot0415','Groot0507','Groot0509',...
    'Groot0513','Groot0514','Groot0515',...
    'Groot0516','Groot0520','Groot0521',...
    'Groot0522','Groot0523'};

spike_dir = 'E:\ZC_data\groot_data_switch\Groot_switch_sorting\sorting';

loaddir = 'E:\ZC_data\groot_data_switch';
savedir = 'E:\ZC_data\groot_data_switch\prelim_analyses\SFC_CrossChan\Error';
if ~exist(savedir,'dir')
    mkdir(savedir);
end

load(sprintf('%s\\new_area_codes_Mar2024.mat',loaddir),'good_chans')

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


eval_labels = {'Base','T1','T2'};
img_size = 500;
event_col = [1,2,4]; % pre-CueFrame for Baseline
img_shift_set = [-500,250,250];
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

error_labels = {'AllErr','OrderErr','ItemErr'};

num_eval = length(eval_labels);
num_files = length(file_tags);
for iFile = 1:num_files
% for iFile = 3:num_files

    curr_file = file_names{iFile};
    curr_tag = file_tags{iFile};
    disp(['Processing ',curr_tag,'...']);
    
    load(sprintf('%s\\trial_cuts\\%s_CueFrameTrialData_LB1500_RB5500_denoise.mat',loaddir,curr_file),'LFP_trial_data','trial_events');
    load(sprintf('%s\\trial_cuts\\%s_TrialInfo_split.mat',loaddir,curr_file),'TrueCorrect','Rule','Response','Target');
    load(sprintf('%s\\%s\\Neuron.mat',spike_dir,curr_file),'Spk','Spk_channel');
    Spk_channel = Spk_channel(:,1);
    
    Error_idx_all = find(TrueCorrect==0 & ~isnan(Rule) & ~isnan(Response(:,1)) & isnan(Target(:,3)));
    Error_idx_order = find(TrueCorrect==0 & ~isnan(Rule) & isnan(Target(:,3)) & ~isnan(Response(:,1)) & (Target(:,1) == Response(:,2) | Target(:,2) == Response(:,1)));
    Error_idx_item = setdiff(Error_idx_all,Error_idx_order);
    trial_events_orig = trial_events;
    LFP_trial_data_orig = LFP_trial_data;
    
    Error_idx_set = cell(1,3);
    Error_idx_set{1} = Error_idx_all;
    Error_idx_set{2} = Error_idx_order;
    Error_idx_set{3} = Error_idx_item;
    
% for Chronux (coherencycpt.m)
%       data1        (continuous data in time x trials form) -- required
%       data2        (structure array of spike times with dimension trials; 
%                     also accepts 1d array of spike times) -- required
    for iTT = 1:3
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

%% SFC cross chan, sliding over time, Ocean 2-item
addpath(genpath('C:\Users\DELL\Documents\scripts\chronux_2_12.v03\chronux_2_12\chronux_2_12'));

file_names = {'ocean2021-07-09','ocean2021-07-13','ocean2021-07-14',...
    'ocean2021-07-20','ocean2021-07-21','ocean2021-07-22',...
    'ocean2021-07-26','ocean2021-07-27','ocean2021-07-29','ocean2021-07-30'};
file_tags = {'Ocean0709','Ocean0713','Ocean0714',...
    'Ocean0720','Ocean0721','Ocean0722','Ocean0726',...
    'Ocean0727','Ocean0729','Ocean0730'};

spike_dir = 'G:\FW_data\NewJuly\ocean_data\sorting';
loaddir = 'C:\Users\DELL\Documents\FW_data\ocean_data\preprocessed';
loaddir2 = 'G:\FW_data\NewJuly\ocean_data\preprocessed\trial_cuts';
savedir = 'C:\Users\DELL\Documents\FW_data\ocean_data\preprocessed\prelim_analyses\SpikeWork\SFC_CrossChan';
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

eval_labels = {'T1_Trig'};
img_bounds_set = {[250,1000,1300]};
event_col = [2]; % T1/T2 trig for Encoding
bin_centers = 0:50:1550;


% eval_labels = {'Base'};
% img_bounds_set = {[750,1000,250]};
% event_col = [1]; % CueFrame for Baseline
% bin_centers = 0:50:500;

num_ts = length(bin_centers);
eval_size = 250; %+-X ms
% eval_starts = [-1000,-250,250,800];
    fscorr = 1; % correction for finite number of spikes per trial
    sfreq = 1000;
    trial_LB = 1500;
    num_chans = length(good_chans); 
    
freq_str = 'Theta';
% freq_band = [4,8];

params.Fs = sfreq;
params.tapers = [3,5];
% params.trialave = 1;
params.trialave = 0;

params.fpass = [2 10];
% num_freq = 4; % changes with tapers
f_select = 2; % 4-8 Hz, assume 4Hz bandwidth;

num_eval = length(eval_labels);
num_files = length(file_tags);
for iFile = 1:num_files
    curr_file = file_names{iFile};
    curr_tag = file_tags{iFile};
    disp(['Processing ',curr_tag,'...']);
    
    load(sprintf('%s\\%s_CueFrameTrialData_LB1500_RB5500_denoise.mat',loaddir2,curr_file),'LFP_trial_data','trial_events');
    load(sprintf('%s\\%s_TrialInfo_split.mat',loaddir2,curr_file),'TrueCorrect','Rule');
%     load(sprintf('%s\\new_area_codes_July2021.mat',loaddir),'area_code','good_chans');
    load(sprintf('%s\\%s\\Neuron.mat',spike_dir,curr_file),'Spk','Spk_channel');
    
    Correct_idx = find(TrueCorrect==1 & ~isnan(Rule));
    trial_events = trial_events(Correct_idx,:);
    num_trials = length(Correct_idx);
    LFP_trial_data = LFP_trial_data(good_chans,:,Correct_idx);
    
% for Chronux (coherencycpt.m)
%       data1        (continuous data in time x trials form) -- required
%       data2        (structure array of spike times with dimension trials; 
%                     also accepts 1d array of spike times) -- required
    for iEval = 1:num_eval
        event_str = eval_labels{iEval};
        if params.trialave == 0
            if exist(sprintf('%s\\%s_%s_%s_SFC_trialFull_temp.mat',savedir,curr_tag,event_str,freq_str),'file')
                load(sprintf('%s\\%s_%s_%s_SFC_trialFull_temp.mat',savedir,curr_tag,event_str,freq_str),'C_set','iBin');
                start_bin = iBin+1;
            else
                start_bin = 1;
                C_set = nan(num_chans,num_chans,num_trials,num_ts);
            end
        else
            if exist(sprintf('%s\\%s_%s_%s_SFC_trialAve_temp.mat',savedir,curr_tag,event_str,freq_str),'file')
                load(sprintf('%s\\%s_%s_%s_SFC_trialAve_temp.mat',savedir,curr_tag,event_str,freq_str),'C_set','iBin')
                start_bin = iBin+1;
            else
                start_bin = 1;
                C_set = nan(num_chans,num_chans,num_ts);
            end
        end

        for iBin = start_bin:num_ts
            tic;
            for iCL = 1:num_chans
                event_var = round((trial_events(:,event_col(iEval)) - trial_events(:,1)).*sfreq)+trial_LB-img_bounds_set{iEval}(1) + bin_centers(iBin); 
                event_ts = trial_events(:,event_col(iEval))-((img_bounds_set{iEval}(1)+bin_centers(iBin))/sfreq);
            
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
                        C_set(iCL,iChan,:,iBin) = C(f_select,:); % 4-8 Hz, assume 4Hz bandwidth;
                    else
                        C_set(iCL,iChan,iBin) = C(f_select);
                    end
                end
            end
            if params.trialave == 0
                save(sprintf('%s\\%s_%s_%s_SFC_trialFull_temp.mat',savedir,curr_tag,event_str,freq_str),'C_set','iBin','-v7.3')
            else
                save(sprintf('%s\\%s_%s_%s_SFC_trialAve_temp.mat',savedir,curr_tag,event_str,freq_str),'C_set','iBin','-v7.3')
            end
            toc;
            disp(['Time Bin ',num2str(iBin),' of ',num2str(num_ts),' processed']);
        end
        disp(['Event Trig ',num2str(iEval),' Done'])
        disp('Saving...')   
        if params.trialave == 0
            save(sprintf('%s\\%s_%s_%s_SFC_trialFull.mat',savedir,curr_tag,event_str,freq_str),'C_set','-v7.3')
        else
            save(sprintf('%s\\%s_%s_%s_SFC_trialAve.mat',savedir,curr_tag,event_str,freq_str),'C_set','-v7.3')
        end
        
        disp('Done.')
    end
end

%% SFC cross chan, sliding over time, Groot 2-item
addpath(genpath('D:\Scripts\chronux_2_12.v03\chronux_2_12\chronux_2_12'));

file_names = {'groot2024-04-15','groot2024-05-07','groot2024-05-09',...
    'groot2024-05-13','groot2024-05-14','groot2024-05-15',...
    'groot2024-05-16','groot2024-05-20','groot2024-05-21','groot2024-05-22'};
file_tags = {'Groot0415','Groot0507','Groot0509',...
    'Groot0513','Groot0514','Groot0515',...
    'Groot0516','Groot0520','Groot0521',...
    'Groot0522'};

spike_dir = 'E:\Groot_switch\sorting';
loaddir = 'E:\Groot_switch_Apr';
loaddir2 = 'E:\Groot_switch_Apr\preprocessed\trial_cuts';
savedir = 'E:\Groot_switch_Apr\preprocessed\prelim_analyses\SFC_CrossChan_Slide';
if ~exist(savedir,'dir')
    mkdir(savedir);
end

load(sprintf('%s\\new_area_codes_Mar2024.mat',loaddir),'good_chans')

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

eval_labels = {'T1_Trig'};
img_bounds_set = {[250,1000,1300]};
event_col = [2]; % T1/T2 trig for Encoding
bin_centers = 0:50:1550;


% eval_labels = {'Base'};
% img_bounds_set = {[750,1000,250]};
% event_col = [1]; % CueFrame for Baseline
% bin_centers = 0:50:500;

num_ts = length(bin_centers);
eval_size = 250; %+-X ms
% eval_starts = [-1000,-250,250,800];
    fscorr = 1; % correction for finite number of spikes per trial
    sfreq = 1000;
    trial_LB = 1500;
    num_chans = length(good_chans); 
    
freq_str = 'Theta';
% freq_band = [4,8];

params.Fs = sfreq;
params.tapers = [3,5];
params.trialave = 1;
% params.trialave = 0;

params.fpass = [2 10];
% num_freq = 4; % changes with tapers
f_select = 2; % 4-8 Hz, assume 4Hz bandwidth;

num_eval = length(eval_labels);
num_files = length(file_tags);
for iFile = 3:num_files
    curr_file = file_names{iFile};
    curr_tag = file_tags{iFile};
    disp(['Processing ',curr_tag,'...']);
    
    load(sprintf('%s\\%s_CueFrameTrialData_LB1500_RB5500_denoise.mat',loaddir2,curr_file),'LFP_trial_data','trial_events');
    load(sprintf('%s\\%s_TrialInfo_split.mat',loaddir2,curr_file),'TrueCorrect','Rule','CueFrameOnTime');
%     load(sprintf('%s\\new_area_codes_July2021.mat',loaddir),'area_code','good_chans');
    load(sprintf('%s\\%s\\Neuron.mat',spike_dir,curr_file),'Spk','Spk_channel');
    
    Correct_idx = find(TrueCorrect==1 & ~isnan(Rule) & ~isnan(CueFrameOnTime));
    trial_events = trial_events(Correct_idx,:);
    num_trials = length(Correct_idx);
    LFP_trial_data = LFP_trial_data(good_chans,:,Correct_idx);
    
% for Chronux (coherencycpt.m)
%       data1        (continuous data in time x trials form) -- required
%       data2        (structure array of spike times with dimension trials; 
%                     also accepts 1d array of spike times) -- required
    for iEval = 1:num_eval
        event_str = eval_labels{iEval};
        if params.trialave == 0
            if exist(sprintf('%s\\%s_%s_%s_SFC_trialFull_temp.mat',savedir,curr_tag,event_str,freq_str),'file')
                load(sprintf('%s\\%s_%s_%s_SFC_trialFull_temp.mat',savedir,curr_tag,event_str,freq_str),'C_set','iBin');
                start_bin = iBin+1;
            else
                start_bin = 1;
                C_set = nan(num_chans,num_chans,num_trials,num_ts);
            end
        else
            if exist(sprintf('%s\\%s_%s_%s_SFC_trialAve_temp.mat',savedir,curr_tag,event_str,freq_str),'file')
                load(sprintf('%s\\%s_%s_%s_SFC_trialAve_temp.mat',savedir,curr_tag,event_str,freq_str),'C_set','iBin')
                start_bin = iBin+1;
            else
                start_bin = 1;
                C_set = nan(num_chans,num_chans,num_ts);
            end
        end

        for iBin = start_bin:num_ts
            tic;
            for iCL = 1:num_chans
                event_var = round((trial_events(:,event_col(iEval)) - trial_events(:,1)).*sfreq)+trial_LB-img_bounds_set{iEval}(1) + bin_centers(iBin); 
                event_ts = trial_events(:,event_col(iEval))-((img_bounds_set{iEval}(1)+bin_centers(iBin))/sfreq);
            
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
                        C_set(iCL,iChan,:,iBin) = C(f_select,:); % 4-8 Hz, assume 4Hz bandwidth;
                    else
                        C_set(iCL,iChan,iBin) = C(f_select);
                    end
                end
            end
            if params.trialave == 0
                save(sprintf('%s\\%s_%s_%s_SFC_trialFull_temp.mat',savedir,curr_tag,event_str,freq_str),'C_set','iBin','-v7.3')
            else
                save(sprintf('%s\\%s_%s_%s_SFC_trialAve_temp.mat',savedir,curr_tag,event_str,freq_str),'C_set','iBin','-v7.3')
            end
            toc;
            disp(['Time Bin ',num2str(iBin),' of ',num2str(num_ts),' processed']);
        end
        disp(['Event Trig ',num2str(iEval),' Done'])
        disp('Saving...')   
        if params.trialave == 0
            save(sprintf('%s\\%s_%s_%s_SFC_trialFull.mat',savedir,curr_tag,event_str,freq_str),'C_set','-v7.3')
        else
            save(sprintf('%s\\%s_%s_%s_SFC_trialAve.mat',savedir,curr_tag,event_str,freq_str),'C_set','-v7.3')
        end
        
        disp('Done.')
    end
end

%% SFC, cross chan, search over freqs, Ocean 2-item
addpath(genpath('C:\Users\DELL\Documents\scripts\chronux_2_12.v03\chronux_2_12\chronux_2_12'));

file_names = {'ocean2021-07-09','ocean2021-07-10','ocean2021-07-13','ocean2021-07-14',...
    'ocean2021-07-15','ocean2021-07-20','ocean2021-07-21','ocean2021-07-22',...
    'ocean2021-07-23','ocean2021-07-26','ocean2021-07-27','ocean2021-07-29','ocean2021-07-30'};
file_tags = {'Ocean0709','Ocean0710','Ocean0713','Ocean0714','Ocean0715',...
    'Ocean0720','Ocean0721','Ocean0722','Ocean0723','Ocean0726',...
    'Ocean0727','Ocean0729','Ocean0730'};

spike_dir = 'G:\FW_data\NewJuly\ocean_data\sorting';
loaddir = 'G:\FW_data\NewJuly\ocean_data\preprocessed';

savedir = 'G:\FW_data\NewJuly\ocean_data\preprocessed\prelim_analyses\SFC_CrossChan_WB'; % Wide Band

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


eval_labels = {'Base','T1','T2'};
img_size = 500;
event_col = [1,2,4]; % pre-CueFrame for Baseline
img_shift_set = [-500,250,250,250,250];
num_ts = length(eval_labels);

eval_size = 250; %+-X ms
    fscorr = 1; % correction for finite number of spikes per trial
    sfreq = 1000;
    trial_LB = 1500;
    num_chans = length(good_chans); 
    
freq_str = 'WB';

params.Fs = sfreq;
params.tapers = [3,5];
params.trialave = 1;

% for WB
params.fpass = [2 150];

num_eval = length(eval_labels);
num_files = length(file_tags);
for iFile = 1:num_files

    curr_file = file_names{iFile};
    curr_tag = file_tags{iFile};
    disp(['Processing ',curr_tag,'...']);
    
    load(sprintf('%s\\trial_cuts\\%s_CueFrameTrialData_LB1500_RB5500_denoise.mat',loaddir,curr_file),'LFP_trial_data','trial_events');
    load(sprintf('%s\\trial_cuts\\%s_TrialInfo_split.mat',loaddir,curr_file),'TrueCorrect','Rule');
    load(sprintf('%s\\%s\\Neuron.mat',spike_dir,curr_file),'Spk','Spk_channel');
    
    Correct_idx = find(TrueCorrect==1 & ~isnan(Rule));
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
                C_set = cell(num_chans,num_chans,num_trials,num_ts);
            end
        else
            if exist(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve_temp.mat',savedir,curr_tag,freq_str),'file')
                load(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve_temp.mat',savedir,curr_tag,freq_str),'C_set','iEval')
                start_bin = iEval+1;
            else
                start_bin = 1;
                C_set = cell(num_chans,num_chans,num_ts);
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
                        C_set{iCL,iChan,:,iEval} = C; % 4-8 Hz, assume 4Hz bandwidth;
                    else
                        C_set{iCL,iChan,iEval} = C;
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

%% SFC, cross chan, search over freqs, Groot 2-item
addpath(genpath('C:\Users\DELL\Documents\scripts\chronux_2_12.v03\chronux_2_12\chronux_2_12'));

file_names = {'groot2024-04-15','groot2024-05-07','groot2024-05-09',...
    'groot2024-05-13','groot2024-05-14','groot2024-05-15',...
    'groot2024-05-16','groot2024-05-20','groot2024-05-21',...
    'groot2024-05-22'};
file_tags = {'Groot0415','Groot0507','Groot0509',...
    'Groot0513','Groot0514','Groot0515',...
    'Groot0516','Groot0520','Groot0521',...
    'Groot0522'};

spike_dir = 'E:\ZC_data\groot_data_switch\Groot_switch_sorting\sorting';
loaddir = 'E:\ZC_data\groot_data_switch';

savedir = 'E:\ZC_data\groot_data_switch\prelim_analyses\SFC_CrossChan_WB_new'; % Wide Band

if ~exist(savedir,'dir')
    mkdir(savedir);
end

load(sprintf('%s\\new_area_codes_Mar2024.mat',loaddir),'good_chans')

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


eval_labels = {'Base','T1','T2'};
img_size = 500;
event_col = [1,2,4]; % pre-CueFrame for Baseline
img_shift_set = [-500,250,250,250,250];
num_ts = length(eval_labels);

eval_size = 250; %+-X ms
    fscorr = 1; % correction for finite number of spikes per trial
    sfreq = 1000;
    trial_LB = 1500;
    num_chans = length(good_chans); 
    
freq_str = 'WB';

params.Fs = sfreq;
params.tapers = [3,5];
params.trialave = 1;

% for WB
params.fpass = [2 150];

num_eval = length(eval_labels);
num_files = length(file_tags);
for iFile = 1:num_files

    curr_file = file_names{iFile};
    curr_tag = file_tags{iFile};
    disp(['Processing ',curr_tag,'...']);
    
    load(sprintf('%s\\trial_cuts\\%s_CueFrameTrialData_LB1500_RB5500_denoise.mat',loaddir,curr_file),'LFP_trial_data','trial_events');
    load(sprintf('%s\\trial_cuts\\%s_TrialInfo_split.mat',loaddir,curr_file),'TrueCorrect','Rule','CueFrameOnTime');
    load(sprintf('%s\\%s\\Neuron.mat',spike_dir,curr_file),'Spk','Spk_channel');
    Spk_channel = Spk_channel(:,1);
    
    Correct_idx = find(TrueCorrect==1 & ~isnan(Rule) & ~isnan(CueFrameOnTime));
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
                C_set = cell(num_chans,num_chans,num_trials,num_ts);
            end
        else
            if exist(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve_temp.mat',savedir,curr_tag,freq_str),'file')
                load(sprintf('%s\\%s_%s_TimeSeg_SFC_trialAve_temp.mat',savedir,curr_tag,freq_str),'C_set','iEval')
                start_bin = iEval+1;
            else
                start_bin = 1;
                C_set = cell(num_chans,num_chans,num_ts);
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
                        C_set{iCL,iChan,:,iEval} = C; % 4-8 Hz, assume 4Hz bandwidth;
                    else
                        C_set{iCL,iChan,iEval} = C;
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

%% Frontal-u SFC vectors, TimeSeg500, no mask, Ocean 2-item
SFC_dir = 'G:\FW_data\NewJuly\ocean_data\preprocessed\prelim_analyses\SFC_CrossChan';
file_tags = {'0709','0710','0713','0714','0715','0720','0721','0722','0723','0726','0727','0729','0730'};
savedir = 'G:\FW_data\NewJuly\ocean_data\preprocessed\prelim_analyses\FW_analyses\Frontal';
spike_dir = 'G:\FW_data\NewJuly\ocean_data\sorting';

% idx_dir = 'C:\Users\DELL\Documents\FW_data\sharing\Unit_PEV\NewJuly\channel_info'; % uPFC
idx_dir = 'C:\Users\DELL\Documents\FW_data\sharing\Unit_PEV\NewJuly\spk_ch_info'; % Frontal

load('C:\Users\DELL\Documents\FW_data\ocean_data\preprocessed\new_area_codes_July2021.mat','good_chans','area_code');

% PM_idx = find(area_code==3);
% PFC_idx = intersect(setdiff(good_chans,12),find(area_code==2));
F_idx = intersect(setdiff(good_chans,12),find(area_code==2 | area_code==3));

SFC_T1r_Comb = []; SFC_T2r_Comb = []; SFC_T1c_Comb = []; SFC_T2c_Comb = [];
SFC_T1r_CombZ = []; SFC_T2r_CombZ = []; SFC_T1c_CombZ = []; SFC_T2c_CombZ = [];
SFC_EntR_Comb = []; SFC_EntC_Comb = []; SFC_EntR_CombZ = []; SFC_EntC_CombZ = [];
SFC_T1d_Comb = []; SFC_T2d_Comb = []; SFC_EntD_Comb = [];
SFC_T1d_CombZ = []; SFC_T2d_CombZ = []; SFC_EntD_CombZ = [];
num_files = length(file_tags);
for iFile = 1:num_files
    curr_tag = file_tags{iFile};
    
%     load(sprintf('%s\\Ocean%s_Theta_TimeSeg_SFC_trialFull.mat',SFC_dir,curr_tag),'C_set')
    load(sprintf('%s\\Ocean%s_Theta_TimeSeg_SFC_trialAve.mat',SFC_dir,curr_tag),'C_set')

    load(sprintf('%s\\ocean2021-07-%s\\Neuron.mat',spike_dir,curr_tag(end-1:end)),'Spk_channel')
    
    load(sprintf('%s\\%s_Frontal_idx.mat',idx_dir,curr_tag),'chIn256');
    
    [~, uniqueIdx] = unique(Spk_channel);
    dupeIdx = ismember( Spk_channel, Spk_channel( setdiff( 1:numel(Spk_channel), uniqueIdx ) ) );
    dupes = unique(Spk_channel(dupeIdx)); dupe_locs = find(dupeIdx);
    Spk_chan_trim = setdiff(chIn256,dupes);
    [~,uPFC_sel] = intersect(good_chans,Spk_chan_trim);
    
    mat_B_mean = squeeze(C_set(uPFC_sel,uPFC_sel,1));
    mat_T1_mean = squeeze(C_set(uPFC_sel,uPFC_sel,2));
    mat_T2_mean = squeeze(C_set(uPFC_sel,uPFC_sel,3));
    
%     T1_eval = (mat_T1_mean - mat_B_mean) > 0;
%     T2_eval = (mat_T2_mean - mat_B_mean) > 0;
    
%     T1_eval_adj = (T1_eval - T2_eval) > 0;
%     T2_eval_adj = (T2_eval - T1_eval) > 0;
    
%     mat_T1_mean = squeeze(mean(mat_T1,3));
%     mat_T2_mean = squeeze(mean(mat_T2,3));
%     mat_B_mean = squeeze(mean(mat_B,3));

    T1_mat_trim = mat_T1_mean;
    T2_mat_trim = mat_T2_mean;
    Ent_mat_trim = mat_T1_mean.*mat_T2_mean;
    T1_diag = diag(T1_mat_trim); T2_diag = diag(T2_mat_trim);
    Ent_diag = diag(Ent_mat_trim);
    
    row_vec_T1 = sum(T1_mat_trim,2); row_vec_T2 = sum(T2_mat_trim,2);
    col_vec_T1 = sum(T1_mat_trim,1); col_vec_T2 = sum(T2_mat_trim,1);
    ent_vec_r = sum(Ent_mat_trim,2); ent_vec_c = sum(Ent_mat_trim,1);
    
    SFC_T1r_Comb = [SFC_T1r_Comb;row_vec_T1]; SFC_T2r_Comb = [SFC_T2r_Comb;row_vec_T2]; 
    SFC_T1c_Comb = [SFC_T1c_Comb;col_vec_T1']; SFC_T2c_Comb = [SFC_T2c_Comb;col_vec_T2'];
    SFC_EntR_Comb = [SFC_EntR_Comb;ent_vec_r]; SFC_EntC_Comb = [SFC_EntC_Comb;ent_vec_c'];
    SFC_T1d_Comb = [SFC_T1d_Comb;T1_diag]; SFC_T2d_Comb = [SFC_T2d_Comb;T2_diag]; 
    SFC_EntD_Comb = [SFC_EntD_Comb;Ent_diag];

%     SFC_T1r_CombZ = [SFC_T1r_CombZ;anscombe(row_vec_T1)]; SFC_T2r_CombZ = [SFC_T2r_CombZ;anscombe(row_vec_T2)]; 
%     SFC_T1c_CombZ = [SFC_T1c_CombZ;anscombe(col_vec_T1')]; SFC_T2c_CombZ = [SFC_T2c_CombZ;anscombe(col_vec_T2')];
%     SFC_EntR_CombZ = [SFC_EntR_CombZ;anscombe(ent_vec_r)]; SFC_EntC_CombZ = [SFC_EntC_CombZ;anscombe(ent_vec_c')];
%     SFC_T1d_CombZ = [SFC_T1d_CombZ;anscombe(T1_diag)]; SFC_T2d_CombZ = [SFC_T2d_CombZ;anscombe(T2_diag)]; 
%     SFC_EntD_CombZ = [SFC_EntD_CombZ;anscombe(Ent_diag)];

    SFC_T1r_CombZ = [SFC_T1r_CombZ;zscore(row_vec_T1)]; SFC_T2r_CombZ = [SFC_T2r_CombZ;zscore(row_vec_T2)]; 
    SFC_T1c_CombZ = [SFC_T1c_CombZ;zscore(col_vec_T1')]; SFC_T2c_CombZ = [SFC_T2c_CombZ;zscore(col_vec_T2')];
    SFC_EntR_CombZ = [SFC_EntR_CombZ;zscore(ent_vec_r)]; SFC_EntC_CombZ = [SFC_EntC_CombZ;zscore(ent_vec_c')];
    SFC_T1d_CombZ = [SFC_T1d_CombZ;zscore(T1_diag)]; SFC_T2d_CombZ = [SFC_T2d_CombZ;zscore(T2_diag)]; 
    SFC_EntD_CombZ = [SFC_EntD_CombZ;zscore(Ent_diag)];
end

save(sprintf('%s\\Ocean13Days_SFCvec_concat_BaseComp_NoSubtract_ZC_AddEnt_TrueIdxVerify_Frontal_NoMask.mat',savedir),'*_Comb','*_CombZ');

% save(sprintf('%s\\Ocean13Days_SFCvec_concat_BaseComp_NoSubtract_AddEnt_TrueIdxVerify_Frontal_NoMask.mat',savedir),'*_Comb','*_CombZ');
% save(sprintf('%s\\Ocean13Days_SFCvec_concat_BaseComp_NoSubtract_AddEnt_TrueIdxVerify_uPFC.mat',savedir),'*_Comb','*_CombZ');

%% Beta SFC x FW: get SFC strength vectors first, Ocean
freq_sel = 10:16;
SFC_dir = 'G:\FW_data\NewJuly\ocean_data\preprocessed\prelim_analyses\SFC_CrossChan_WB';

file_tags = {'0709','0710','0713','0714','0715','0720','0721','0722','0723','0726','0727','0729','0730'};
savedir = 'G:\PaperPrep\CurrBiolRevision';
spike_dir = 'G:\FW_data\NewJuly\ocean_data\sorting';

% idx_dir = 'C:\Users\DELL\Documents\FW_data\sharing\Unit_PEV\NewJuly\channel_info'; % uPFC
idx_dir = 'C:\Users\DELL\Documents\FW_data\sharing\Unit_PEV\NewJuly\spk_ch_info'; % Frontal

load('C:\Users\DELL\Documents\FW_data\ocean_data\preprocessed\new_area_codes_July2021.mat','good_chans','area_code');

% PM_idx = find(area_code==3);
% PFC_idx = intersect(setdiff(good_chans,12),find(area_code==2));
F_idx = intersect(setdiff(good_chans,12),find(area_code==2 | area_code==3));

SFC_T1r_Comb = []; SFC_T2r_Comb = []; SFC_T1c_Comb = []; SFC_T2c_Comb = [];
SFC_T1r_CombZ = []; SFC_T2r_CombZ = []; SFC_T1c_CombZ = []; SFC_T2c_CombZ = [];
SFC_EntR_Comb = []; SFC_EntC_Comb = []; SFC_EntR_CombZ = []; SFC_EntC_CombZ = [];
SFC_T1d_Comb = []; SFC_T2d_Comb = []; SFC_EntD_Comb = [];
SFC_T1d_CombZ = []; SFC_T2d_CombZ = []; SFC_EntD_CombZ = [];
num_files = length(file_tags);
for iFile = 1:num_files
    curr_tag = file_tags{iFile};
    
    load(sprintf('%s\\Ocean%s_WB_TimeSeg_SFC_trialAve.mat',SFC_dir,curr_tag),'C_set')
    C_set_full = C_set;
    C_set = nan(size(C_set_full));
    for iC1 = 1:size(C_set,1)
        for iC2 = 1:size(C_set,1)
            for iTime = 1:3
                if ~isempty(C_set_full{iC1,iC2,iTime})
                    C_set(iC1,iC2,iTime) = nanmean(C_set_full{iC1,iC2,iTime}(freq_sel));
                end
            end
        end
    end

    load(sprintf('%s\\ocean2021-07-%s\\Neuron.mat',spike_dir,curr_tag(end-1:end)),'Spk_channel')
    
    load(sprintf('%s\\%s_Frontal_idx.mat',idx_dir,curr_tag),'chIn256');
    
    [~, uniqueIdx] = unique(Spk_channel);
    dupeIdx = ismember( Spk_channel, Spk_channel( setdiff( 1:numel(Spk_channel), uniqueIdx ) ) );
    dupes = unique(Spk_channel(dupeIdx)); dupe_locs = find(dupeIdx);
    Spk_chan_trim = setdiff(chIn256,dupes);
    [~,uPFC_sel] = intersect(good_chans,Spk_chan_trim);
    
    mat_B_mean = squeeze(C_set(uPFC_sel,uPFC_sel,1));
    mat_T1_mean = squeeze(C_set(uPFC_sel,uPFC_sel,2));
    mat_T2_mean = squeeze(C_set(uPFC_sel,uPFC_sel,3));
    
%     T1_eval = (mat_T1_mean - mat_B_mean) > 0;
%     T2_eval = (mat_T2_mean - mat_B_mean) > 0;
    
%     T1_eval_adj = (T1_eval - T2_eval) > 0;
%     T2_eval_adj = (T2_eval - T1_eval) > 0;
    
%     mat_T1_mean = squeeze(mean(mat_T1,3));
%     mat_T2_mean = squeeze(mean(mat_T2,3));
%     mat_B_mean = squeeze(mean(mat_B,3));

    T1_mat_trim = mat_T1_mean;
    T2_mat_trim = mat_T2_mean;
    Ent_mat_trim = mat_T1_mean.*mat_T2_mean;
    T1_diag = diag(T1_mat_trim); T2_diag = diag(T2_mat_trim);
    Ent_diag = diag(Ent_mat_trim);
    
    row_vec_T1 = sum(T1_mat_trim,2); row_vec_T2 = sum(T2_mat_trim,2);
    col_vec_T1 = sum(T1_mat_trim,1); col_vec_T2 = sum(T2_mat_trim,1);
    ent_vec_r = sum(Ent_mat_trim,2); ent_vec_c = sum(Ent_mat_trim,1);
    
    SFC_T1r_Comb = [SFC_T1r_Comb;row_vec_T1]; SFC_T2r_Comb = [SFC_T2r_Comb;row_vec_T2]; 
    SFC_T1c_Comb = [SFC_T1c_Comb;col_vec_T1']; SFC_T2c_Comb = [SFC_T2c_Comb;col_vec_T2'];
    SFC_EntR_Comb = [SFC_EntR_Comb;ent_vec_r]; SFC_EntC_Comb = [SFC_EntC_Comb;ent_vec_c'];
    SFC_T1d_Comb = [SFC_T1d_Comb;T1_diag]; SFC_T2d_Comb = [SFC_T2d_Comb;T2_diag]; 
    SFC_EntD_Comb = [SFC_EntD_Comb;Ent_diag];

    SFC_T1r_CombZ = [SFC_T1r_CombZ;zscore(row_vec_T1)]; SFC_T2r_CombZ = [SFC_T2r_CombZ;zscore(row_vec_T2)]; 
    SFC_T1c_CombZ = [SFC_T1c_CombZ;zscore(col_vec_T1')]; SFC_T2c_CombZ = [SFC_T2c_CombZ;zscore(col_vec_T2')];
    SFC_EntR_CombZ = [SFC_EntR_CombZ;zscore(ent_vec_r)]; SFC_EntC_CombZ = [SFC_EntC_CombZ;zscore(ent_vec_c')];
    SFC_T1d_CombZ = [SFC_T1d_CombZ;zscore(T1_diag)]; SFC_T2d_CombZ = [SFC_T2d_CombZ;zscore(T2_diag)]; 
    SFC_EntD_CombZ = [SFC_EntD_CombZ;zscore(Ent_diag)];
end

save(sprintf('%s\\Ocean13Days_BetaSFCvec_concat.mat',savedir),'*_Comb','*_CombZ');

%% Beta SFC x FW: get SFC strength vectors first, Groot
freq_sel = 10:16;
SFC_dir = 'E:\ZC_data\groot_data_switch\prelim_analyses\SFC_CrossChan_WB_new';

file_tags = {'0415','0507','0509','0513','0514','0515','0516','0520','0521','0522'};
savedir = 'G:\PaperPrep\CurrBiolRevision';
spike_dir = 'E:\ZC_data\groot_data_switch\Groot_switch_sorting\sorting';

idx_dir = 'E:\ZC_data\groot_data_switch\sharing\groot_switch_ch_idx';

load('C:\Users\DELL\Documents\ZC_data\groot_data\new_area_codes_Mar2024.mat','good_chans','area_code');

SFC_T1r_Comb = []; SFC_T2r_Comb = []; SFC_T1c_Comb = []; SFC_T2c_Comb = [];
SFC_T1r_CombZ = []; SFC_T2r_CombZ = []; SFC_T1c_CombZ = []; SFC_T2c_CombZ = [];
SFC_EntR_Comb = []; SFC_EntC_Comb = []; SFC_EntR_CombZ = []; SFC_EntC_CombZ = [];
SFC_T1d_Comb = []; SFC_T2d_Comb = []; SFC_EntD_Comb = [];
SFC_T1d_CombZ = []; SFC_T2d_CombZ = []; SFC_EntD_CombZ = [];
num_files = length(file_tags);
for iFile = 1:num_files
    curr_tag = file_tags{iFile};
    
    load(sprintf('%s\\Groot%s_WB_TimeSeg_SFC_trialAve.mat',SFC_dir,curr_tag),'C_set')
    C_set_full = C_set;
    C_set = nan(size(C_set_full));
    for iC1 = 1:size(C_set,1)
        for iC2 = 1:size(C_set,1)
            for iTime = 1:3
                if ~isempty(C_set_full{iC1,iC2,iTime})
                    C_set(iC1,iC2,iTime) = nanmean(C_set_full{iC1,iC2,iTime}(freq_sel));
                end
            end
        end
    end

    load(sprintf('%s\\groot2024-%s-%s\\Neuron.mat',spike_dir,curr_tag(1:2),curr_tag(end-1:end)),'Spk_channel')
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
    
%     T1_eval = (mat_T1_mean - mat_B_mean) > 0;
%     T2_eval = (mat_T2_mean - mat_B_mean) > 0;
    
%     T1_eval_adj = (T1_eval - T2_eval) > 0;
%     T2_eval_adj = (T2_eval - T1_eval) > 0;
    
%     mat_T1_mean = squeeze(mean(mat_T1,3));
%     mat_T2_mean = squeeze(mean(mat_T2,3));
%     mat_B_mean = squeeze(mean(mat_B,3));

    T1_mat_trim = mat_T1_mean;
    T2_mat_trim = mat_T2_mean;
    Ent_mat_trim = mat_T1_mean.*mat_T2_mean;
    T1_diag = diag(T1_mat_trim); T2_diag = diag(T2_mat_trim);
    Ent_diag = diag(Ent_mat_trim);
    
    row_vec_T1 = sum(T1_mat_trim,2); row_vec_T2 = sum(T2_mat_trim,2);
    col_vec_T1 = sum(T1_mat_trim,1); col_vec_T2 = sum(T2_mat_trim,1);
    ent_vec_r = sum(Ent_mat_trim,2); ent_vec_c = sum(Ent_mat_trim,1);
    
    SFC_T1r_Comb = [SFC_T1r_Comb;row_vec_T1]; SFC_T2r_Comb = [SFC_T2r_Comb;row_vec_T2]; 
    SFC_T1c_Comb = [SFC_T1c_Comb;col_vec_T1']; SFC_T2c_Comb = [SFC_T2c_Comb;col_vec_T2'];
    SFC_EntR_Comb = [SFC_EntR_Comb;ent_vec_r]; SFC_EntC_Comb = [SFC_EntC_Comb;ent_vec_c'];
    SFC_T1d_Comb = [SFC_T1d_Comb;T1_diag]; SFC_T2d_Comb = [SFC_T2d_Comb;T2_diag]; 
    SFC_EntD_Comb = [SFC_EntD_Comb;Ent_diag];

    SFC_T1r_CombZ = [SFC_T1r_CombZ;zscore(row_vec_T1)]; SFC_T2r_CombZ = [SFC_T2r_CombZ;zscore(row_vec_T2)]; 
    SFC_T1c_CombZ = [SFC_T1c_CombZ;zscore(col_vec_T1')]; SFC_T2c_CombZ = [SFC_T2c_CombZ;zscore(col_vec_T2')];
    SFC_EntR_CombZ = [SFC_EntR_CombZ;zscore(ent_vec_r)]; SFC_EntC_CombZ = [SFC_EntC_CombZ;zscore(ent_vec_c')];
    SFC_T1d_CombZ = [SFC_T1d_CombZ;zscore(T1_diag)]; SFC_T2d_CombZ = [SFC_T2d_CombZ;zscore(T2_diag)]; 
    SFC_EntD_CombZ = [SFC_EntD_CombZ;zscore(Ent_diag)];
end

save(sprintf('%s\\Groot10Days_BetaSFCvec_concat.mat',savedir),'*_Comb','*_CombZ');

%% SFC phase distrib. trimmed via ANOVA, two monkeys
Ocean_SFCP_dir = 'G:\FW_data\NewJuly\ocean_data\preprocessed\prelim_analyses\SFC_CrossChan_WithPhase';
Ocean_tags = {'Ocean0709','Ocean0710','Ocean0713','Ocean0714','Ocean0715',...
    'Ocean0720','Ocean0721','Ocean0722','Ocean0723','Ocean0726',...
    'Ocean0727','Ocean0729','Ocean0730'};
spike_dir_O = 'G:\FW_data\NewJuly\ocean_data\sorting';
idx_dir_O = 'C:\Users\DELL\Documents\FW_data\sharing\Unit_PEV\NewJuly\spk_ch_info';
load('G:\PaperPrep\PostNeuronRevision\Ocean13Days_ANOVAvec_concat_Frontal.mat','ANOVA_vec_SPK_D1');
ANOVA_vec_O = ANOVA_vec_SPK_D1;
load('G:\PaperPrep\PostNeuronRevision\Ocean13Days_ROI_idx_concat.mat','day_idx');
day_idx_O = day_idx;
load('C:\Users\DELL\Documents\FW_data\ocean_data\preprocessed\new_area_codes_July2021.mat','good_chans','area_code');
F_idx_O = intersect(setdiff(good_chans,12),find(area_code==2 | area_code==3));
[~,usable_idx_O] = intersect(good_chans,F_idx_O);

num_files_O = length(Ocean_tags);
phase_T1_full_O = []; phase_T1_SigMem_O = []; phase_T1_NSMem_O = []; phase_T1_S1_O = []; phase_T1_S2_O = []; phase_T1_S1O_O = []; phase_T1_S2O_O = [];
phase_T2_full_O = []; phase_T2_SigMem_O = []; phase_T2_NSMem_O = []; phase_T2_S1_O = []; phase_T2_S2_O = []; phase_T2_S1O_O = []; phase_T2_S2O_O = [];
phase_mats_O = cell(2,num_files_O); phase_vecs_O = cell(2,num_files_O,7);

Groot_SFCP_dir = 'E:\ZC_data\groot_data_switch\prelim_analyses\SFC_CrossChan_WithPhase';
Groot_tags = {'Groot0415','Groot0507','Groot0509','Groot0513','Groot0514',...
    'Groot0515','Groot0516','Groot0520','Groot0521','Groot0522'};
spike_dir_G = 'E:\ZC_data\groot_data_switch\Groot_switch_sorting\sorting';
idx_dir_G = 'E:\ZC_data\groot_data_switch\sharing\groot_switch_ch_idx';
load('G:\PaperPrep\PostNeuronRevision\Groot10Days_ANOVAvec_concat_Frontal.mat','ANOVA_vec_SPK_D1');
ANOVA_vec_G = ANOVA_vec_SPK_D1;
load('G:\PaperPrep\PostNeuronRevision\Groot10Days_ROI_idx_concat.mat','day_idx');
day_idx_G = day_idx;
load('C:\Users\DELL\Documents\ZC_data\groot_data\new_area_codes_Mar2024.mat','good_chans','area_code');
usable_idx_G = intersect(good_chans,find(area_code==2 | area_code==3));

num_files_G = length(Groot_tags);
phase_T1_full_G = []; phase_T1_SigMem_G = []; phase_T1_NSMem_G = []; phase_T1_S1_G = []; phase_T1_S2_G = []; phase_T1_S1O_G = []; phase_T1_S2O_G = [];
phase_T2_full_G = []; phase_T2_SigMem_G = []; phase_T2_NSMem_G = []; phase_T2_S1_G = []; phase_T2_S2_G = []; phase_T2_S1O_G = []; phase_T2_S2O_G = [];
phase_mats_G = cell(2,num_files_G); phase_vecs_G = cell(2,num_files_G,7);

for iFile = 1:num_files_O
    curr_tag = Ocean_tags{iFile};
    day_sel = find(day_idx_O == iFile);
    ANOVA_vec = ANOVA_vec_O(day_sel);
    ANOVA_sel_S = find(~isnan(ANOVA_vec));
    ANOVA_sel_NS = find(isnan(ANOVA_vec));
    ANOVA_sel_1 = find(ANOVA_vec == 1 | ANOVA_vec == 3);
    ANOVA_sel_2 = find(ANOVA_vec == 2 | ANOVA_vec == 3);
    ANOVA_sel_1O = find(ANOVA_vec == 1);
    ANOVA_sel_2O = find(ANOVA_vec == 2);
    
    load(sprintf('%s\\%s_Theta_TimeSeg_SFC_trialAve.mat',Ocean_SFCP_dir,curr_tag),'Pha_set');
    load(sprintf('%s\\ocean2021-07-%s\\Neuron.mat',spike_dir_O,curr_tag(end-1:end)),'Spk_channel')
    load(sprintf('%s\\%s_Frontal_idx.mat',idx_dir_O,curr_tag(end-3:end)),'chIn256');
    Pha_set = Pha_set(usable_idx_O,usable_idx_O,:);
    
    [~, uniqueIdx] = unique(Spk_channel);
    dupeIdx = ismember( Spk_channel, Spk_channel( setdiff( 1:numel(Spk_channel), uniqueIdx ) ) );
    dupes = unique(Spk_channel(dupeIdx)); dupe_locs = find(dupeIdx);

    Spk_chan_trim = setdiff(chIn256,dupes);
    [~,uF_sel] = intersect(F_idx_O,Spk_chan_trim);
    
    temp_T1_full = squeeze(Pha_set(uF_sel,uF_sel,2)); temp_T2_full = squeeze(Pha_set(uF_sel,uF_sel,3)); 
    phase_mats_O{1,iFile} = temp_T1_full;
    phase_mats_O{2,iFile} = temp_T2_full;
    temp_T1_SigMem = temp_T1_full(ANOVA_sel_S,ANOVA_sel_S); temp_T2_SigMem = temp_T2_full(ANOVA_sel_S,ANOVA_sel_S);
    temp_T1_NSMem = temp_T1_full(ANOVA_sel_NS,ANOVA_sel_NS); temp_T2_NSMem = temp_T2_full(ANOVA_sel_NS,ANOVA_sel_NS);
    temp_T1_S1 = temp_T1_full(ANOVA_sel_1,ANOVA_sel_1); temp_T2_S1 = temp_T2_full(ANOVA_sel_1,ANOVA_sel_1);
    temp_T1_S2 = temp_T1_full(ANOVA_sel_2,ANOVA_sel_2); temp_T2_S2 = temp_T2_full(ANOVA_sel_2,ANOVA_sel_2);
    temp_T1_S1O = temp_T1_full(ANOVA_sel_1O,ANOVA_sel_1O); temp_T2_S1O = temp_T2_full(ANOVA_sel_1O,ANOVA_sel_1O);
    temp_T1_S2O = temp_T1_full(ANOVA_sel_2O,ANOVA_sel_2O); temp_T2_S2O = temp_T2_full(ANOVA_sel_2O,ANOVA_sel_2O);
    
    phase_T1_full_O = [phase_T1_full_O;wrapTo2Pi(temp_T1_full(:))]; 
    phase_T1_SigMem_O = [phase_T1_SigMem_O;wrapTo2Pi(temp_T1_SigMem(:))]; 
    phase_T1_NSMem_O = [phase_T1_NSMem_O;wrapTo2Pi(temp_T1_NSMem(:))];
    phase_T1_S1_O = [phase_T1_S1_O;wrapTo2Pi(temp_T1_S1(:))];
    phase_T1_S2_O = [phase_T1_S2_O;wrapTo2Pi(temp_T1_S2(:))];
    phase_T1_S1O_O = [phase_T1_S1O_O;wrapTo2Pi(temp_T1_S1O(:))];
    phase_T1_S2O_O = [phase_T1_S2O_O;wrapTo2Pi(temp_T1_S2O(:))];

    phase_T2_full_O = [phase_T2_full_O;wrapTo2Pi(temp_T2_full(:))]; 
    phase_T2_SigMem_O = [phase_T2_SigMem_O;wrapTo2Pi(temp_T2_SigMem(:))]; 
    phase_T2_NSMem_O = [phase_T2_NSMem_O;wrapTo2Pi(temp_T2_NSMem(:))];
    phase_T2_S1_O = [phase_T2_S1_O;wrapTo2Pi(temp_T2_S1(:))];
    phase_T2_S2_O = [phase_T2_S2_O;wrapTo2Pi(temp_T2_S2(:))];
    phase_T2_S1O_O = [phase_T2_S1O_O;wrapTo2Pi(temp_T2_S1O(:))];
    phase_T2_S2O_O = [phase_T2_S2O_O;wrapTo2Pi(temp_T2_S2O(:))];
    
    phase_vecs_O{1,iFile,1} = wrapTo2Pi(temp_T1_full(:));
    phase_vecs_O{1,iFile,2} = wrapTo2Pi(temp_T1_SigMem(:));
    phase_vecs_O{1,iFile,3} = wrapTo2Pi(temp_T1_NSMem(:));
    phase_vecs_O{1,iFile,4} = wrapTo2Pi(temp_T1_S1(:));
    phase_vecs_O{1,iFile,5} = wrapTo2Pi(temp_T1_S2(:));
    phase_vecs_O{1,iFile,6} = wrapTo2Pi(temp_T1_S1O(:));
    phase_vecs_O{1,iFile,7} = wrapTo2Pi(temp_T1_S2O(:));
    
    phase_vecs_O{2,iFile,1} = wrapTo2Pi(temp_T2_full(:));
    phase_vecs_O{2,iFile,2} = wrapTo2Pi(temp_T2_SigMem(:));
    phase_vecs_O{2,iFile,3} = wrapTo2Pi(temp_T2_NSMem(:));
    phase_vecs_O{2,iFile,4} = wrapTo2Pi(temp_T2_S1(:));
    phase_vecs_O{2,iFile,5} = wrapTo2Pi(temp_T2_S2(:));
    phase_vecs_O{2,iFile,6} = wrapTo2Pi(temp_T2_S1O(:));
    phase_vecs_O{2,iFile,7} = wrapTo2Pi(temp_T2_S2O(:));
end

for iFile = 1:num_files_G
    curr_tag = Groot_tags{iFile};
    day_sel = find(day_idx_G == iFile);
    ANOVA_vec = ANOVA_vec_G(day_sel);
    ANOVA_sel_S = find(~isnan(ANOVA_vec));
    ANOVA_sel_NS = find(isnan(ANOVA_vec));
    ANOVA_sel_1 = find(ANOVA_vec == 1 | ANOVA_vec == 3);
    ANOVA_sel_2 = find(ANOVA_vec == 2 | ANOVA_vec == 3);
    ANOVA_sel_1O = find(ANOVA_vec == 1);
    ANOVA_sel_2O = find(ANOVA_vec == 2);
    
    load(sprintf('%s\\%s_Theta_TimeSeg_SFC_trialAve.mat',Groot_SFCP_dir,curr_tag),'Pha_set');
    load(sprintf('%s\\groot2024-%s-%s\\Neuron.mat',spike_dir_G,curr_tag(end-3:end-2),curr_tag(end-1:end)),'Spk_channel')
    Spk_channel = Spk_channel(:,1);
    load(sprintf('%s\\%s_ch_idx.mat',idx_dir_G,curr_tag(end-3:end)),'chIn256');
    Pha_set = Pha_set(usable_idx_G,usable_idx_G,:);
    
    [~, uniqueIdx] = unique(Spk_channel);
    dupeIdx = ismember( Spk_channel, Spk_channel( setdiff( 1:numel(Spk_channel), uniqueIdx ) ) );
    dupes = unique(Spk_channel(dupeIdx)); dupe_locs = find(dupeIdx);

    Spk_chan_trim = setdiff(chIn256,dupes);
    [~,uF_sel] = intersect(usable_idx_G,Spk_chan_trim);
    
    temp_T1_full = squeeze(Pha_set(uF_sel,uF_sel,2)); temp_T2_full = squeeze(Pha_set(uF_sel,uF_sel,3)); 
    phase_mats_G{1,iFile} = temp_T1_full;
    phase_mats_G{2,iFile} = temp_T2_full;
    
    temp_T1_SigMem = temp_T1_full(ANOVA_sel_S,ANOVA_sel_S); temp_T2_SigMem = temp_T2_full(ANOVA_sel_S,ANOVA_sel_S);
    temp_T1_NSMem = temp_T1_full(ANOVA_sel_NS,ANOVA_sel_NS); temp_T2_NSMem = temp_T2_full(ANOVA_sel_NS,ANOVA_sel_NS);
    temp_T1_S1 = temp_T1_full(ANOVA_sel_1,ANOVA_sel_1); temp_T2_S1 = temp_T2_full(ANOVA_sel_1,ANOVA_sel_1);
    temp_T1_S2 = temp_T1_full(ANOVA_sel_2,ANOVA_sel_2); temp_T2_S2 = temp_T2_full(ANOVA_sel_2,ANOVA_sel_2);
    temp_T1_S1O = temp_T1_full(ANOVA_sel_1O,ANOVA_sel_1O); temp_T2_S1O = temp_T2_full(ANOVA_sel_1O,ANOVA_sel_1O);
    temp_T1_S2O = temp_T1_full(ANOVA_sel_2O,ANOVA_sel_2O); temp_T2_S2O = temp_T2_full(ANOVA_sel_2O,ANOVA_sel_2O);
    
    phase_T1_full_G = [phase_T1_full_G;wrapTo2Pi(temp_T1_full(:))]; 
    phase_T1_SigMem_G = [phase_T1_SigMem_G;wrapTo2Pi(temp_T1_SigMem(:))]; 
    phase_T1_NSMem_G = [phase_T1_NSMem_G;wrapTo2Pi(temp_T1_NSMem(:))];
    phase_T1_S1_G = [phase_T1_S1_G;wrapTo2Pi(temp_T1_S1(:))];
    phase_T1_S2_G = [phase_T1_S2_G;wrapTo2Pi(temp_T1_S2(:))];
    phase_T1_S1O_G = [phase_T1_S1O_G;wrapTo2Pi(temp_T1_S1O(:))];
    phase_T1_S2O_G = [phase_T1_S2O_G;wrapTo2Pi(temp_T1_S2O(:))];

    phase_T2_full_G = [phase_T2_full_O;wrapTo2Pi(temp_T2_full(:))]; 
    phase_T2_SigMem_G = [phase_T2_SigMem_G;wrapTo2Pi(temp_T2_SigMem(:))]; 
    phase_T2_NSMem_G = [phase_T2_NSMem_G;wrapTo2Pi(temp_T2_NSMem(:))];
    phase_T2_S1_G = [phase_T2_S1_G;wrapTo2Pi(temp_T2_S1(:))];
    phase_T2_S2_G = [phase_T2_S2_G;wrapTo2Pi(temp_T2_S2(:))];
    phase_T2_S1O_G = [phase_T2_S1O_G;wrapTo2Pi(temp_T2_S1O(:))];
    phase_T2_S2O_G = [phase_T2_S2O_G;wrapTo2Pi(temp_T2_S2O(:))];
    
    phase_vecs_G{1,iFile,1} = wrapTo2Pi(temp_T1_full(:));
    phase_vecs_G{1,iFile,2} = wrapTo2Pi(temp_T1_SigMem(:));
    phase_vecs_G{1,iFile,3} = wrapTo2Pi(temp_T1_NSMem(:));
    phase_vecs_G{1,iFile,4} = wrapTo2Pi(temp_T1_S1(:));
    phase_vecs_G{1,iFile,5} = wrapTo2Pi(temp_T1_S2(:));
    phase_vecs_G{1,iFile,6} = wrapTo2Pi(temp_T1_S1O(:));
    phase_vecs_G{1,iFile,7} = wrapTo2Pi(temp_T1_S2O(:));
    
    phase_vecs_G{2,iFile,1} = wrapTo2Pi(temp_T2_full(:));
    phase_vecs_G{2,iFile,2} = wrapTo2Pi(temp_T2_SigMem(:));
    phase_vecs_G{2,iFile,3} = wrapTo2Pi(temp_T2_NSMem(:));
    phase_vecs_G{2,iFile,4} = wrapTo2Pi(temp_T2_S1(:));
    phase_vecs_G{2,iFile,5} = wrapTo2Pi(temp_T2_S2(:));
    phase_vecs_G{2,iFile,6} = wrapTo2Pi(temp_T2_S1O(:));
    phase_vecs_G{2,iFile,7} = wrapTo2Pi(temp_T2_S2O(:));
end

save('G:\PaperPrep\PostNeuronRevision\TwoMonkeys_SFC_Phase_vecs_by_unit_ANOVAsplit.mat','phase_*');
