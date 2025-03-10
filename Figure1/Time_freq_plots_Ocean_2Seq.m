%% TF analysis: parameter setup
% requires presence of trial variables
        params.trial_subset = 'all';
        params.cycles = [1,0.5];
        params.sfreq = 1000;
        
        params.alpha = 0.05;
        params.winsize = params.sfreq;
        params.maxfreq = 200;
        params.freqs = [2, 200];
        params.freqscale = 'log';
        params.plotphase = 'off';
        params.plotitc = 'off';
        params.padratio = 1;
        params.mcorrect = 'fdr';

        params.verbose = 'off';
        params.baseline = [-1*params.sfreq,-0.5*params.sfreq];    % CueFrameOnTime
        params.tlimits = [-1500,5500]; % actually -1s to +4s % CueFrameOnTime
        params.timewarp_flag = 1;
        
%% TF analysis (requires EEGLAB newtimef())
file_names = {'ocean2021-07-09','ocean2021-07-13','ocean2021-07-14',...
    'ocean2021-07-20','ocean2021-07-21','ocean2021-07-22',...
    'ocean2021-07-26','ocean2021-07-27','ocean2021-07-29',...
    'ocean2021-07-30'};
file_tags = {'0709','0713','0714','0720','0721','0722','0726','0727','0729','0730'};
savedir = 'G:\FW_data\NewJuly\ocean_data\preprocessed';

fig_dir = 'G:\FW_data\NewJuly\ocean_data\preprocessed\prelim_analyses\TF_plots\TFplots_tw_2-200Hz_log_SepDays_NewJuly';
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end

trial_subset_types = {'all_correct'};
title_strings = {'All Correct'}; 
sfreq = 1000; bad_chans = [144,159,160];

trial_LB = 1500; trial_RB = 5500;

cue_type = 'CueFrame';

trial_limit_flag = 0; 
trial_limit = 51;

num_files = length(file_tags);
num_types = length(trial_subset_types);

for iFile = 1:num_files
    file_name = file_names{iFile}; file_tag = file_tags{iFile};
    
    load(sprintf('%s\\trial_cuts\\%s_%sTrialData_LB%d_RB%d_denoise.mat',savedir,file_name,cue_type,trial_LB,trial_RB),'LFP_trial_data','trial_events');

    load(sprintf('%s\\chan_labels_rough.mat',savedir), 'chan_labels');
    num_chans = length(chan_labels);

    load(sprintf('%s\\trial_cuts\\%s_TrialInfo_split.mat',savedir,file_name),...
        'TrueCorrect','Rule');
    Correct_idx = find(TrueCorrect==1 & ~isnan(Rule));

% 2nd dim for events: 1, CueFrameOnTime; 2,3/4,5/6,7, Target(Off)Time; 
% 8, RuleCueTime; 9, RuleOffTime; 10, GoCueTime; 11-13, ResponseTime; 
% 14, ReleseBarTime; 15/16: TrialStart(End)Time
    LFP_events_set_cut = round((trial_events(:,[1,2,4,8,10]) - trial_events(:,1)).*sfreq);
    
    tic;
    for iType = 1:num_types
        curr_type = trial_subset_types{iType};
        curr_title = title_strings{iType};
        disp(['Initiating Type ',num2str(iType),': ',curr_type])
        temp_dir1 = sprintf('%s\\%s\\%s\\figs',fig_dir,file_tag,curr_type);
        if ~exist(temp_dir1,'dir')
            mkdir(temp_dir1)
        end
        temp_dir2 = sprintf('%s\\%s\\%s\\pngs',fig_dir,file_tag,curr_type);
        if ~exist(temp_dir2,'dir')
            mkdir(temp_dir2)
        end
        switch lower(curr_type)
            case {'all'} % has response (and GoCue)
                trim_idx = 1:size(LFP_trial_data,3);
            case {'all_correct'} % correct/error trial always has response (and GoCue)
                trim_idx = Correct_idx;
            case {'all_error'}
                trim_idx = ErrorResp_idx;
            case {'error_oe'}
                trim_idx = OE_idx;
            case {'error_ie'}
                trim_idx = IE_idx;
        end
        params.timewarp = LFP_events_set_cut(trim_idx,:);
        
        for iChan = 1:num_chans
            if ismember(iChan,bad_chans)
                continue
            end
            %         subset = setdiff(trim_idx,find(bad_trial_idx(iChan,:)==1))';
            subset = trim_idx;
            if trial_limit_flag == 1 && length(subset) > trial_limit
                subset = subset(randi(length(subset),trial_limit)); % random subselection to balance against low trial count conditions
            end
            curr_data = LFP_trial_data(iChan,:,subset);
            
            [~,~,~] = TFplots_MonkeySeq(curr_data,params);
            suptitle([file_tag,', ',chan_labels{iChan},', ',cue_type,'-triggered, n=',num2str(length(subset)),', ',curr_title]);
            saveas(gcf,sprintf('%s\\%s_%s_TF_LB%d_RB%d_%s_%s.fig',temp_dir1,...
                file_tag,cue_type,trial_LB,trial_RB,curr_type,chan_labels{iChan}));
            saveas(gcf,sprintf('%s\\%s_%s_TF_LB%d_RB%d_%s_%s.png',temp_dir2,...
                file_tag,cue_type,trial_LB,trial_RB,curr_type,chan_labels{iChan}));
            close(gcf)
            disp(['Chan ',num2str(iChan),' finished'])
        end
    end
    toc;

end
