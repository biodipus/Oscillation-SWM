% Time-frequency analysis workhorse
% requires EEGLAB function "newtimef"
% trial_data: 2-D array (frames,trials)
function [ersp,itc,tfdata] = TFplots_MonkeySeq(trial_data,params)
    if isempty(params)
        params.trial_subset = 'all';
        params.cycles = [1,0.5];
        params.sfreq = 1000;
        
%         params.baseline = 0;
        params.baseline = [-1*params.sfreq,-0.5*params.sfreq];  % GoCue & RuleCue
        params.tlimits = [-1000,1000]; % actually -0.5s to +0.5s % GoCue & RuleCue
%         params.baseline = [-3.5*params.sfreq,-3*params.sfreq];  % GoCue
%         params.tlimits = [-3500,3500]; % actually -3s to +3s % GoCue
%         params.baseline = [-0.5*params.sfreq,0];    % CueFrameOnTime
%         params.tlimits = [-1000,5500]; % actually -0.5s to +5s % CueFrameOnTime
        
        params.alpha = 0.05;
        params.winsize = params.sfreq;
        params.maxfreq = 48;
        params.freqs = [2, 48];
        params.freqscale = 'linear';
        params.plotphase = 'off';
        params.plotitc = 'off';
        params.padratio = 1;
        params.mcorrect = 'fdr';
        params.verbose = 'off';
        params.timewarp_flag = 0;
        params.timewarp = []; % input should be event times
    end

    if isfield(params, 'timewarp_flag')
        if params.timewarp_flag == 1
            if isfield(params, 'alpha')
                [ersp,itc,~,~,~,~,~,tfdata] = newtimef(trial_data,sum(abs(params.tlimits))+1,params.tlimits,...
                    params.sfreq,params.cycles,'winsize',params.winsize,'maxfreq',params.maxfreq,...
                    'baseline',params.baseline,'alpha',params.alpha,'freqs',params.freqs,...
                    'freqscale',params.freqscale,'plotphase',params.plotphase,...
                    'padratio',params.padratio,'mcorrect',params.mcorrect,'plotitc',params.plotitc,'verbose',params.verbose,...
                    'timewarp',params.timewarp);
            else
                [ersp,itc,~,~,~,~,~,tfdata] = newtimef(trial_data,sum(abs(params.tlimits))+1,params.tlimits,...
                    params.sfreq,params.cycles,'winsize',params.winsize,'maxfreq',params.maxfreq,...
                    'baseline',params.baseline,'freqs',params.freqs,...
                    'freqscale',params.freqscale,'plotphase',params.plotphase,...
                    'padratio',params.padratio,'plotitc',params.plotitc,'verbose',params.verbose,...
                    'timewarp',params.timewarp);
            end
        else
            if isfield(params, 'alpha')
                [ersp,itc,~,~,~,~,~,tfdata] = newtimef(trial_data,sum(abs(params.tlimits))+1,params.tlimits,...
                    params.sfreq,params.cycles,'winsize',params.winsize,'maxfreq',params.maxfreq,...
                    'baseline',params.baseline,'alpha',params.alpha,'freqs',params.freqs,...
                    'freqscale',params.freqscale,'plotphase',params.plotphase,...
                    'padratio',params.padratio,'mcorrect',params.mcorrect,'plotitc',params.plotitc,'verbose',params.verbose);
            else
                [ersp,itc,~,~,~,~,~,tfdata] = newtimef(trial_data,sum(abs(params.tlimits))+1,params.tlimits,...
                    params.sfreq,params.cycles,'winsize',params.winsize,'maxfreq',params.maxfreq,...
                    'baseline',params.baseline,'freqs',params.freqs,...
                    'freqscale',params.freqscale,'plotphase',params.plotphase,...
                    'padratio',params.padratio,'plotitc',params.plotitc,'verbose',params.verbose);
            end
        end
    else
        [ersp,itc,~,~,~,~,~,tfdata] = newtimef(trial_data,sum(abs(params.tlimits))+1,params.tlimits,...
            params.sfreq,params.cycles,'winsize',params.winsize,'maxfreq',params.maxfreq,...
            'baseline',params.baseline,'alpha',params.alpha,'freqs',params.freqs,...
            'freqscale',params.freqscale,'plotphase',params.plotphase,...
            'padratio',params.padratio,'mcorrect',params.mcorrect,'plotitc',params.plotitc,'verbose',params.verbose);
    end
    
    disp('Ook')
    
end