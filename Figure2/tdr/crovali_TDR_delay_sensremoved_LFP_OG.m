%%

%% Paths
clear all;
close all;

% monkey = 'ocean';
monkey = 'grootsw';

if strcmp(monkey, 'ocean')
    alldays = {'0709','0710','0713','0714','0715','0720','0721','0722','0723','0726','0727','0729','0730'};
end
if strcmp(monkey, 'groot')
    alldays = {'03-11','03-14','03-15','03-16','03-19','03-21','03-22','03-23','03-25','03-26','03-27'};
end
if strcmp(monkey, 'grootsw')
    alldays = {'04-15','05-07','05-09','05-13','05-14','05-15','05-16','05-20','05-21','05-22'};
end
% TDR path
tdrDir = 'C:\Wen\OneDrive\Project\Osc_control\Code\';
addpath(fullfile(tdrDir,'tdr'));
addpath(fullfile(tdrDir,'tdr','tools'));
for dayid = 1:length(alldays)
    clearvars -except alldays dayid exp_ratio enc_weight_all monkey
    %% Parameters
    days = alldays{dayid};
    ruletype = ['all']; 
    plotflag = 0;
    removesens = 0;
    sensdim = 2;
    osc_band = {'80_120','4_8'};
    bandid = 1;
    cv = 0;

    if strcmp(monkey, 'ocean')
        time = [50:55]; 

        datafile = ['C:\Wen\Project\Osc_control\Data\Data_ocean\LFP\power_remove_erp\', ...
            'frontal_tf200_', osc_band{bandid}, '_denoise_morlet_', days(1:end), '.mat'];

        infofile = ['C:\Wen\Project\Osc_control\Data\Data_ocean\LFP\power_remove_erp\', ...
            'frontal_tf200_', osc_band{bandid}, '_denoise_morlet_', days(1:end), ];

        if sensdim == 5
            savefile = ['C:\Wen\Project\Osc_control\Data\Result_ocean\LFP\tdr_frontal_L2_croval\tdr_mem_delay1_5dim\', days(1:end),'_weight_mem_delay_', osc_band{bandid}, '.mat'];
        else
            savefile = ['C:\Wen\Project\Osc_control\Data\Result_ocean\LFP\tdr_frontal_L2_croval\tdr_mem_delay1\', days(1:end),'_weight_mem_delay_', osc_band{bandid}, '.mat'];
        end

        if strcmp(ruletype,'repeat') || strcmp(ruletype,'mirror')
            time = [70:75];
            savefile = ['C:\Wen\Project\Osc_control\Data\Result_ocean\LFP\tdr_frontal_L2_croval\tdr_mem_delay2_5dim\',ruletype, days(1:end),'_',osc_band{bandid},...
                '_weight_mem_delay.mat'];
        end
    end

    if strcmp(monkey, 'groot')
        time = [50:65];

        datafile = ['C:\Wen\Project\Osc_control\Data\Data_groot\LFP\power_remove_erp\', ...
            'frontal_L2_tf200_', osc_band{bandid}, '_morlet_', days(1:end), '.mat'];

        infofile = ['C:\Wen\Project\Osc_control\Data\Data_groot\LFP\power_remove_erp\', ...
            'frontal_L2_tf200_', osc_band{bandid}, '_morlet_', days(1:end), ];

        if sensdim == 5
            savefile = ['C:\Wen\Project\Osc_control\Data\Result_groot\LFP\tdr_frontal_L2_croval\tdr_delay1_erp_remove\', days(1:end),'_weight_mem_delay(750ms)', osc_band{bandid}, '.mat'];
        else
            savefile = ['C:\Wen\Project\Osc_control\Data\Result_groot\LFP\tdr_frontal_L2_croval\tdr_delay1_erp_remove\', days(1:end),'_weight_mem_delay(750ms)', osc_band{bandid}, '.mat'];
        end
    end

    %%%% grootsw
    if strcmp(monkey, 'grootsw')
        time = [50:55]; 
        datafile = ['C:\Wen\Project\Osc_control\Data\Data_grootsw\LFP_sw\power_remove_erp\', ...
            'frontal_L2_tf200_', osc_band{bandid}, '_morlet_', days(1:end), '.mat'];

        infofile = ['C:\Wen\Project\Osc_control\Data\Data_grootsw\LFP_sw\power_remove_erp\', ...
            'frontal_L2_tf200_', osc_band{bandid}, '_morlet_', days(1:end), ];

        if sensdim == 5
            savefile = ['C:\Wen\Project\Osc_control\Data\Result_grootsw\LFP_sw\tdr_frontal_L2_croval\tdr_mem_delay1_5dim\', days(1:end),'_weight_mem_delay_', osc_band{bandid}, '.mat'];
        else
            savefile = ['C:\Wen\Project\Osc_control\Data\Result_grootsw\LFP_sw\tdr_frontal_L2_croval\tdr_mem_delay1\', days(1:end),'_weight_mem_delay_', osc_band{bandid}, '.mat'];
        end
        if strcmp(ruletype,'repeat') || strcmp(ruletype,'mirror')
            time = [70:75];
            savefile = ['C:\Wen\Project\Osc_control\Data\Result_groot\LFP_sw\tdr_frontal_L2_croval\tdr_mem_delay2_5dim\',ruletype, days(1:end),'_',osc_band{bandid},...
                '_weight_mem_delay.mat'];
        end
    end
    
    if cv == 1
        savefile = [savefile(1:end-4), '_LO.mat'];
    end

    %% organize Simultanously data for memory
    ch_info_all = readtable(infofile,'Sheet','ch_info_all', 'VariableNamingRule', 'preserve');
    info = readtable(infofile,'Sheet','info');
    load(datafile, 'fr_all', 'mkts_merge','ch_info_all') % neuron*trial*time
    
    %%%% variables to be save
    project_r1_2D_norm_all = nan(size(fr_all,1), sensdim, size(fr_all,2));
    project_r2_2D_norm_all = nan(size(fr_all,1), sensdim, size(fr_all,2));
    r1_resp_pr1_2D_alltime_all = nan(sensdim, size(fr_all,2), size(fr_all,3));
    r2_resp_pr2_2D_alltime_all = nan(sensdim, size(fr_all,2), size(fr_all,3));
    r1_weight_all = nan(size(fr_all,1), size(fr_all,2));
    r2_weight_all = nan(size(fr_all,1), size(fr_all,2));
    
    if cv == 1
        repeat = size(info,1)
    else
        repeat=1
    end

    for rep = 1:repeat
        fr = permute(fr_all, [1,3,2]); % 
        fr = fr+0.00001; 

        % Normalization parameters
        dd.response = fr;
        dd.time = [1:size(fr,2)];
        [meanT,stdT] = tdrMeanAndStd(dd);
        nrmlpars = [];
        nrmlpars.ravg = meanT;
        nrmlpars.rstd = stdT;
        dd = tdrNormalize(dd,nrmlpars);
        fr = dd.response;
        %%% select rule
        if strcmp(ruletype,'repeat')
            idx = info.Rule == 1;
            info = info(idx,:);
            fr = fr(:,:,idx);
        elseif strcmp(ruletype,'mirror')
            idx = info.Rule == 2;
            info = info(idx,:);
            fr = fr(:,:,idx);
        else
            idx = info.Rule > 0;
            info = info(idx,:);
            fr = fr(:,:,idx);
        end
        if cv==1
            [trainIdx, ~] = setdiff([1:size(info,1)], rep);
            testIdx = rep;
            infotest = info(testIdx,:);
            infotrain = info(trainIdx,:);
            frtest = fr(:,:,testIdx);
%             frtest = reshape(frtest, [size(frtest, 1), size(frtest, 2), 1]);
            frtrain = fr(:,:,trainIdx);
        else
            infotest = info;
            infotrain = info;
            frtest = fr;
            frtrain = fr;
        end

        task_index.rank1 = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5 6 6 6 6 6]';
        task_index.rank2 = [2 3 4 5 6 1 3 4 5 6 1 2 4 5 6 1 2 3 5 6 1 2 3 4 6 1 2 3 4 5]';
        task_variable_num = 2;

        dataT.task_variable.correct = infotrain.Correct;
        dataT.task_variable.rank1 = infotrain.Target_1;
        dataT.task_variable.rank2 = infotrain.Target_2;

        dataT.response = mean(frtrain(:,time,:),2); %T1

        dataT.time = time;
        for i = 1:size(frtrain,1)
            dataT.dimension{i,1} = ['unit_', num2str(i)];
        end

        %% remove sensory subspace
        % Averaging parameters
        avgpars = [];
        avgpars.trial = [];
        avgpars.time = [];

        % Normalization parameters
        nrmlpars = [];
        nrmlpars.ravg = meanT;
        nrmlpars.rstd = stdT;                                       

        % Normalize
        dataT_nrml = dataT;

        if removesens == 1 % remove sensory subspace
            %%%% load sens subspace
            tmp = load([sensfile,days(1:end),'_weight_enc250.mat']);
            project_enc_2D_norm = tmp.project_enc_2D_norm;
            tmp = squeeze(dataT_nrml.response)';
            tmp2 = tmp*project_enc_2D_norm*project_enc_2D_norm';
            tmp3 = squeeze(dataT_nrml.response(:,1,:)) - tmp2';
            dataT_nrml.response(:,1,:) = tmp3;
        end

        %% data for plot the whole trajectory in projection subspace
        dataRaw.response = frtest(:,:,:);
        dataRaw.time = [1:1:size(frtest,2)];

        % Normalization parameters
        nrmlpars = [];
        nrmlpars.ravg = meanT;
        nrmlpars.rstd = stdT;

        dataRaw_nrml = dataRaw;
        %% Linear regression
        regpars = [];

        regpars.regressor = {'b0';'rank1';'rank2';'rank1*rank2'}; % memory
        regpars.regressor_normalization = 'none'; % max_abs

        % Linear regression
        coef_fulUN = tdrRegression_cos(dataT_nrml, regpars, plotflag, task_variable_num, 10);
        %% define memory factor axis
        project_r1 = squeeze(mean(coef_fulUN.response(:,:,1:5),2));
        project_r2 = squeeze(mean(coef_fulUN.response(:,:,6:10),2));

        %%% define r1 subspace
        [pc2un,scores,variances] = pca(project_r1,'Economy',false); % pca on response in beta space
        project_r1_2D = project_r1 * pc2un(:,1:sensdim); % select first 2 dim of pca space

        for i = 1:size(project_r1_2D,2) % loop on dims
            project_r1_2D_norm(:,i) = project_r1_2D(:,i)./norm(project_r1_2D(:,i)); % the final projection matrix
        end

        %%% define r2 subspace
        [pc2un,scores,variances] = pca(project_r2,'Economy',false); % pca on response in beta space
        project_r2_2D = project_r2 * pc2un(:,1:sensdim); % select first 2 dim of pca space

        for i = 1:size(project_r2_2D,2)
            project_r2_2D_norm(:,i) = project_r2_2D(:,i)./norm(project_r2_2D(:,i));
        end

        kk=[project_r1_2D_norm,project_r2_2D_norm];
        [qq,rr] = qr(kk);
        project_r1_2D_norm = qq(:,1:sensdim);
        project_r2_2D_norm = qq(:,sensdim+1:sensdim*2);

        %%%% projection whole trial
        time_num = size(dataRaw_nrml.response, 2);
        trial_num = size(dataRaw_nrml.response, 3);

        %%%% R1 subspace
        r1_resp_pr1_2D_alltime = zeros(sensdim, trial_num, time_num);
        for i = 1:time_num
            r1_resp_pr1_2D_alltime(:,:,i) = project_r1_2D_norm' * squeeze(dataRaw_nrml.response(:,i,:));
        end
        %%%% R2 subspace
        r2_resp_pr2_2D_alltime = zeros(sensdim, trial_num, time_num);
        for i = 1:time_num
            r2_resp_pr2_2D_alltime(:,:,i) = project_r2_2D_norm' * squeeze(dataRaw_nrml.response(:,i,:));
        end

        %%%% explained variance
        var_sub = nan(time_num,2);
        var_ori = nan(time_num,1);
        for i = 1:time_num
            tmp = r1_resp_pr1_2D_alltime(:,:,i)';
            var_sub(i,1) = sum(var(tmp));
            tmp = r2_resp_pr2_2D_alltime(:,:,i)';
            var_sub(i,2) = sum(var(tmp));
            tmp = squeeze(dataRaw_nrml.response(:,i,:))';
            var_ori(i) = sum(var(tmp));
        end
        exp_ratio{dayid} = var_sub./var_ori;

        %%% feature weight
        for i =1:size(project_r1_2D_norm,1)
            r1_weight(i,1) = norm(project_r1_2D_norm(i,:));
            r2_weight(i,1) = norm(project_r2_2D_norm(i,:));
        end
        if cv == 1
            project_r1_2D_norm_all(:,:,rep) = project_r1_2D_norm;
            project_r2_2D_norm_all(:,:,rep) = project_r2_2D_norm;
            r1_resp_pr1_2D_alltime_all(:,rep,:) = r1_resp_pr1_2D_alltime;
            r2_resp_pr2_2D_alltime_all(:,rep,:) = r2_resp_pr2_2D_alltime;
            r1_weight_all(:,rep) = r1_weight;
            r2_weight_all(:,rep) = r2_weight;
        end
    end
    if cv == 1
        save(savefile, 'r1_weight', 'r2_weight','project_r1_2D_norm_all','r1_resp_pr1_2D_alltime_all',...
            'project_r2_2D_norm_all', 'r2_resp_pr2_2D_alltime_all', 'ch_info_all','var_sub','var_ori');
    else
        save(savefile, 'r1_weight', 'r2_weight','project_r1_2D_norm','r1_resp_pr1_2D_alltime',...
            'project_r2_2D_norm', 'r2_resp_pr2_2D_alltime', 'ch_info_all','var_sub','var_ori');
    end
    save(savefile)

end


