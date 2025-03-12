function coef = tdrRegression(data,regpars,plotflag, task_variable_num)
% tdrRegression linear regression on population responses
%
% Inputs:
%  data: population response (sequential or simultaneous)
%  regpars.regressor: cell array of regressors {nrg 1}
%     Each entry is a string, corresponding to the name of one of the fields
%     in data.task_variable.
%  regpars.regressor_normalization:
%     (1) normalize regressors by maximum ('max_abs')
%     (2) use raw regressors ('none')
%  plotflag: summary plot, yes (1) or no (0)
%
% Outputs:
%  coef.name: regressor names (same as regpars.regressor)
%  coef.response: regression coefficients over units and times [nun npt nrg]
%     nun: n of units
%     npt: n of time samples
%     nrg: n of regressors
%  coef.time: time axis
%  coef.dimension: dimension labels
%
% coef = tdrRegression(data,regpars,plotflag)


% Default inputs
if nargin<3
    plotflag = 0;
end

% Initialize
coef.name = regpars.regressor;

% Check if data is sequentially or simultaneously recorded
if isfield(data,'unit') && ~isempty(data.unit)

    
    %--- Sequential recordings ---    
    % Dimensions
    npt = length(data.time);
    nun = length(data.unit);
    nrg = length(regpars.regressor);
    
    % Initialize
    bb = zeros(nun,npt,nrg); % regression coefficients
    
    % Loop over units
    for iun = 1:nun
        
        % Construct the regressor values from the task variables
        ntr = size(data.unit(iun).response,1);
        regmod = ones(ntr,nrg);
        
        % Get the regressors
        for irg = 1:nrg
            
            % Loop over task variables
            isep = [0 strfind(regpars.regressor{irg},'*') length(regpars.regressor{irg})+1];
            for jsep = 1:length(isep)-1
                
                % The task variable
                varname = regpars.regressor{irg}(isep(jsep)+1:isep(jsep+1)-1);
                
                % Check that not constant term
                if ~strcmp(varname,'b0');
                    % Update regressor
                    regmod(:,irg) = ...
                        regmod(:,irg).*data.unit(iun).task_variable.(varname);
                end
            end
        end
        
        % Normalize the regressors
        switch regpars.regressor_normalization
            case 'none'
                regnrm = regmod;
                
            case 'max_abs'
                regnrm = zeros(size(regmod));
                for irg = 1:nrg
                    regnrm(:,irg) = regmod(:,irg) / max(abs(regmod(:,irg)));
                end
        end
        
        % Loop over time points
        bbt = zeros(npt,nrg);
        for ipt = 1:npt
            % Responses to predict
            yy = squeeze(data.unit(iun).response(:,ipt));
            
            % Linear regression
            bbt(ipt,:) = regress(yy,regnrm);
        end
        
        % Keep coefficients
        bb(iun,:,:) = bbt;
    end
    
    % The state space dimensions
    dimension = {data.unit(:).dimension}';
    

else
    
    
    %--- Simultaneous recordings ---
    % Dimensions
    [nun npt ntr] = size(data.response);
    nrg = length(regpars.regressor);
    dummy_num = task_variable_num*5;
    
    % Initialize
%     bb = zeros(nun,npt,nrg); % regression coefficients
    bb = zeros(nun,npt,dummy_num); % regression coefficients
    
    % Construct the regressor values from the task variables
    regmod = ones(ntr,nrg);
    
    % Get the regressors
    for irg = 1:nrg
        
        % Loop over task variables
        isep = [0 strfind(regpars.regressor{irg},'*') length(regpars.regressor{irg})+1];
        for jsep = 1:length(isep)-1
            
            % The task variable
            varname = regpars.regressor{irg}(isep(jsep)+1:isep(jsep+1)-1);
            
            % Check that not constant term
            if ~strcmp(varname,'b0')
                % Update regressor
                regmod(:,irg) = ...
                    regmod(:,irg).*data.task_variable.(varname);
            end
        end
    end
        
    % Normalize the regressors
    switch regpars.regressor_normalization
        case 'none'
            regnrm = regmod;
            
        case 'max_abs'
            regnrm = zeros(size(regmod));
            for irg = 1:nrg
                regnrm(:,irg) = regmod(:,irg) / max(abs(regmod(:,irg)));
            end            
    end

    % 回归自变量中心化
%     regnrm(:,[2,3]) = regnrm(:,[2,3])-mean(regnrm(:,[2,3]),1);
%     regnrm(:,4) = regnrm(:,2).*regnrm(:,3);

    % 按分类变量生成回归design matrix

    regnrm = x2fx(regnrm(:,2:2+task_variable_num-1), 'linear', [1:1+task_variable_num-1]); % 不加入交互项

    % Loop over units
    for iun = 1:nun
        % Loop over time points
%         bbt = zeros(npt,nrg);
        bbt = zeros(npt,dummy_num);
        for ipt = 1:npt
            % Responses to predict
            yy = squeeze(data.response(iun,ipt,:));
            
            % Linear regression            
%             bbt(ipt,:) = regress(yy,regnrm);
%             bbt(ipt,1:2) = regress(yy,regnrm(:,[1,4]));
            tmp = regress(yy,regnrm);
            bbt(ipt,:) = tmp(2:dummy_num+1);

%             b0 = regnrm(:,1);
%             r1 = regnrm(:,2);
%             r2 = regnrm(:,3);
% %             r1r2 = regnrm(:,4);
%             tbl = table(yy,b0,r1,r2);
%             tbl.r1 = categorical(tbl.r1);
%             tbl.r2 = categorical(tbl.r2);
%             lm = fitlm(tbl,'yy~r1+r2');
%             bbt(ipt,:) = lm.Coefficients.Estimate(2:11);
%             bbt(ipt,1:4) = fitlm(tbl,'interactions',{'r1','r2'}, 'ResponseVar','yy',...
%                 'PredictorVars',{'r1','r2'});
        end
        
        % Keep coefficients
        bb(iun,:,:) = bbt;
    end
    
    % The state space dimensions
    dimension = data.dimension;
    
end


% The norm of the raw regression vectors
nrg = dummy_num;
bnorm = zeros(npt,nrg);
for irg = 1:nrg
    for ipt = 1:npt
        bnorm(ipt,irg) = norm(squeeze(bb(:,ipt,irg)));
    end
end

% Plot the norms
plotflag = 0;
if plotflag
    % Plot
    figure; hp=plot(data.time,bnorm);
    
    % Labels
    for irg = 1:nrg
%         ht=text(data.time(end),bnorm(end,irg),['  ' regpars.regressor{irg}],...
%             'interpreter','none');
%         set(ht,'horizontalalignment','left','verticalalignment','middle',...
%             'color',get(hp(irg),'color'));
    end
    set(gca,'xlim',[data.time(1) data.time(end)],'ylim',[0 inf]);
    xlabel('time (s)'); ylabel('regressor norm');
    
end


% Keep what you need
coef.response = bb;
coef.time = data.time;
coef.dimension = dimension;
% regpars.norm = bnorm;


