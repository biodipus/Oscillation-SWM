function [coef_fix,coef_nrm] = tdrVectorTimeAverage(coef_time,vecpars,plotflag)
% tdrVectorTimeAverage Define fixed regression vectors by time averaging coefficients
%
% Inputs:
%  coef_time: regression coefficients 
%     .name: regressor names [nvc 1]
%     .response: coefficients [nun npt nvc]
%     .time: time axis [1 npt]
%  vecpars: times to average. Each field of vecpars contains the time
%     window(s) to average for a given coefficient. The fieldname
%     corresponds to the coefficient name.
%     e.g. vecpars.coefname(i).time_win = [tmin tmax] averages the 
%          coefficient 'coefname' over the time window [tmin tmax]. 
%     Can specify many time windows (i=1,2,...).
%     If the time window is empty, uses the time of maximum norm.
%  plotflag: summary plot yes (1) or no (0)
%
% Outputs:
%  coef_fix: fixed (non time-varying) regression vectors
%     .name: regression vector name {nvc 1}
%     .response: regression vector coefficients [ndm 1 nvc]
%     .dimension: dimension labels {ndm 1}
%  coef_nrm(i): time-varying norm of the regression coefficients for vector i
%     .bnorm: vector norm [npt 1]
%     .jpt: time samples that were averaged
%     .irg: index of coefficient
%     .name: name of coefficient
%
% [coef_fix,coef_nrm] = tdrVectorTimeAverage(coef_time,vecpars,plotflag)


if nargin < 3 || isempty(plotflag)
    plotflag = 0;
end

% The vectors
vecname = fieldnames(vecpars);

% Initialize
coef_fix = coef_time;
coef_fix.name = {};

% The number of vectors
nve = 0;
nvc = length(vecname);
for irg = 1:nvc
    nve = nve + length(vecpars.(vecname{irg}));
end

% Dimensions
[nun npt ncf] = size(coef_time.response);

% The norm of the vectors
bnorm = zeros(npt,ncf);
for icf = 1:ncf
    for ipt = 1:npt
        bnorm(ipt,icf) = norm(squeeze(coef_time.response(:,ipt,icf)));
    end
end

% Initialize
response = zeros([nun 1 nve]);
clear plotinfo;

% Loop over vectors
ive = 1;
nvc = length(vecname);
for irg = 1:nvc
    
    % The matching vector in the input
    imatch = find(strcmp(coef_time.name,vecname{irg}));
    
    % Loop over temporal windows
    nwn = length(vecpars.(vecname{irg}));
    for iwn = 1:nwn
        
        if isempty(vecpars.(vecname{irg})(iwn).time_win)
            % Find maximum of the norm
            [~,jpt] = nanmax(bnorm(:,imatch),1);
        else
            % Time points to average
            jpt = find(...
                coef_time.time>=vecpars.(vecname{irg})(iwn).time_win(1) & ...
                coef_time.time<=vecpars.(vecname{irg})(iwn).time_win(end));
            
            % Find closest if no match
            if ~any(jpt)
                [~,jpt] = min(abs(coef_time.time - mean(vecpars.(vecname{irg})(iwn).time_win)));
            end
        end
        
        % Average vectors over time
        response(:,1,ive) = mean(coef_time.response(:,jpt,imatch),2);
        
        % Vector name
        if nwn==1
            coef_fix.name{ive,1} = vecname{irg};
        elseif nwn>1
            coef_fix.name{ive,1} = [vecname{irg} '_' num2str(iwn)];
        end
        
        % Info for plots
        coef_nrm(ive).bnorm = bnorm(:,imatch);
        coef_nrm(ive).jpt = jpt;
        coef_nrm(ive).irg = irg;
        coef_nrm(ive).name = coef_fix.name{ive};
        
        % Update
        ive = ive + 1;
        
    end
end

% Keep what you need
coef_fix.response = response;
coef_fix.time = [];

% PLOT a summary with the norm of the vectors and the times when they were
% averaged
if plotflag
    % Line colors
    lc = linecolors(nvc,'jet');
    
    figure; hold on;
    
    % Loop over vectors
    for ive = 1:nve
        % Times that were averaged
        jpt = coef_nrm(ive).jpt;
        
        % The raw vectors
        hp1=plot(coef_time.time,coef_nrm(ive).bnorm,'-');
        set(hp1,'color',lc(coef_nrm(ive).irg,:));
        
        % The values that were averaged
        hp2=plot(coef_time.time(jpt),coef_nrm(ive).bnorm(jpt),'o');
        set(hp2,'color',lc(coef_nrm(ive).irg,:));
        
        % Axis labels
        xlabel('time (s)'); ylabel('coefficient norm');
        
        % The vector name
        ht=text(coef_time.time(jpt(1)),coef_nrm(ive).bnorm(jpt(1)),coef_fix.name{ive});
        set(ht,'horizontalalignment','left','verticalalignment','top','interpreter','none',...
            'color',lc(coef_nrm(ive).irg,:));
    end
end



