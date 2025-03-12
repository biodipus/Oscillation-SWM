function [data_subSUB,data_subORG,varnrm] = tdrSubspaceProjection(data_fulORG,subpars,plotflag)
% tdrSubspaceProjection Project population response into subspace
%
% Inputs:
%  data_fulORG: population response in orignal space (original basis)
%  subpars.subSUB_fulORG: projection matrix from original space (original
%     basis) into subspace (subspace basis)
%  subpars.dimension: labels of subspace dimensions {ndm 1}
% 
%  Outputs:
%   data_subSUB: population response in subspace (subspace basis)
%   data_subORG: population response in subspace (original basis)
%   varnrm: total variance along subspace axes (normalized)
%
% [data_subSUB,data_subORG,varSUB,varORG] = tdrSubspaceProjection(data_fulORG,subpars,plotflag)

% Basis transformation
subSUB_fulORG = subpars.subSUB_fulORG;

% Inputs
if nargin<3 || isempty(plotflag)
    plotflag = 0;
end

% Dimensions
[~,npt,ntr] = size(data_fulORG.response);

% Fold time and trials/conditions
Resp_fulORG = data_fulORG.response(:,:);

% Initialize
data_subSUB = data_fulORG;
data_subORG = data_fulORG;

% Keep ntrial (if simultaneous recordings)
if isfield(data_fulORG,'n_trial') && ~isempty(data_fulORG.n_trial)
    if size(data_fulORG.n_trial,1)==1
        data_subSUB.n_trial = data_fulORG.n_trial;
        data_subORG.n_trial = data_fulORG.n_trial;
    else
        % For sequential recordings n_trial is only meaningful at the level of
        % individual neurons, not along arbitrary state space dimensions.
        data_subSUB.n_trial = [];
        data_subORG.n_trial = [];
    end
end

% Project data into subspace
Resp_subSUB = subSUB_fulORG * Resp_fulORG;

% Represent subspace data in original basis
subORG_subSUB = subSUB_fulORG';
Resp_subORG = subORG_subSUB * Resp_subSUB;

% Unfold time and trials/conditions
data_subSUB.response = reshape(Resp_subSUB,[size(Resp_subSUB,1) npt ntr]);
data_subORG.response = reshape(Resp_subORG,[size(Resp_subORG,1) npt ntr]);

% Variance along subspace axes
varSUB = nanvar(Resp_subSUB,2)';

% Total variance
varORG = nanvar(Resp_fulORG,2)';

% Normalize variance
varnrm = varSUB / sum(varORG);

% The dimension names
data_subSUB.dimension = subpars.dimension;
data_subORG.dimension = data_fulORG.dimension;

% Plot
if plotflag
    
    figure; ha = [];
    
    ha(1) = subplot(1,2,1); 
    plot(100* cumsum(varSUB / sum(varORG)),'o-');
    xlabel('subspace axis'); ylabel('cumulative variance [%tot]');
    
    ha(2) = subplot(1,2,2);
    plot(log10(100* varSUB / sum(varORG)),'o-');
    xlabel('subspace axis'); ylabel('variance explained (log10[%tot])');
        
    set(ha,'plotbox',[1 1 1],'xlim',[-inf inf]);
    
end


