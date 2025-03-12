function [coef_ort,lowUN_lowTA] = tdrVectorOrthogonalize(coef_fix,ortpars)
% tdrVectorOrthogonalize orthogonalize regression vectors
%
% Inputs:
%  coef_fix: regression vectors
%     .name: vector names {nvc 1}
%     .response: vector coefficients [ndm 1 nvc]
%     .dimension: dimension labels {ndm 1}
%  ortpars.name: vectors to orthogonalize {nrg 1}
% 
% Outputs:
%  coef_ort: orthogonalized regression vectors
%     .name: vector names {nrg 1}
%     .response: vector coefficients [ndm 1 nrg]
%     .dimension: dimension labels {ndm 1}
%  lowUN_lowTA: projection matrix from task-related subspace (subspace
%     basis) into original state space (original basis)
%  
% [coef_ort,lowUN_lowTA] = tdrVectorOrthogonalize(coef_fix,ortpars)

% Dimensions
nrg = length(ortpars.name);
nun = size(coef_fix.response,1);

% Initialize
coef_ort = coef_fix;
coef_ort.name = ortpars.name;
coef_ort.response = zeros(nun,1,nrg);

raw = zeros(nun,nrg);
for irg = 1:nrg
    % Find vector
    jmatch = strcmp(coef_fix.name,ortpars.name{irg});
    
    % Keep vector
    raw(:,irg) = coef_fix.response(:,1,jmatch);
end

% Orthogonalize
[qq,rr] = myqr(raw);
ort = qq(:,1:nrg);

% Keep what you need
coef_ort.response(:,1,:) = ort;
lowUN_lowTA = ort;
