function [hh,tt,jj] = boxresample(rr,R,S,Fs,t0)
% boxresample resamples a signal with a sliding box
%
% Inputs are:
%  rr: raw signal [ntr ntt]
%  R: size of the box (in samples)
%  S: step size (in samples)
%  Fs: sampling frequency
%  t0: time of first sample in raw signal
%
% Outputs:
%  hh: resampled signal
%  tt: center time of sliding box
%  jj: sample indeces in each each window
%
% [hh,tt,jj] = boxresample(rr,R,S,Fs,t0)


if nargin<3 || isempty(S)
    S = R;
end
if nargin<4
    Fs = [];
end
if nargin<5
    t0 = 0;
end

% Dimensions
[nr,nt] = size(rr);

% Start index for the sliding box
ib = 1:S:nt-R+1;
nb = length(ib);

% Initialize
hh = zeros(nr,nb);
jj = cell(1,nb);

% Loop over trials
% for ir = 1:nr
%     for jb = 1:nb
%         hh(ir,jb) = nanmean(rr(ir,ib(jb):ib(jb)+R-1),2);
%     end
% end
for jb = 1:nb
    hh(:,jb) = nanmean(rr(:,ib(jb):ib(jb)+R-1),2);
    jj{1,jb} = ib(jb):ib(jb)+R-1;
end

% Center of sliding box
if isempty(Fs)
    tt = [];
else
    tR = R/Fs;
    tt = (ib-1)/Fs + tR/2 + t0;
end



