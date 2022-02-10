function [maxCor,n,grf,gsc,gsr,grfu,gscu,gsru] = fmcw_alignbulk(frc,fsc,grc,grf,gsc,gsr,dr)
% Script for ALIGN BULK
% Bulk Align upper internals (to keep the search window smaller for the
% fine scale correlation)
% Several configurations need to be set before running the script:
% bulkAlignRange, maxOffsetM
%
% Inputs:
% - frc/fsc: coarseRange and specCor for Shot 1
% - grc/grf/gsc/gsr: coarseRange, fineRange, specCor, specRaw for Shot 2
% - dr: Bin width
%
% Outputs:
% - maxCor: Best correlation found
% - n: Lag (from f) at which best correlation was found
% - gsc/gsr: Shifted specCor and specRaw for Shot 2
% - grfu/gscu/gsru: Unshifted rangeFine, specCor, specRaw for Shot 2
%
% Craig Stewart, TJ Young
% 05 August 2016

%% Parameters

global cfg

% Allign internal layers to account for changes in cable length surface
% accumulation and firn compaction using xcor. Note this only offsets to
% the closest integer range bin.
Disp(['Co-registering profiles using amplitude cross-correlation'])
Disp(['> depth range: ' mat2str(cfg.bulkAlignRange)])
Disp(['> max offset : ' mat2str(cfg.maxOffsetM)])

% Find depth bins to be used (f)
fi = find((frc>=min(cfg.bulkAlignRange) & frc<max(cfg.bulkAlignRange)));
maxlag = ceil(cfg.maxOffsetM/dr); % max bin lags

% Cross-correlate f and g
if strcmp(cfg.setup,'MIMO') == 1
    % This is only used if I am using this script to shift chirps in cables on MIMO radar
    [~,ampCor,~,lags] = fmcw_xcorr_array(fsc,gsc,fi,maxlag);
    % Assume that cable shift is positive, i.e. cable(g) is longer/equal to cable(f)
    ampCor = ampCor(length(ampCor)/2:end);
    lags = lags(length(lags)/2:end);
else
    [~,ampCor,~,lags] = fmcw_xcorr(fsc,gsc,fi,maxlag);
end

% Obtain best amplitude correlation
[maxCor,ii] = max(ampCor); % get index (mci) of best amplitude correlation
n = lags(ii); % n is the number of steps g should be shifed right to match f
Disp(['correlation =  ' num2str(maxCor) ])
if maxCor < 0.8
    Disp(['Warning: poor correlation in bulk internal allignment - check files'])
end

% Create copies of g for plotting prior to range shift
grfu = grf;
gsru = gsr; % keep a copy of g.sr for plotting
gscu = gsc; % keep a copy of g.sc for plotting

% Apply the offset to shot 2 to make this match shot 1
if n==0
    Disp('Internals match - no offset required')
    Disp(' ')
    
else
    Disp(['Shifting profile 2, ' int2str(n) ' steps left to align internals. (' num2str(n*dr) 'm)'])
    Disp(' ')
    gsr = circshift(gsr,[0 -n]); % lagg offset
    gsc = circshift(gsc,[0 -n]); % lagg offset
    grf = circshift(grf,[0 -n]); % lagg offset
    
    if cfg.doPlotAlignBulk || cfg.doPlotAll
        % plot before and after offset
        figure
        plot(frc,dB(abs(fsc)),'r');
        hold on
        plot(grc,dB(abs(gscu)),'b');
        plot(grc,dB(abs(gsc)),'c');
        legend('shot1','shot2','shot2 shifted')
        title(['Profile bulk co-registration, depth range: ' mat2str(cfg.bulkAlignRange) ' m'])
        set(gca,'xlim',cfg.bulkAlignRange)
        set(gcf,'pos',[232 554 560 420])
    end
end


