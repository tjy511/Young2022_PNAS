% Script to track the bed through convergence of bursts over a day. 
% 
% Input: Range-processed data 'StoreApRES14_RangeProcess.mat', which is the
% time series of range-processed data (i.e. the output of running the
% 'fmcw_range.m' script from the fmcw package given in Stewart (2018)). 
% 
% The script is set up in three parts, of which you will need to run each
% part consecutively (i.e. use the 'Run Section' button, or CMD+Enter): 
% 0. Set up script. Establishes configurations of script and bins bursts
%    into days. If the first set of bursts does not have enough bursts to
%    determine the bed, it is scrapped. 
% 1. Run script to determine range to bed: 
%   A. Align Bulk to remove movement at surface. All bursts (g) are 
%      compared to the first burst. 
%   B. Locate bed window. Manually locates bed of starting day, and then uses
%      this knowledge to track the bed location through time. This is done by
%      identifying a conservative bed "window" that restricts the search
%      range. 
%   C.  Determine range of bed. Within the bed window, there exists a number
%      of clusters that could represent the bed. The number of clusters are 
%      first identified manually during the first day, then clusters are
%      formed through 1-D K-means analysis. The bed is then identified as the
%      mean of the cluster closest to the previous day's bed depth. If the
%      bed is lost, the window is moved up or down, and step 3 is repeated. 
% 2. 
%
% TJ Young <tjy22@cam.ac.uk>
% 31 August 2016

%% 0. Set up script

% Load data
ccc
fileIn = '/Volumes/GoogleDrive/My Drive/Manuscripts/2020Young_StoreBM/repository/StoreApRES14_RangeProcess.mat.';
load(fileIn)

% Set up and define config
global cfg
cfg.bulkAlignRange = [20 40];
cfg.maxOffsetM = 0.5;
cfg.verbose = 0;
cfg.setup = 'qMono';
cfg.doPlotAll = 0;
cfg.doPlotAlignBulk = 0;
cfg.bedSearchRange = [608 620]; % Bed search range (m)
cfg.coarseChunkWidth = 4; % Segments to partition bed search (m)

% Bin bursts into days
[tc(:,1),tc(:,2),tc(:,3),tc(:,4),tc(:,5),tc(:,6)] = datevec(t); % YMDHMS
tciy = accumarray(tc(:,2:3),1); % Observations per day to year matrix
tci = reshape(tciy',[],1); % Reshape to linear vector
tci = tci(tci~=0);

% Manually get rid of first chirp if it is too short
if tci(1) < round(mean(tci))
    t = t(1+tci(1):end);
    rangeCoarse = rangeCoarse(1+tci(1):end,:);
    rangeFine = rangeFine(1+tci(1):end,:);
    specCor = specCor(1+tci(1):end,:);
    specRaw = specRaw(1+tci(1):end,:);
    tci = tci(2:end);
end

dr = mean(diff(rangeCoarse(1,:))); % Bin width (m)

%% 1. Run script by day, to determine range to bed
%clear bedDepth bedDepths bedCheck
wincount = 0;
counter = 1; 
t0 = datenum(tc(1,1),tc(1,2),tc(1,3)); % Starting date
for dn = 1:length(tci);
    % Update indices
    win = [wincount+1 wincount+tci(dn)];
    ind = win(1):win(2);
    wincount = win(2);
    tr(dn) = t0 + dn;
    
    % Assign bursts within day to variables to be used within this script
    rcday = rangeCoarse(ind,:);
    rfday = rangeFine(ind,:);
    scday = specCor(ind,:);
    srday = specRaw(ind,:);
    
    % A. Align Bulk to remove movement at surface
    
    % Set up f
    f.rc = rcday(1,:);
    f.rf = rfday(1,:);
    f.sc = scday(1,:);
    f.sr = srday(1,:);

    for ii = 2:size(rcday,1)
        % Set up f
        g.rc(ii,:) = rcday(ii,:);
        g.rf(ii,:) = rfday(ii,:);
        g.sc(ii,:) = scday(ii,:);
        g.sr(ii,:) = srday(ii,:);
        
        % Run Align Bulk script
        [AB.maxCor(ii),AB.n(ii),g.rf(ii,:),g.sc(ii,:),...
            g.sr(ii,:),g.rfu(ii,:),g.scu(ii,:),g.sru(ii,:)] = ...
            fmcw_alignbulk(f.rc,f.sc,g.rc(ii,:),g.rf(ii,:),g.sc(ii,:),g.sr(ii,:),dr);
    end
    
    % B. Locate bed window
   
    % Obtain starting window (and associated bins) of first bed depth
    if dn == 1
        
        % Manually pick starting location of bed
        %Text = 'Pick window to zoom.' 
        [AC.depthStart,AC.winStart] = pickBedConverge(rcday,scday,srday);
        
        % Obtain indices of bin locations
        AC.winWidth = floor(cfg.coarseChunkWidth/dr);
        [~,AC.winSearchRange(1)] = closest(f.rc,cfg.bedSearchRange(1)); % Bounds of bed search area
        [~,AC.winSearchRange(2)] = closest(f.rc,cfg.bedSearchRange(2)-cfg.coarseChunkWidth); % Bounds of search area
        AC.binStarts = AC.winSearchRange(1):ceil(AC.winWidth/1):AC.winSearchRange(2); % Indices of bin start locations
        AC.depthStarts = cfg.bedSearchRange(1):cfg.coarseChunkWidth:cfg.bedSearchRange(2)-cfg.coarseChunkWidth; % Indices of depth start locations
        
        for ii = 1:numel(AC.binStarts)
            winRange = [AC.binStarts(ii) AC.binStarts(ii)+AC.winWidth]; % Bounds of window search area
            AC.fi(ii,:) = winRange(1):winRange(2); % Index of range bins within search area
        end
        
        % Find window and bin that bed depth is located in
        %tmp = AC.depthStart - AC.depthStarts;
        tmp = median(AC.fi,2) - AC.winStart;
        tmp = abs(tmp);
        %tmp(tmp<0) = NaN;
        %AC.fiLoc(dn) = find(AC.depthStarts == floor(AC.depthStart));
        [~,AC.fiLoc(dn)] = min(tmp);
        AC.bins(dn,:) = AC.fi(AC.fiLoc(dn),:); % Bed bins to use
    % Use the same window as the previous burst
    % If the bed is lost, the bin will move at a later stage. 
    else
        AC.fiLoc(dn) = AC.fiLoc(dn-1); % Window number to use
        AC.bins(dn,:) = AC.fi(AC.fiLoc(dn),:); % Bed bins to use
        
        %% Manually pick bed every n days to "train" script
        bedCheck.depth(1) = AC.depthStart;
        bedCheck.win(1) = AC.winStart;
        bedCheck.t(1) = tr(1);
        if mod(dn-1,28) == 0 % Run every 21 days
            counter = counter+1;
            bedCheck.t(counter) = tr(dn);
            [bedCheck.depth(counter),bedCheck.win(counter)] = pickBedConverge(rcday,scday,srday);
        end
    end
    
    % C. Determine range of bed
    
    % Determine possible bed depth(s)
    bedsday = rcday(:,AC.bins(dn,:)) + rfday(:,AC.bins(dn,:));
    bedsday = reshape(bedsday,[numel(bedsday),1]);
    
    if dn == 1
        % Manually identify number of clusters. 
        % Only need to do this at the start. 
        figure, hold on
        title('Type in the number of clusters that you see in this figure.');
        plot(bedsday, '.')
        clust.num = input('Enter the number of clusters that you see on the figure: ');
        close(gcf);
        
        bedday = median(rcday(:,AC.winStart) + rfday(:,AC.winStart)); % Starting bed depth (only for dn == 1)
        disp(['Obtained bed range for day: 1/132'])
    end
    
    % Find total ranges of possible bed depth clusters
    bedDepths(dn,:) = bedCluster(bedsday, clust.num);
    
    % Identify cluster closest to location of chosen  bed
    if dn == 1
        clust.separation = nanmean(diff(sort(bedDepths(dn,:))));
        [~,clust.pick] = closest(bedDepths(dn,:),bedday);
        bedDepth(dn,:) = bedDepths(dn,clust.pick);
    else
        % Track bed by matching closest cluster with the bed depth of the last burst
        [~,clust.pick] = closest(bedDepths(dn,:),bedDepth(dn-1,:));
        bedDepth(dn,:) = bedDepths(dn,clust.pick);
        bedShift = bedDepth(dn,:)-bedDepth(dn-1,:);
        
        % Check if we need to move bins--the bed has moved out of the bin
        if abs(bedShift) > clust.separation/2 
            try
            if bedShift > 0
                % Move up a window: bed has decreased out of window range
                AC.fiLoc(dn) = AC.fiLoc(dn-1)-1;
                disp('Bed depth moved out of window: Shifting search window right...')
                
            elseif bedShift < 0
                % Move down a window: bed has increased out of window range
                AC.fiLoc(dn) = AC.fiLoc(dn-1)+1;
                disp('Bed depth moved out of window: Shifting search window left...')
            end
            catch
              print('Warning: Bed has jumped over margin of error.')
            end
            AC.bins(dn,:) = AC.fi(AC.fiLoc(dn),:); % New bed bins to use
            
            % Repeat C using new bed bins
            bedsday = rcday(:,AC.bins) + rfday(:,AC.bins);
            bedsday = reshape(bedsday,[numel(bedsday),1]);
            
            bedDepths(dn,:) = bedCluster(bedsday, clust.num); % Find total ranges of possible bed depth clusters
            
            [~,clust.pick] = closest(bedDepths(dn,:),bedDepth(dn-1,:));
            bedDepth(dn,:) = bedDepths(dn,clust.pick); % Track closest cluster
        end
        disp(['Obtained bed range for day: ',num2str(dn),'/',num2str(numel(tci))])
    end
end

%% Plotting fancies
F1 = figure; hold on, box on
set(F1,'position',[100,100,600,400])
plot(tr,bedDepths,'.')
plot(tr,bedDepth,'k')
plot(bedCheck.t,bedCheck.depth,'ko')
dynamicDateTicks

% Manually pick bed curve

% Old bedDepth
F2 = figure; hold on, box on
set(F2,'position',[800,100,600,400])
plot(tr,bedDepths,'.')
plot(tr,bedDepth,'k')
plot(bedCheck.t,bedCheck.depth,'ko')
dynamicDateTicks 
ylim([612 616])

% Note: the dynamicDateTicks script is available on MatLab file exchange: 
% https://uk.mathworks.com/matlabcentral/fileexchange/27075-intelligent-dynamic-date-ticks

[ibx,iby] = ginput; % Press return key to exit

close(F2);

% Convert hand-drawn curve to real data
ibiy = interp1(ibx,iby,tr);

% Pick nearest points from drawn line
for ii = 1:numel(ibiy)
    idx(ii) = knnsearch(bedDepths(ii,:)',ibiy(ii));
    bd(ii) = bedDepths(ii,idx(ii));
end

% New bedDepth
figure, hold on
plot(tr,bedDepths,'.')
plot(tr,bedDepth,'k')
plot(bedCheck.t,bedCheck.depth,'ko')
plot(tr,bd,'k','lineWidth',2)
dynamicDateTicks
ylim([612 616])

%% Export data (comment out if running entire script)
cd('/Users/tjy511/Documents/School/PhD/radar/results/basalmelt/convergence')
save('bed_convergence.mat','tr','bedDepth','bedDepths','bedCheck','bd')
