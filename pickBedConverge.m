function [loc,bin] = pickBedConverge(rangeCoarse,specCor,specRaw,Text)

% Quick function to locate bed range bin.
%
% TJ Young
% 30 August 2016

% Plot figure
h = figure; hold on, box on
sp(1) = subplot(2,1,1); hold on
sp(2) = subplot(2,1,2); hold on
for ii = 1:size(rangeCoarse,1)
    plot(sp(1),rangeCoarse(ii,:),db(specCor(ii,:)));
    plot(sp(2),rangeCoarse(ii,:),angle(specRaw(ii,:)));
end
xlabel('Depth (m)')
ylabel('Vrms (dB)')
linkaxes(sp,'x')
%if nargin == 4
    %title(strcat('Burst date: ',datestr(t)))
%end

% Identify area to zoom in (2x)
title(sp(1),'Click TWICE in upper panel, to set lower/upper bounds for the next window to zoom in to.');
[zm,~] = ginput(2); % Zoom in
xlim([zm(1) zm(end)]) % Zoom in
title(sp(1),'Click TWICE in upper panel, to set lower/upper bounds for the next window to zoom in to.');
[bd,~] = ginput(2); % Select bounds for peak
xlim([bd(1) bd(end)]); % Zoom in

title(sp(1),'Click ONCE in upper panel, as close to peak that will be used for analysis.');

% Pick location of convergence
[loc,~] = ginput(1);

yLimits = get(gca,'yLim');
close(gcf);

[loc,bin] = closest(rangeCoarse(1,:),loc);