function exportReadmeFigures(outputDir)
%EXPORTREADMEFIGURES Run the MATLAB examples and save the three README PNGs.
%   Add MATLAB scripts to the path, then call exportReadmeFigures. Pass
%   an output directory to save somewhere other than the repository's figures folder.

if nargin < 1
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    outputDir = fullfile(repoRoot,'figures');
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

examples = {'mergingExample_the','mergingExample_syn','TCvsNC'};
filenames = {'merging-covariance.png','merging-time-series.png','tc-vs-nc.png'};
for k = 1:numel(examples)
    fig = runExample(examples{k});
    closeFigure = onCleanup(@() delete(fig));
    fig.Visible = 'off';
    styleFigure(fig);
    exportgraphics(fig,fullfile(outputDir,filenames{k}), ...
        'Resolution',150,'Padding',12,'BackgroundColor','white');
    clear closeFigure
end
end

function fig = runExample(example)
% The existing example scripts clear their workspace.
run(which(example));
fig = gcf;
end

function styleFigure(fig)
theme(fig,'light');
fig.Units = 'pixels';
fig.Position = [100 100 1040 375];
fig.Color = 'white';
axesList = findobj(fig,'Type','axes');
[~, order] = sort(arrayfun(@(ax) ax.Position(1),axesList));
ax = axesList(order);

if numel(ax) == 2
    positions = [0.09 0.16 0.385 0.66; 0.59 0.16 0.385 0.66];
    delete(findobj(fig,'Type','legend'));
    labels = {'Mean squared error','Squared correlation (R^2)'};
    for j = 1:2
        ylabel(ax(j),labels{j});
        title(ax(j),'');
        ax(j).XTickLabel = {'WA','SNR-opt','maxR','EW'};
        ax(j).YGrid = 'on';
    end
    bars = flipud(findobj(ax(1),'Type','bar'));
    lgd = legend(ax(1),bars,{'Mean individual product','True parameters','NC estimates'});
else
    positions = [0.079 0.18 0.244 0.59; 0.408 0.18 0.244 0.59; 0.737 0.18 0.244 0.59];
    titles = {'Variance difference','Variance MSE','Fit failures'};
    for j = 1:3
        title(ax(j),titles{j});
        ax(j).XLim = [0 1];
        ax(j).XTick = 0:0.25:1;
        ax(j).XGrid = 'off';
    end
    ax(1).YLimMode = 'auto';
    set(findobj(ax(1),'Type','line'),'Color',[0.25 0.28 0.33]);
    set(findobj(ax(1),'Type','patch'),'FaceColor',[0.70 0.74 0.79]);
    ylabel(ax(2),'Error variance MSE');
    ax(2).YLimMode = 'auto';
    ax(2).YTickMode = 'auto';
    ax(3).YLim = [0 100];
    delete(findobj(fig,'Type','legend'));
    lines = flipud(findobj(ax(2),'Type','line'));
    lgd = legend(ax(2),lines,{'TC','NC'});
end

lgd.FontSize = 13;
lgd.Box = 'off';
lgd.Orientation = 'horizontal';
lgd.Location = 'none';
lgd.Units = 'normalized';
lgd.Position(1:2) = [(1-lgd.Position(3))/2, 0.87];

for j = 1:numel(ax)
    axis(ax(j),'normal');
    ax(j).Units = 'normalized';
    ax(j).Position = positions(j,:);
    ax(j).FontSize = 14;
    ax(j).Box = 'off';
    ax(j).TickDir = 'out';
    ax(j).GridAlpha = 0.15;
end
end
