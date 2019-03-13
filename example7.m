% Example script on how to setup and launch the data explore function
clear; clc;

% First generate the main figure to serve as the ui window.
mainfig = figure;
mainfig.OuterPosition(1:2) = [30, 50];
mainfig.OuterPosition(3:4) = mainfig.OuterPosition(3:4).*[1.8 2];
x = 0:0.1:10;

cnt = gobjects(0);
for plt = 1:8
    if plt == 1
        ax = subplot(4,2,[plt 2]);
    elseif plt == 3
        ax = subplot(4,2,[plt 7]);
    elseif plt == 2 || plt == 5 || plt == 7
        continue
    else
        ax = subplot(4,2,plt);
    end

    [~,cnt(length(cnt)+1)] = contourf(160*membrane(1,100));
    cb = colorbar;
    if plt == 3
        cb.Location ='northoutside';
    end
%     cb.AxisLocationMode = 'manual';

    title(['Plot ' num2str(plt)])
    xlabel('x')
    ylabel('y')

end

%% Launch exploreData
% exploreResults has an option for display boxes to show information on the
% selected point. It can get this information from the line selected and/or
% the user can provide the information in a cell array.
% uispec = {{ 'Selected Index', 1:length(x), 'y2', y2/max(y2) },...
%           { 'y2', y2/max(y2) },...
%           { 'y2', y2/max(y2) }};

% Push button call back
% This must be a cell array where the first entry is an anonymous function
% and the following are arguments for that function. See the documentation
% on this argument for details on the requirements of this function.
cutoff = 1; % just an argument for the example function
pbtn_callback = {@ userCallback, cnt(1), [], []};

% Set two aditional options
usefigdat = true; % use data from figure for display boxes
linkselect = true; % if there are multiple plots select the same point on all of them (requires each line has the same number of points)

% Finally launch the ui figure
exploreResults( mainfig, pbtn_callback, cnt, [], usefigdat, linkselect );

%% --- Example User Function --- %%
function [ newfig ] = userCallback(src, event, slct, ui, ln, extrapnt, newfig)
% An example of a function to assign to the ui push button. See
% documentatoin for details on the reserved input arguments for src, event,
% slct and ui.

% I recomend placing a breakpoint in here and running the script to take a
% look at the structure of the arguments passed to this function

ind = slct(1).ind; % for this case all indices should be the same

% --- Do something for fun --- %
if ~isempty(extrapnt) && isvalid(extrapnt)
    delete(extrapnt);
end
% Place a marker on the selected axes
extrapnt = line(ln.XData(ind), ln.YData(ind), 'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 50);

% If this point on the third plot is below the cutoff
if ~isempty(newfig) && isvalid(newfig)
    close(newfig);
end
newfig = figure;
surf(160*membrane(1,100), 'EdgeColor', 'none');
axis tight

% --------------------------- %

% We can optionally update the arguments for this callback like so. The 
% user data should be organized in a cell array as with the initial call to
% exploreResults. The first three arguments to the callback are set by 
% exploreResults and normally don't need to be editted here.

user_args = {@ userCallback, ln, extrapnt, newfig};
src.Callback{end} = user_args;

end
