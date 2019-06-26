% Example script on how to setup and launch the data explore function

% First generate the main figure to serve as the ui window.
mainfig = figure;
mainfig.OuterPosition(1:2) = [30, 50];
mainfig.OuterPosition(3:4) = mainfig.OuterPosition(3:4)*1.8;
x = 0:0.1:10;

for plt = 2:4
    if plt == 3
        ax = subplot(2,2,[1 plt]);
    else
        ax = subplot(2,2,plt);
    end

    if plt == 2
        plot(x, cos(x), 'o');
    elseif plt == 3
        plot(x, sin(x), 'o');
    elseif plt == 4
        plot(x, tan(x), 'o');
    end

    hold on
    y2 = (x).^2;
    plot(x, y2/max(y2), 'x');
    hold off

    title(['Plot ' num2str(plt)])
    xlabel('x')
    ylabel('y')

end

%% Launch exploreResults
% Push button call back
% This must be a cell array where the first entry is an anonymous function
% and the following are arguments for that function. See the documentation
% on this argument for details on the requirements of this function.
cutoff = 1; % just an argument for the example function
pbtn_callback = {'User Function', @ userCallback, lines2(1), [], cutoff, []};

% Finally launch the ui figure
exploreResults( mainfig, pbtn_callback, 'DataFromAxes', true, 'SelectionLinkAxes', true );

%% --- Example User Function --- %%
function [ newfig ] = userCallback(src, event, ui, slct, ln, extrapnt, cutoff, newfig)
% An example of a function to assign to the ui push button. See
% documentatoin for details on the reserved input arguments for src, event,
% slct and ui.

% I recomend placing a breakpoint in here and running the script to take a
% look at the structure of the arguments passed to this function

ind = slct(1).index; % for this case all indices should be the same

% --- Do something for fun --- %
if ~isempty(extrapnt) && isvalid(extrapnt)
    delete(extrapnt);
end
% Place a marker on the selected axes
extrapnt = line(ln.XData(ind), ln.YData(ind), 'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 10);

% If this point on the third plot is below the cutoff
if abs(ui.xpl(3).data{2,2}(ind)) < cutoff
    if ~isempty(newfig) && isvalid(newfig)
        close(newfig);
    end
    newfig = figure;
    surf(160*membrane(1,100), 'EdgeColor', 'none');
    axis tight
else
    disp(['Cutoff is ' num2str(cutoff)])
end

% --------------------------- %

% We can optionally update the arguments for this callback like so. The 
% user data should be organized in a cell array as with the initial call to
% exploreResults. The first three arguments to the callback are set by 
% exploreResults and normally don't need to be editted here.+

user_args = {@ userCallback, ln, extrapnt, cutoff, newfig};
src.Callback{end} = user_args;

end
