% Example script on how to setup and launch the data explore function

% First generate the main figure to serve as the ui window.
mainfig = figure;
sf = surf(160*membrane(1,100), 'EdgeColor','none');
title('Surf Plot')
xlabel('x')
ylabel('y')
zlabel('z')

% Push button call back
% This must be a cell array where the first entry is an anonymous function
% and the following are arguments for that function. See the documentation
% on this argument for details on the requirements of this function.
pbtn_callback = { {'View x', @ userCallback, 'x'}, ...
                  {'View y', @ userCallback, 'y'}, ...
                  {'View z', @ userCallback, 'z'} };

% Finally launch the ui figure
exploreData( mainfig, pbtn_callback, sf, 'DataFromAxes', true, 'SelectionLinkCharts', false );

%% --- Example User Function --- %%
function userCallback(src, event, ui, slct, dim)
% An example of a function to assign to the ui push button. See
% documentation for details on the reserved input arguments for src, event,
% slct and ui.

if strcmp(dim, 'x')
    disp(['x = ' num2str(slct(1).point(1))]);
elseif strcmp(dim, 'y')
    disp(['y = ' num2str(slct(1).point(2))]);
elseif strcmp(dim, 'z')
    disp(['z = ' num2str(slct(1).point(3))]);
else
    disp('Please choose either ''x'', ''y'', ''z''');
end

end
