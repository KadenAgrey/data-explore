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
phrase = 'Display Me!'; % just an argument for the example function
pbtn_callback = {{'User Function', @ userCallback, phrase}};

% Finally launch the ui figure
exploreData( mainfig, pbtn_callback, sf, 'DataFromAxes', true, 'SelectionLinkCharts', false );

%% --- Example User Function --- %%
function [ newfig ] = userCallback(src, event, ui, slct, phrase)
% An example of a function to assign to the ui push button. See
% documentation for details on the reserved input arguments for src, event,
% slct and ui.

disp(phrase);

end
