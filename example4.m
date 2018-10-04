% Example script on how to setup and launch the data explore function

% First generate the main figure to serve as the ui window.
mainfig = figure;
plt = plot3(sin(1:0.1:10), cos(1:0.1:10), 1:0.1:10);
legend('Legend')
title('3D Plot')
xlabel('x')
ylabel('y')
zlabel('z')

ax = mainfig.Children;
ax(end).Tag = 'Explorable';
plt.Tag = 'Explorable'; 

% Push button call back
% This must be a cell array where the first entry is an anonymous function
% and the following are arguments for that function. See the documentation
% on this argument for details on the requirements of this function.
phrase = 'Display Me!'; % just an argument for the example function
pbtn_callback = {@ userCallback, phrase};

% Set two aditional options
usefigdat = true; % use data from figure for display boxes
linkselect = true; % if there are multiple plots select the same point on all of them (requires each line has the same number of points)

% Finally launch the ui figure
exploreResults( mainfig, pbtn_callback, [], usefigdat, linkselect );

%% --- Example User Function --- %%
function userCallback(src, event, slct, ui, phrase)
% An example of a function to assign to the ui push button. See
% documentation for details on the reserved input arguments for src, event,
% slct and ui.

disp(phrase);

end
