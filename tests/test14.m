% Example script on how to setup and launch the data explore function

% First generate the main figure to serve as the ui window.
fig = figure;
fig.Position(3) = fig.Position(3)*1.7;

subplot(1,2,1);
ln = plot3(sin(linspace(1,10,100)), cos(linspace(1,10,100)), linspace(1,10,100));
legend('Legend')
title('3D Plot 1')
xlabel('x')
ylabel('y')
zlabel('z')

subplot(1,2,2);
ln = plot3(sin(linspace(1,10,100)), tan(linspace(0,pi/3,100)), cos(2*linspace(1,10,100)));
legend('Legend')
title('3D Plot 2')
xlabel('x')
ylabel('y')
zlabel('z')

% Push button call back
% This must be a cell array where the first entry is an anonymous function
% and the following are arguments for that function. See the documentation
% on this argument for details on the requirements of this function.
phrase = 'Display Me!'; % just an argument for the example function
pbtn_callback = {{'User Function', @ userCallback, phrase}};
markerprops = {'MarkerSize',5, 'MarkerFaceColor','none', 'MarkerEdgeColor','k', 'Marker','o'};

% Finally launch the ui figure
exploreData( fig, pbtn_callback, 'DataFromAxes', true, 'SelectionLinkCharts', true, 'SelectionPerChart', 1, 'SelectionProperties', markerprops );

%% --- Example User Function --- %%
function userCallback(src, event, ui, slct, phrase)
% An example of a function to assign to the ui push button. See
% documentation for details on the reserved input arguments for src, event,
% slct and ui.

disp(phrase);

end
