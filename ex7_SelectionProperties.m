%% ex7_SelectionProperties
% ----------------------------------------------------------------------- %
% With the 'SelectionProperties' option we can set visual properties for
% the selected point markers.
% ----------------------------------------------------------------------- %

%% Build a figure
FIG = figure;

subplot(1,2,1)

x = 0.05:0.05:2;
y = log(x).^2;

CHARTS(1) = plot(x,y,'-o');

title('Example Plot 1')
xlabel('x')
ylabel('log(x)^2')

subplot(1,2,2)

x = 0.05:0.05:2;
y = log(x).^2;

CHARTS(2) = plot(x,y,'-o');

title('Example Plot 2')
xlabel('x')
ylabel('log(x)^2')

yyaxis right
y = sqrt(x);

CHARTS(3) = plot(x,y,'-s');

ylabel('x^(1/2)')

%% Setup exploreResults
PUSHBUTTONS = {{'Display', @myFunc}};

USERDATA = {{}, ... % chart 1
            {'y^2', sqrt(x).^2; 'Index', 1:length(x)}}; % chart 3

% ----------------------------------------------------------------------- %
% This cell array should be a list of name/value pairs of marker visual
% properties. 
%   MARKERPROPS = {'name1', value1, 'name2', value2};
% The properties that can be changed are:
%   Marker, MarkerEdgeColor, MarkerFaceColor, MarkerSize
MARKERPROPS = {'Marker', 'd', ...
               'MarkerEdgeColor', [1 0 0], ...
               'MarkerFaceColor', 'none', ...
               'MarkerSize', 15};
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% Set 'SelectionProperties' to the cell array MARKERPROPS to control the
% visual settings of the marker.
exploreResults(FIG, PUSHBUTTONS, CHARTS([1 3]), 'DataFromUser', USERDATA, ...
    'SelectionLinkCharts', true, 'SelectionPerChart', 2, 'SelectionProperties', MARKERPROPS);
% ----------------------------------------------------------------------- %

%% Callback Functions
function myFunc(src, event, ui, slct)
% Display the data structs passed to the function by exploreResults. These
% contain information on the ui elements and selected points.
disp('The structure of the arguments holding selected point information')
disp('slct:')
disp(slct)
disp('ui:')
disp(ui)

disp(['The number of selectable charts is: ' num2str(length(ui.xpl))])
disp(['The number of selected points is: ' num2str(length(slct))])
fprintf('\n')

end
