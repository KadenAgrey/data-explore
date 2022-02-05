%% ex5_SelectionLinkCharts
% ----------------------------------------------------------------------- %
% With the 'SelectionLinkCharts' option we can force the same point to be
% selected on all charts. 
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

%% Setup exploreData
PUSHBUTTONS = {{'Display', @myFunc}};

USERDATA = {{}, ... % chart 1
            {'y^2', sqrt(x).^2; 'Index', 1:length(x)}}; % chart 3

% ----------------------------------------------------------------------- %
% Setting 'SelectionLinkCharts' to true will force all charts to have the
% same index selected. WARNING: this only works when all charts plotted
% have the same number of points!
exploreData(FIG, PUSHBUTTONS, CHARTS([1 3]), 'DataFromUser', USERDATA, ...
    'SelectionLinkCharts', true);
% ----------------------------------------------------------------------- %

%% Callback Functions
function myFunc(src, event, ui, slct)
% Display the data structs passed to the function by exploreData. These
% contain information on the ui elements and selected points.
disp('The structure of the arguments holding selected point information')
disp('slct:')
disp(slct)
disp('ui:')
disp(ui)

disp(['The number of selectable charts is: ' num2str(length(ui.xpl))])
disp(['The number of selected points is: ' num2str(length(slct))])
fprintf('\n')

% ----------------------------------------------------------------------- %
% Display information on slct.links. This indicates which selected points
% are linked together via indicies.
if ~isempty(slct)
    disp('''slct.links'' is now filled out for each selected point.')
    disp('The first point is linked to the point(s) slct(1).links')
    disp(['    Point 1 is linked to point(s): ' num2str(slct(1).links)])
    disp(['    Point 2 is linked to point(s): ' num2str(slct(2).links)])
else
    disp('Please select a data point!')
end
% ----------------------------------------------------------------------- %

end
