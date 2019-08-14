%% ex4_DataFromX
% ----------------------------------------------------------------------- %
% This example demonstrates how to control which data is displayed on the
% figure and passed through to the callback function. This is really
% helpful when you want to access information on a point not associated
% with one of the plot axes.
% 
% This is done with the name/value pair arguments 'DataFromUser' and
% 'DataFromAxes'.
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

% ----------------------------------------------------------------------- %
% A cell array with name/data pairs that will be displayed with the 
% datatips. The cell array should be nested so that the first index 
% corresponds to the chart list, and the next index to name/data pairs. The
% name/data pairs should be a cell array where the first column is the
% names and the second the data array.
% 
% Data arrays must have the same number of elements and dimension as the 
% chart from which it will be selected.
USERDATA = {{}, ... % chart 1
            {'y^2', sqrt(x).^2; 'Index', 1:length(x)}}; % chart 2
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% The 'DataFromUser' argument adds data series to be displayed along with
% the plotted data. 
% The 'DataFromAxes' argument disables including the plotted data in the
% display.

exploreResults(FIG, PUSHBUTTONS, CHARTS([1 3]), 'DataFromUser', USERDATA);
% exploreResults(FIG, PUSHBUTTONS, CHARTS([1 3]), 'DataFromAxes', false);
% exploreResults(FIG, PUSHBUTTONS, CHARTS([1 3]), 'DataFromUser', USERDATA, 'DataFromAxes', false);
% ----------------------------------------------------------------------- %

%% Callback Functions
function myFunc(src, event, ui, slct)
% Display information on the first selected point.
disp(['The first selected point comes from chart ' ...
      num2str(slct(1).chartnum) ' in ui.xpl with:'])
disp(['    point = ' num2str(slct(1).point')])
disp(['    index = ' num2str(slct(1).index)])

disp('The data cell in ui.xpl.data:')
disp(ui.xpl(slct(1).chartnum).data)

end
