%% ex3_ChoosingSelectableCharts
% ----------------------------------------------------------------------- %
% This example demonstrates how be selective about which charts (lines) can
% be selected from a plot.
% ----------------------------------------------------------------------- %

%% Build a figure
% ----------------------------------------------------------------------- %
% Here we are storing the chart objects returned by plot() in an array. By 
% passing this as the third argument to exploreResults we can select which
% lines we want to be "selectable".
% ----------------------------------------------------------------------- %

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
% Including this arugment allows us to control which charts we want to be
% selectable. In this case passing all of CHARTS through is equivalent to
% the default functionality of exploreResults, which is to allow all lines
% to be selected.
exploreResults(FIG, PUSHBUTTONS, CHARTS([1 3]));
% exploreResults(FIG, PUSHBUTTONS);
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
disp('')

% Display information on the first selected point.
disp(['The number of selectable charts is: ' num2str(length(ui.xpl))])
disp(['The number of selected points is: ' num2str(length(slct))])
disp(['The first selected point comes from chart ' ...
      num2str(slct(1).chartnum) ' in ui.xpl with:'])
disp(['    point = ' num2str(slct(1).point')])
disp(['    index = ' num2str(slct(1).index)])
disp('')

disp('The data labels for this point are:')
str = '    ';
for d = 1:size(ui.xpl(slct(1).chartnum).data,1)
    str = [str '  ' ui.xpl(slct(1).chartnum).data{d,1}];
end
disp(str)
disp('')

disp('The full data arrays can be accessed by ''ui.xpl( slct(1).chartnum ).data''')

end
