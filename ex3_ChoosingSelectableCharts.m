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
myarg1 = 'My String';
PUSHBUTTONS = {{'Display', @myFunc1}, ...
               {'Another Button', @myFunc2, myarg1}};

% ----------------------------------------------------------------------- %
% Including this arugment allows us to control which charts we want to be
% selectable. In this case passing all of CHARTS through is equivalent to
% the default functionality of exploreResults, which is to allow all lines
% to be selected.
exploreResults(FIG, PUSHBUTTONS, CHARTS([1 3]));
% exploreResults(FIG, PUSHBUTTONS);
% ----------------------------------------------------------------------- %

%% Callback Functions
function myFunc1(src, event, ui, slct)

display([slct.chart])

end

function myFunc2(src, event, ui, slct, myarg)

disp(myarg)

end