%% ex1_Introduction
% ----------------------------------------------------------------------- %
% This example demonstrates two things required to run the function. 
%   1. the format to define pushbuttons and their callback functions 
%      (functions executed when they are pressed). 
%      See section: Setup exploreResults.
%   2. the order and structure of arguments that exploreResults passes to 
%      those functions and how to pass your own arguments.
%      See section: Callback Function
% ----------------------------------------------------------------------- %

%% Build a figure
% Just an example plot demonstrate some functionality
FIG = figure;

x = 0.05:0.05:2;
y = log(x).^2;

plot(x,y,'-o')

title('Example Plot')
xlabel('x')
ylabel('log(x)^2')

%% Setup exploreResults
% ----------------------------------------------------------------------- %
% This argument defines the number, label, and callback functions plus
% arguments of pushbuttons on the figure window. exploreResults will the 
% pushbuttons automatically. 
% 
% Each element of the PUSHBUTTONS cell corresponds to one pushbutton and
% defines the button name, callback function handle, and user defined input
% arguments to the function as a nested cell array.
myarg1 = 'Button Pressed';
PUSHBUTTONS = {{'Display', @myFunc1}, ...
               {'Another Button', @myFunc2, myarg1}};
% ----------------------------------------------------------------------- %

exploreResults(FIG, PUSHBUTTONS)

%% Callback Functions
function myFunc1(src, event, ui, slct)
% ----------------------------------------------------------------------- %
% This is the function that will execute when the pushbutton is pressed. 
% The first 4 arguments to the function MUST be reserved for 
% exploreResults. The latter arguments should correspond to those in the 
% cell array. Eg:
%   pbtnfcn = {'My Button', @myFunc, arg1, arg2};
% will correspond to the user function,
%   function myFunc(src, event, ui, slct, arg1, arg2 )
% 
% The reserved arguments are:
%   src: matlab variable pointing to the src of the callback, which is
%   the pushbutton object.
%   event: matlab variable giving event information.
%   ui: exploreResults variable with handles for all ui elements and
%   explorable chart objects.
%       ui.dcm: data cursor manager object
%       ui.pbtn: push button object
%       ui.xpl(i).chart: explorable chart object (a line or contour etc.)
%       ui.xpl(i).data: cell array containing axes and/or user data
%       associated with the chart object.
%   slct: exploreResults variable with information on the selected
%   points. 
%       slct(j).chart: chart object point is on
%       slct(j).chartnum: index of chart object in ui.xpl(i).chart
%       slct(j).links: if SelectionLinkCharts is true, index of other 
%       selected points, slct(j), linked to this one.
%       slct(j).index: linear index of selected point in associated 
%       chart data, ui.xpl(chartnum).data.
%       slct(j).point: value of associated data at index.
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% Display the data structs passed to the function by exploreResults. These
% contain information on the ui elements and selected points.
display(slct)
display(ui)
% ----------------------------------------------------------------------- %

end

function myFunc2(src, event, ui, slct, myarg)
% ----------------------------------------------------------------------- %
% This is the function that will execute when the pushbutton is pressed. 
% The first 4 arguments to the function MUST be reserved for 
% exploreResults. The latter arguments should correspond to those in the 
% cell array. Eg:
%   pbtnfcn = {'My Button', @myFunc, arg1, arg2};
% will correspond to the user function,
%   function myFunc(src, event, ui, slct, arg1, arg2 )
% 
% The reserved arguments are:
%   src: matlab variable pointing to the src of the callback, which is
%   the pushbutton object.
%   event: matlab variable giving event information.
%   ui: exploreResults variable with handles for all ui elements and
%   explorable chart objects.
%       ui.dcm: data cursor manager object
%       ui.pbtn: push button object
%       ui.xpl(i).chart: explorable chart object (a line or contour etc.)
%       ui.xpl(i).data: cell array containing axes and/or user data
%       associated with the chart object.
%   slct: exploreResults variable with information on the selected
%   points. 
%       slct(j).chart: chart object point is on
%       slct(j).chartnum: index of chart object in ui.xpl(i).chart
%       slct(j).links: if SelectionLinkCharts is true, index of other 
%       selected points, slct(j), linked to this one.
%       slct(j).index: linear index of selected point in associated 
%       chart data, ui.xpl(chartnum).data.
%       slct(j).point: value of associated data at index.
% ----------------------------------------------------------------------- %

disp(myarg)

end