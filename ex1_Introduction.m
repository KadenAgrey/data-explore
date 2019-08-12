%% ex1_Introduction
% ----------------------------------------------------------------------- %
% exploreResults is a helpful function designed to allow the user to select
% data points from a plot (or multiple plots) and run further analysis or
% visualization on them. 
% This example demonstrates two things required to run the function. 
%   1. the format to define pushbuttons and their callback functions 
%      (functions executed when they are pressed). 
%      See section: Setup exploreResults.
%   2. the order and structure of arguments that exploreResults passes to 
%      those functions and how to pass your own arguments.
%      See section: Callback Functions
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
% More about the reserved arguments can be found in the initial 
% documentation for exploreResults under the description of the input
% argument "pbtnfcn".
% 
% This function can be defined to do anything you'd like; compute results
% based on the selected points, display figures of data taken from the
% selected points, launch another instance of explore results, etc. 
% 
% This function cannot, unfortunately, return any data back to the 
% workspace easily. This is not how MATLAB has designed their graphics 
% objects to work. To return data to the workspace the user has two
% options. First is to save the data into a matlab file and load it again 
% later. Second is to place the data in the figures 'UserData' property. 
%   FIG.UserData = mydata;
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
% This is just a placeholder function used as an example for defining a
% pushbutton whose function takes user input (myarg).
% ----------------------------------------------------------------------- %

disp(myarg)

end