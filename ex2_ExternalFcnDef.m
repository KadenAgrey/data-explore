%% ex2_ExternalFcnDef
% ----------------------------------------------------------------------- %
% This example simply demonstrates that the callback functions do not need
% to be defined within the same script that exploreData is called. They
% can be external as well.
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

%% Setup exploreData
% ----------------------------------------------------------------------- %
% Here we are using a callback function defined in an external function
% file. Going forward we will keep all callback functions defined within
% the script that launches exploreData simply for convenience.
PUSHBUTTONS = {{'Display', @ex2myFunc}};
% ----------------------------------------------------------------------- %

exploreData(FIG, PUSHBUTTONS);
