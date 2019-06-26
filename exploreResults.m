function [ fig ] = exploreResults( fig, pbtnfcn, varargin)
% exploreResults launches an interactive plot from which data points can be
% selected and detailed results for that case viewed.
% By: Kaden Agrey
% v1.1 2018.07.11
% 
% Input
%   mainplot: Figure handle on which to place ui elements. The script will 
%   automatically add space for the for the ui elements.
% 
% 
%   pbtnfcn: This argument is used to initialize the push buttons. It 
%   should be a cell array with the button name as the first column, a 
%   function handle as the second column and arguments to the function as 
%   later columns. Each row corresponds to a different pushbutton, which 
%   will be added to the figure as needed.
%   The first 4 arguments to the function MUST be reserved for
%   exploreResults. The latter arguments should correspond to those in the
%   cell array. Eg.
%       function argout = myFunc(src, event, slct, ui, arg1, arg2, ... )
%       src: matlab variable pointing to the src of the callback.
%       event: matlab variable giving event information
%       slct: exploreResults variable with information on the selected
%       points. 
% TODO: fields have been changed, update comment here
%           slct(i).xplobj: graphics object with explorable data
%           slct(i).pnt: list of graphics objects for selected points
%           slct(i).ind: list of selected linear indices from explorable 
%           data
%       ui: exploreResults variable with pointers to all ui elements
% 
% 
%   lines: A list of the graphics objects that will be selectable. This is
%   not a list of axes, but of the axes children that will be selectable.
%   Supported types are: lines (2D, 3D), scatters (2D, 3D), contours and
%   surfaces. Other types may work. These object options are set correctly 
%   by default in MATLAB.
%       ln.HitTest = 'on'
%       ln.PickableParts = 'visible'
% 
%   TODO: Set fig, axes, and line units carefully
%   TODO: manually positioned colorbars shuoldn't always be aligned with the
%   row they belong to - switch to an algorithm that offsets row members
%   instead of aligning them.
%   TODO: Find out how to force Matlab to wait for the renderer before
%   executing the mainplot formatting. - drawnow; doesn't seem to work
% 
% 
%   txtedtdat: A cell array with name/data pairs that will display
%   information on the selected points under each axes. The cell array
%   should be nested so that the first index corresponds to the axes number
%   and the second to the name/data pairs.
%       txtedtdat = { {name1, data1, name2, data2}, ... % axes 1
%                     {name1, data1} }                  % axes 2
%   Axes numbers should index in the same order that the "Explorable" axes
%   were created and the data given to each name must have the same number
%   of elements as the line from which it will be selected.
% 
% 
%   usefigdat: Additionally the user can allow the function to
%   automatically get the data from the axes and line objects by setting 
%   this option to true. This will create two text/edit pairs for each 
%   "Explorable" axes.
% 
% 
%   linkselect: When this option is set to true it will force all
%   axes with explorable children to have the same index selected. 
%   Selecting a point on one plot will select the same point on all plots. 
% 
%   TODO: Let user specify which plots to link
% 
% Change Log
%   2018.07.31: Added support for multiple rows of subplots
%   2018.08.15: Added support for plots generated from gridded data (like
%   contours) via linear index conversion.
%   2018.09.26: Set unexplorable lines to not be 'pickable'.
%   2019.03.09: Change from using Tag to indicate selectable lines to
%   passing a list of objects, this is more Matlab-like
%   2019.05.05: Added input parsing and param/value arguments
%   2019.05.10: Added support for multiple user push button functions

%% --- Parse Inputs --- %%
in = inputParser(); % initialize parser object

% validDataBoxMode = {'active', 'inactive'};
% checkDataBoxMode = @(s) validateString(s, validDataBoxMode);

% Required
in.addRequired('fig', @isgraphics); % figure handle to build ui on
in.addRequired('pbtnfcn'); % function handle (with arguments) to call when button is pressed
% Optional Positional
in.addOptional('lines',gobjects(0)); % lines to select data points from (make optional later)
% Optional Name/Value Pairs
in.addParameter('DataFromAxes',true); % display data from axes in boxes
in.addParameter('DataFromUser',{}); % display boxes will be added with user data
in.addParameter('SelectionLinkAxes',true); % link selected points accross all selectable objects
% in.addParameter('SelectionPerLine', inf); % number of data points that can be selected per axes
% in.addParameter('SelectionPerAxes', inf); % number of data points that can be selected per axes
% in.addParameter('SelectionProperties', []); % properties of selection marker
% in.addParameter('DataBoxMode', 'inactive', checkDataBoxMode); % allow manual entry of data into display boxes (will be passed to push buttons)

in.parse(fig, pbtnfcn, varargin{:})

% Pull fig out of the object for easier access
fig = in.Results.fig;

%% --- Initialize --- %%
% These variables define the size and spacing of the ui objects
pbtn_h = 30; % [pixels] PushButton height
pbtn_blw_marg = 5; % [pixels] margin below PushButton
pbtn_abv_marg = 3; % [pixels] margin above PushButton
% txtedt_h = 20*2; % [pixels] Text & Edit box height
% txtedtmarg = [3 10 0 10]; % [pixels] Text-Edit box margins
slctopt_fields = {'SelectionLinkAxes'};
% slctopt_fields = {'SelectionLinkAxes','SelectionPerLine',...
%                   'SelectionPerAxes'};

% Initialize xplr struct with figure objects and slct struct with graphics 
% objects to select points from. ui struct is also initialized, each field
% holds the ui objects by associated graphics object. (push buttons for the
% figure, text/edit boxes for each axes they are placed under etc.)

% Get all axes parenting explorable lines
xplr = struct('ax', [], 'ln', getExplorableLines(fig, in.Results.lines));
xplr.ax = getExplorableAx(fig, xplr.ln);

% Initialize the struct to reference ui objects, this is passed to 
% callbacks.
%   ui.pbtn is an array holding all pbtn objects
%   ui.dcm holds the datacursormode object
% Calling datacursormode here creates the mode object, whose default
% functions we will change to add our own functionality.
ui = struct('pbtn', [], 'dcm', datacursormode(fig));

% Get objects selected points will be associated with. This is passed to
% callbacks.
if isempty(in.Results.DataFromUser)
    user_data = cell(length(xplr.ln),1);
else
    user_data = in.Results.DataFromUser;
end

slct = struct('xplobj', [], 'data', []);
for p = 1:length(xplr.ln)
    slct(p).xplobj = xplr.ln(p);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ugly
    if isempty(slct(p).xplobj.ZData)
        slct(p).data = {slct(p).xplobj.Parent.XLabel.String, slct(p).xplobj.XData; ...
                        slct(p).xplobj.Parent.YLabel.String, slct(p).xplobj.YData};
    else
        slct(p).data = {slct(p).xplobj.Parent.XLabel.String, slct(p).xplobj.XData; ...
                        slct(p).xplobj.Parent.YLabel.String, slct(p).xplobj.YData; ...
                        slct(p).xplobj.Parent.ZLabel.String, slct(p).xplobj.ZData};
    end
    slct(p).data = [slct(p).data; user_data{p}];
end

% Define a struct containing selection options to pass to the callbacks
slctopt = struct();
for f = slctopt_fields
    slctopt.(f{:}) = in.Results.(f{:});
end

%% --- Prepair Main Figure --- %%
% --- Figure size and plot positioning --- %
% % TODO: clean up algorithm
% % Resize the figure and position the plots to fit the ui
% children = findobj(fig.Children, 'flat', '-not', 'AxisLocationMode', 'auto');
% nch = length(children); % number of children in figure
% if nch > 1
%     cbmap = cellfun(@(C) strcmp(C,'colorbar'), get(children, 'Type'));
% else
%     cbmap = strcmp(get(children, 'Type'),'colorbar');
% end
% fig.Units = 'pixels';
% set(fig.Children, 'Units', 'pixels');
% 
% % Get size of all the figure children
% ax_tins = zeros(nch,4);
% for ch = 1:nch
%     % Some objects don't have a TightInset
%     if isprop(children(ch), 'TightInset')
%         ax_tins(ch,:) = children(ch).TightInset;
%     else
%         ax_tins(ch,:) = [0 0 0 0];
%     end
% end
% 
% ax_pos = getAsMat(children, 'Position');
% ax_siz = ax_pos(:,3:4); % save plot sizes now - they can automagically change
% ax_sub = getSubPlotInd(ax_pos(~cbmap,:)); % subplot-like indices for children (cell array)
% ax_sub = assignCBSubPlotInd(children, ax_sub, cbmap);
% ax_row = cellfun(@max, ax_sub(:,1)); % for our purposes we only want the lowest row each plot is part of
% nrow = max(ax_row);
% 
% % Get rows that require padding (columns won't require padding)
% pad_row = []; % list of rows requiring padding
% for row = 1:nrow
%     rowch = find(row==ax_row)';
%     for ch = rowch
% 
%         % If the child is Explorable then add the row to the pad list
%         ind = children(ch)==xplr.ax;
%         if any(ind) && ( in.Results.DataFromAxes || ~isempty(in.Results.DataFromUser(ind)) )
%             pad_row(length(pad_row)+1,1) = row;
%             break
%         end
% 
%     end
% end
% 
% % Reposition children
% % TODO: manually positioned colorbars shouldn't always be aligned with the
% % row they belong to - switch to an algorithm that offsets row members
% % equally instead of aligning them.
% ax_pos_new = ax_pos; % initially nothing has moved
% for row = nrow:-1:1 % start from the bottom
%     % Add text/edit padding if this row requires it
%     if ismember(row, pad_row)
%         axpad = txtedt_h + sum(txtedtmarg([2 4]));
%     else
%         axpad = txtedtmarg(2);
%     end
% 
%     % Find the axes from which to get the height this row should start at
%     rowch = find(row==ax_row); % children in this row
%     rowcol = unique( cell2mat(ax_sub(rowch,2)) );
% 
%     rowt = max(ax_tins(rowch,:),[],1); % the largest TightInsets from this row (we only use the bottom and top)
%     rowy = ax_pos(rowch(1),2) - rowt(2); % lowest y-position of this row
% 
%     % Note: here a 10 pixel overlap is allowed when checking if a figure is
%     % below the row
%     allsharecol = cellfun(@(C) any(ismember(C, rowcol)), ax_sub(:,2));
%     allblw = ( sum(ax_pos(:,[2 4]), 2) + ax_tins(:,4) - 10 < rowy );
%     chblw = find( allsharecol.*allblw ); % children below and sharing a column with this row
% 
%     if ~isempty(chblw)
%         rowy = max( sum(ax_pos_new(chblw,[2 4]), 2) + ax_tins(chblw,4) ); % new y-position of this row
%     else
%         rowy = pbtn_h - txtedtmarg(2) + 5; % if this is the bottom row place it 5 pixels above the PushButtons
%     end
% 
%     % Move all the axes of this row accordingly
%     for ch = rowch'
%         children(ch).Position(2) = rowy + rowt(2) + axpad;
%     end
% 
%     % Update the axes positions recorded
%     ax_pos_new = getAsMat(children, 'Position'); % ################## CHNAGE LATER
%     ax_pos_new(:,3:4) = ax_siz; % maintain axes sizes
% 
% end
% 
% % Resize figure window
% fig.OuterPosition(2) = 50; % move to the bottom of the screen
% fig.Position(4) = max( sum(ax_pos_new(:,[2 4]), 2) + ax_tins(:,4) ) + 10;
% 
% % Reset the position sizes incase they have magically changed
% for ch = 1:length(children)
%     children(ch).Position(3:4) = ax_pos_new(ch, 3:4);
% end

fig.Units = 'pixels';
set(fig.Children, 'Units', 'pixels');

fig.Position(4) = fig.Position(4) + pbtn_h + pbtn_blw_marg + pbtn_abv_marg;
children = findobj(fig.Children, 'flat', '-not', 'AxisLocationMode', 'auto');
for ch = children'
    ch.Position(2) = ch.Position(2) + pbtn_h + pbtn_blw_marg + pbtn_abv_marg;
end

% --- Make non-explorable axes children 'unpickable' --- %
% Set the 'PickableParts' property for non-explorable axes children to
% 'none'.
% We only need to consider axes that are marked explorable
for a = 1:length(xplr.ax)
    for ch = 1:length(xplr.ax(a).Children)
        if ~ismember(xplr.ax(a).Children(ch), xplr.ln)
            xplr.ax(a).Children(ch).PickableParts = 'none';
        end
    end
end

%% --- Setup UI --- %%
% uispec = cell(1,length(xplr.ax));

% Make button to view details of selected point
horz = getLeftChild(fig, 'pixels'); % place in line with farthest left plot

pbwidth = 150;
ui.pbtn = gobjects(size(in.Results.pbtnfcn,1),1);
for pb = 1:size(in.Results.pbtnfcn,1)
    ui.pbtn(pb) = uicontrol(fig, 'Style', 'pushbutton',... % make the button
                                 'String', in.Results.pbtnfcn{pb,1},...
                                 'Units', 'pixels',...
                                 'Position', [0 0 pbwidth pbtn_h]);

    ui.pbtn(pb).Position = ui.pbtn(pb).Position + [horz + pbwidth*(pb-1) pbtn_blw_marg 0 0];
    % ui.fig.pbtn(pb).Units = 'normalized'; % reset units
end

% Make text and edit objects for each explorable axes
% nax = length(xplr.ax);
% for a = 1:nax
% 
%     % Initial text and edit box size and options
%     txtopt = {'Position', [0 0 75 txtedt_h/2]};
%     editopt = {'Position', [0 0 75 txtedt_h/2], 'String', 'Empty', 'Enable', 'inactive'};
% 
%     if in.Results.DataFromAxes % user setting to use figure axes data
%         axdat = getAxesData(xplr.ax(a), xplr.ln);
%         if isempty(in.Results.DataFromUser) || isempty(in.Results.DataFromUser{a})
%             uispec{a} = axdat;
%         else
%             uispec{a} = [axdat, in.Results.DataFromUser{a}];
%         end
%     end
% 
%     % Make the text/edit pair
%     [ui.ax(a).txt, ui.ax(a).edt] = makeValueDispBar(xplr.ax(a), uispec{a}, txtedtmarg, txtopt, editopt);
% 
% end

% Assign user callback functions to push buttons
for pb = 1:size(in.Results.pbtnfcn,1)
    ui.pbtn(pb).Callback = {@ pbtnCallback, ui, slct, in.Results.pbtnfcn(pb,2:end)};
end

% Assign listener to dcm 'Enable' property
% plistenE = addlistener(ui.dcm,'Enable','PostSet',{@setTipCallbacks, fig, ui, slct, slctopt, {}});

% Assign default callbacks to dcmode
dcmode = getuimode(fig, 'Exploration.Datacursor');
% Set button down function
fcn = {@lnSelectPnt, dcmode.WindowButtonDownFcn, {}, ui, slct, slctopt};
dcmode.WindowButtonDownFcn = fcn;
% Set button up function
fcn = {@mvLinkedTipsButtonUp, dcmode.WindowButtonUpFcn, ui.dcm, {}};
dcmode.WindowButtonUpFcn = fcn;
% Set key press function
fcn = {@mvLinkedTipsKeyPress, dcmode.KeyPressFcn, ui.dcm, {}};
dcmode.KeyPressFcn = fcn;
% Set stop mode function - not needed if we continuously update mode
% callback args.
% fcn = {@updateArgsOnStop, dcmode.ModeStopFcn, dcmode};
% dcmode.ModeStopFcn = fcn;

% Assign UpdateFcn to data cursor mode object
ui.dcm.UpdateFcn = {@dcmUpdate, ui, slct};

%% --- Finalize Figure Properties --- %%
ui.dcm.Enable = 'on'; % turn data cursor mode on

% Reset units on figure and all figure children
fig.Units = 'normalized';
set(fig.Children, 'Units', 'normalized');

end

%% --- Functions --- %%
% --- Getters --- %
function [lns] = getExplorableLines(fig, lns)
% Gets axes for line/contours/surfaces etc. passed as arguments to
% exploreResults and checks that they all belong to the figure <fig>.

if isempty(lns)
    ch = findobj(fig.Children, 'flat', '-not', 'Type', 'colorbar');
    lns = gobjects(length(ch),1);
    for c = 1:length(ch)
        tmp = findobj(ch(c).Children, 'flat', '-not', 'Type', 'text');
        lns(c) = tmp(end);
    end
end

end

function [axs] = getExplorableAx(fig, lns)
% Gets axes for line/contours/surfaces etc. passed as arguments to
% exploreResults and checks that they all belong to the figure <fig>.

axsall = gobjects(length(lns),1); % holds an axes for each ln, some will be duplicates
for n = 1:length(lns)
    if ~isequal(ancestor(lns(n), 'figure'), fig)
        error(['The gobject ' lns(n).DisplayName 'is not part of the main ' ...
               'figure. All exploreable lines must be on the same figure'])
    end

    axsall(n) = ancestor(lns(n), 'axes');
end
axssort = unique(axsall); % No duplicates but we must order these correctly

% Loop over all children of parent figure
axs = gobjects(length(axssort),1);
ca = 1; % counter
for a = length(fig.Children):-1:1
% We loop backwards because figure children are in reverse order of 
% creation, and it's easiest for the user to require their data is 
% organized in order of creation.
    if ismember(fig.Children(a),axssort)
        axs(ca) = fig.Children(a);
        ca = ca + 1;
    end
end

end

function [figdat] = getAxesData(ax, lns)
% Gets explorable (selectable) data and labels from an axes. <lns> is a
% list of selectable line objects.

c = 1; % counter
figdat = cell(1,1);
for n = 1:length(ax.Children) % loop over lines on the axes
    ln = ax.Children(n);
    if ismember(ln, lns)
        % Get x and y data
        if any(strcmp(ln.Type, {'contour','surface'})) && isvector(ln.XData)
            [x,y] = meshgrid(ln.XData,ln.YData);
            figdat(c:c+3) = {ax.XLabel.String, x, ax.YLabel.String, y};
        else
            figdat(c:c+3) = {ax.XLabel.String, ln.XData, ax.YLabel.String, ln.YData};
        end
        c = c+4;
        % Get z data if present
        if ~isempty(ln.ZData)
            figdat(c:c+1) = {ax.ZLabel.String, ln.ZData};
            c = c+2;
        end
    end
end

end

function [subch] = getSubPlotInd(ch_pos)
% Gets the subplot row and column for the graphical object positions passed
% in <ch_pos>.

nch = size(ch_pos,1);
subch = cell(nch,2); % the [row col] pair each child belongs to

% Find unique rows and the row each child is in
row_pos = []; % position of rows
[~, sortind] = sort(ch_pos(:,2), 'descend'); % descending order is required for the following algorithm
for ch = sortind'
    % A row is 'new' if it hasn't been indexed yet, or if the vertical 
    % position is less than the current row.
    if ~ismembertol(ch_pos(ch,2), row_pos) && ( isempty( row_pos ) || all( ch_pos(ch,2)<row_pos ) )
        row_pos(length(row_pos)+1,1) = ch_pos(ch,2);
        subch{ch,1} = find( sum(ch_pos(ch,[2 4])) > row_pos );
    else
        subch{ch,1} = find( sum(ch_pos(ch,[2 4])) > row_pos );
    end
end

% Find unique rows and the row each child is in
col_pos = []; % position of cols
[~, sortind] = sort(ch_pos(:,1), 'descend'); % descending order is required for the following algorithm
for ch = sortind'
    % A column is 'new' if it hasn't been indexed yet, or if the horizontal 
    % position greater than the current column.
    if ~ismembertol(ch_pos(ch,1), col_pos) && ( isempty( col_pos ) || all( ch_pos(ch,1)<col_pos ) )
        col_pos(length(col_pos)+1,1) = ch_pos(ch,1);
        subch{ch,2} = find( sum(ch_pos(ch,[1 3])) > col_pos );
    else
        subch{ch,2} = find( sum(ch_pos(ch,[1 3])) > col_pos );
    end
end

% Reverse order of columns, the above algorithm requires them to be found
% in the wrong order.
ncol = max(cell2mat(subch(:,2)));
cols = ncol:-1:1;
for ch = 1:size(subch,1)
    subch{ch,2} = cols(subch{ch,2})';
end

end

function [pmat] = getAsMat(gobjs, prop)
% Returns the property requested of all graphics objects passed in as a
% matrix. (only works if the property can be stored as a matrix)

if length(gobjs) > 1
    pmat = cell2mat(get(gobjs, prop));
else
    pmat = get(gobjs, prop);
end

end

function [pos, ind] = getLeftChild(parent, unit)
% Gets leftmost child object
nch = length(parent.Children);
pos = zeros(nch,1);
for ch = 1:nch
    oldunit = parent.Children(ch).Units; % to reset units
    parent.Children(ch).Units = unit; % set units

    pos(ch) = parent.Children(ch).Position(1); % get horizontal position

    parent.Children(ch).Units = oldunit; % reset units
end
[pos,ind] = min(pos);

end

function [pos, ind] = getBottomChild(parent)
% Gets lowest child object
nch = length(parent.Children);
pos = zeros(nch,1);
for ch = 1:nch
    oldunit = parent.Children(ch).Units; % to reset units
    parent.Children(ch).Units = unit; % set units

    pos(ch) = parent.Children(ch).Position(2); % get horizontal position

    parent.Children(ch).Units = oldunit; % reset units
end
[pos,ind] = min(pos);

end

function [lnkdt] = getLinkedTips(fig)
% Gets linked tips cell array from fig window button down function
% arguments.

dcmode = getuimode(fig, 'Exploration.Datacursor');
if ~isempty(dcmode)
    lnkdt = dcmode.WindowButtonDownFcn{3};
else % if the figure is being closed return an empty cell
    lnkdt = {};
end

end

% --- Setters --- %
function [] = setTipCallbacks(src, event, fig, ui, slct, slctopt, lnkdt)
% Assign button down callbacks for figure window
disableFigFcnListener(fig)

% Set button down function
wndwfcn = {@lnSelectPnt, fig.WindowButtonDownFcn, lnkdt, ui, slct, slctopt};
fig.WindowButtonDownFcn = wndwfcn;
% Set button up function
wndwfcn = {@mvLinkedTipsButtonUp, fig.WindowButtonUpFcn, ui.dcm, lnkdt};
fig.WindowButtonUpFcn = wndwfcn;
% Set key press function
wndwfcn = {@mvLinkedTipsKeyPress, fig.KeyPressFcn, ui.dcm, lnkdt};
fig.KeyPressFcn = wndwfcn;

end

function [] = setLinkedTips(fig, lnkdt)
% Update lnkdt argument in figure and mode callbacks

setLinkedTipsFig(fig, lnkdt);
dcmode = getuimode(fig, 'Exploration.Datacursor');
setLinkedTipsMode(dcmode, lnkdt);

end

function [] = setLinkedTipsFig(fig, lnkdt)

disableFigFcnListener(fig)
fig.WindowButtonDownFcn{end}{3} = lnkdt;
fig.WindowButtonUpFcn{end}{4} = lnkdt;
fig.KeyPressFcn{end}{4} = lnkdt;

end

function [] = setLinkedTipsMode(mode, lnkdt)

mode.WindowButtonDownFcn{3} = lnkdt;
mode.WindowButtonUpFcn{4} = lnkdt;
mode.KeyPressFcn{4} = lnkdt;

end

% --- Makers --- %
function [dc] = makeDataCursor(dcm, target, ind, properties)
 
% Get the cursor position in data units
x = target.XData(ind);
y = target.YData(ind);
%     if ~isempty(target.ZData)
%         z = target.ZData(ind);
%     else
%         z = 1;
%     end

% Get cursor position in pixels
units = get( target.Parent, 'Units' ); % store fig units
set( target.Parent, 'Units', 'pixels' ); % set to pixels
[figpnt(1), figpnt(2)] = data2fig(target.Parent, x, y); % get position
set( target.Parent, 'Units', units ); % % reset units

dc = dcm.createDatatip(target, figpnt);

% % Create a copy of the context menu for the datatip:
% set(dc,'UIContextMenu',dcm.UIContextMenu);
% set(dc,'HandleVisibility','off');
% set(dc,'Host',target);
% set(dc,'ViewStyle','datatip');
% 
% % Set the data-tip orientation to top-right rather than auto
% set(dc,'OrientationMode','manual');
% set(dc,'Orientation','top-right');
% 
% % Update the datatip marker appearance
% set(dc, 'MarkerSize',5, 'MarkerFaceColor','none', ...
%     'MarkerEdgeColor','k', 'Marker','o', 'HitTest','off');
% 
% dc.update([x,y,1; x,y,-1]);
% dcm.updateDataCursors
% dcm.editUpdateFcn

end

% --- Utility --- %
function [] = chkExtent(txt)
% Checks that the extent of the string in a text object does not exceed the
% size. Because the 'Extent' property doesn't consider automatic string
% wrap this function doesn't either.

unit = txt.Units;
txt.Units = 'pixels';

% Width first - get size difference
dif3 = txt.Extent(3) - txt.Position(3) + 3; % Include some padding
% Adjust box size accordingly
if dif3 > 0
    txt.Position(3) = txt.Position(3) + dif3;
end
% Set the height likewise
dif4 = txt.Extent(4) - txt.Position(4) + 3;
if dif4 > 0
    txt.Position(4) = txt.Position(4) + dif4;
end

txt.Units = unit; % reset units

end

function [] = disableFigFcnListener(fig)
% This uses undocumented functionality, see link below if it breaks. We
% need to diable some listeners so that we can change the button down
% function of the figure and interject with our own code first. 
% https://undocumentedmatlab.com/blog/enabling-user-callbacks-during-zoom-pan
hManager = uigetmodemanager(fig);
try
    set(hManager.WindowListenerHandles, 'Enable', 'off');  % HG1
catch
    [hManager.WindowListenerHandles.Enabled] = deal(false);  % HG2
end

end

function [pnt] = ind2pnt(line, data, ind)
% Converts linear indices into the corresponding data point. This function
% accounts for gridded data. If the data are from a mesh grid object
% (surface or volumetric) then x, y, and z may be given as vectors. If this
% is the case it's assumed this spatial information comes first, and that
% the first non-vector data is the levels.

pnt = zeros(1,length(data)); % holds the point to return

% If the line doesn't use meshgrids then ind applies directly and we return
% early.
if ~any(strcmp(line.Type, {'contour', 'surface'}))
    pnt = cellfun(@(C) C(ind), data);
    return;
end

% If there is potential for gridded data then we need to check whether
% spatial information is gridded or still a vector. In either case we will
% return a vector with an index for each entry in <data>.
isvec = cellfun(@(C) isvector(C),data);
if any(isvec)
% The spatial data is in vector form. Note that this assumes the user
% provides only data with the same dimensions as the level data.
    % Get subscript indices of point to return.
    levi = find(~isvec,1); % index of 'level' data
    sub = zeros(1,length(data));
    if levi < 4
        [sub(1), sub(2)] = ind2sub(size(data{levi}), ind);
    else
        [sub(1), sub(2), sub(3)] = ind2sub(size(data{levi}), ind);
    end
    sub(levi:end) = ind;

else
% The spatial data is in matrix form. Note that this assumes the user
% provides only data with the same dimensions as the level data.
    sub = zeros(1,length(data));
    sub(:) = ind;
end

for p = 1:length(data)
    pnt(p) = data{p}(sub(p));
end

end

function [ind] = pnt2ind(line, pnt)
% Uses X, Y, and ZData in a line to find the associated index of a given
% point. Will always return the linear index matching the ZData, even if X
% and Y are vectors.

if strcmp(line.Type, {'contour', 'surface'})
% Data is a grid
    % We must search the X and Y data for subscripts because ZData may not 
    % be unique.
    % Get x index
    if isvector(line.XData) % XData is a vector
        x = find(line.XData == pnt(1), 1);
    else % XData is a nd grid
        x = find(line.XData(1,:) == pnt(1), 1);
    end
    % Get y index
    if isvector(line.YData) % XData is a vector
        y = find(line.YData == pnt(2), 1);
    else % XData is a nd grid
        y = find(line.YData(1,:) == pnt(2), 1);
    end

    % Convert to linear index
    ind = sub2ind(size(line.ZData),[y x]);
else
% Data is a 1D sequence.
    ind = find(line.XData == pnt(1), 1);
end

end

function [xf, yf] = data2fig(ax, x, y)
% Converts data coordinates (<x> and <y>) on an axes to equivalent
% coordinates in the figure with the same units as <ax>.
 
pos = ax.Position;
 
xf = pos(1) + pos(3)*(x - ax.XLim(1))/diff(ax.XLim);
yf = pos(2) + pos(4)*(y - ax.YLim(1))/diff(ax.YLim);
 
end

function [] = rmCursors(~, ~, dcm, dc)
% Invokes the datacursormanager.removeCursor function to remove an array of
% cursors.
%
% Input
% ~: 
%   First two arguments are placeholders for when this function is set as a
%   callback.

for c = dc(:)'
    dcm.removeDataCursor(c)
end

end

function [] = rmCursorsNoLink(~, ~, dcm, dt)
% Invokes the datacursormanager.removeCursor function to remove an array of
% cursors after disabling the delete function

set(dt, 'DeleteFcn', []);
rmCursors([], [], dcm, [dt.Cursor])

end

function [] = mvCursors(~,~, dc, curcur)
% Invokes the datacursor.moveTo function to move an array of cursors to the
% same index on their sources as curcur
%
% Input
% ~: 
%   First two arguments are placeholders for when this function is set as a
%   callback.

ind = curcur.DataIndex;
for c = dc(:)'
    x = c.DataSource.XData(ind);
    y = c.DataSource.YData(ind);
    if ~isempty(c.DataSource.ZData)
        z = c.DataSource.ZData(ind);
    end

    units = c.DataSource.Parent.Units;
    c.DataSource.Parent.Units = 'pixels';

    [fx, fy] = data2fig(c.DataSource.Parent, x, y);
    c.moveTo([fx, fy]);

    c.DataSource.Parent.Units = units;
end

end

function [] = mvLinkedTips(linkedtips, curtip)
% Moves all tips linked to the current cursor to the same index on their
% source.

cind = cellfun(@(C) any(ismember(C, curtip)), linkedtips); % cell index of tips linked to curtip

if any(cind)
    % Is DataIndex equal for the tips?
    DataIndices = arrayfun(@(A) A.DataIndex, [linkedtips{cind}.Cursor]);
    if ~all( DataIndices == DataIndices(1) )
        rtind = linkedtips{cind} ~= curtip;
        mvCursors([], [], [linkedtips{cind}(rtind).Cursor], curtip.Cursor); % move tips linked to curtip
    end

end

end

function [] = mvSrcLinkedTips(linkedtips, curtip, xplobj)
% If the source of the current cursor has changed then its linked tips are 
% moved to the remaining sources. This algorithm is fast enough for
% reassigning a few tips, but can be made much more efficient if required.

cind = cellfun(@(C) any(ismember(C, curtip)), linkedtips); % cell index of tips linked to curtip

if any(cind)
    % Are all explorable sources covered?
    Sources = arrayfun(@(A) A.DataSource, [linkedtips{cind}.Cursor], 'UniformOutput', false);
    cursrc = curtip.Cursor.DataSource; % source of current tip

    % Find which xplobj's have been used
    xplUsed = false(length(xplobj),1);
    for m = 1:length(xplobj)
        x = xplobj(m);

        for S = Sources
            if x{1} == S{1}
                xplUsed(m) = true;
                break;
            end
        end
    end

    curxpl = cellfun(@(C) isequal(C, cursrc), xplobj); % current xplobj
    if ~all(xplUsed) % if all xplobj aren't used
        xplUsed = curxpl; % reset xplUsed
        curlnk = linkedtips{cind} == curtip; % find current tip in linkedtips
        % Assign cursors to remaining data sources
        for tip = linkedtips{cind}(~curlnk)
            x = find(~xplUsed,1); % location of first unused src
            tip.DataSource = xplobj{x};
            xplUsed(x) = true;
        end
    end

end

end

% --- Callbacks --- %
function [] = pbtnCallback(src, event, ui, slct, usrcall)
% Callback function for the user push buttons. Takes the currently
% selected point information and the user supplied anonymous function with
% arguments. After checking that points have been properly selected will
% launch the user function with the first four arguments as matlab defined
% <src>, <event>, and exploreResults defined <usrobj> and <usrslct>.
% 
% TODO: What should be done if no points are selected? What about if less
% points than explorable objects have been selected?
% 
% TODO: Currently haven't decided how to handle linked data tips here. I'm
% thinking of changing the whole system over to a list of indices to make
% things easier here and in the figure callbacks.

% Define the external information structures. These are designed to make it
% easier to access selected point and all associated data.

% Contains ui graphics objects and lines associated with displayed data.
usrobj = struct('xpl', slct, 'dcm', ui.dcm, 'pbtn', ui.pbtn);
usrslct = struct('line', [], 'linenum', [], 'links', [], 'index', [], 'point', []);

cinfo = ui.dcm.getCursorInfo;
for p = 1:length(cinfo)
    lnnum = find( cellfun(@(C) isequal(C, cinfo(p).Target), {slct.xplobj}) , 1 ); % index of line in list of explorable objects

    % Get Target .(line)
    usrslct(p).line = cinfo(p).Target;
    usrslct(p).linenum = lnnum;

    % DataIndex .(index) and from cursor info if available. Else find it
    % from cursor info .(Position).
    if isfield(cinfo, 'DataIndex')
        usrslct(p).index = cinfo(p).DataIndex;
    else
        usrslct(p).index = pnt2ind(cinfo(p).Target, cinfo(p).Position);
    end

    % Get .(point) manually
    usrslct(p).point = ind2pnt(cinfo(p).Target, slct(lnnum).data(:,2), usrslct(p).index);
end

% Call the user's function and pass arguments through
usrcall{1}( src, event, usrobj, usrslct, usrcall{2:end} );

end

function [] = lnSelectPnt(fig, event, fcn, lnkdt, ui, slct, slctopt)
% Callback assigned to each selectable object (line, contour, surface
% etc.). This is interjected before the standard matlab datatip mode
% callback is executed and allows us to manage the datatips directly.

% Run the original callback
if isa(fcn,'function_handle')
   fcn(fig, event);
end

% If an explorable object was hit we need to continue with "normal" or 
% "extend" if the modifier is "shift" or "alt" as this means a new data 
% cursor was created.
% Note: undocumented event property "HitObject"
isxplhit = any(arrayfun( @(S) isequal(S.xplobj, event.HitObject), slct ));
% Copied from %matlabroot%/toolbox/matlab/graphics/
%               datacursormanager.m@localWindowButtonDownFcn()
mod = get(fig,'CurrentModifier');
isAddRequest = numel(mod)==1  && (strcmp(mod{1},'shift') || strcmp(mod{1},'alt'));
if ~isxplhit || ~( strcmp(fig.SelectionType,'normal') ...
        || (isAddRequest && strcmp(fig.SelectionType,'extend')) )
    return;
end

% Store the current cursor to reset after
curcur = ui.dcm.CurrentCursor;
alldt = findall(fig.Children, 'Type', 'hggroup');

% Determine if a new cursor was added
if ~isAddRequest && slctopt.SelectionLinkAxes
    if numel(alldt) < numel(slct) || ( numel([lnkdt{:}]) < numel(alldt) )
        isAddRequest = true;
    end
end

% Link Axes
if slctopt.SelectionLinkAxes && isAddRequest
    % Get index of selected object and point
    s = arrayfun( @(S) isequal(S.xplobj, ui.dcm.CurrentCursor.DataSource), slct );
    ind = ui.dcm.CurrentCursor.DataIndex;

    % Make Linked data tips
    cinfo = getCursorInfo(ui.dcm);
    xlns = {slct.xplobj}; % explorable lines
    clns = [cinfo.Target]; % cursor lines
    % This tip will be overwritten, I just want to initialize with the correct data type
    newdt = gobjects(0);
    for ln = xlns(~s)
        n = clns==ln{1};
        % Are there no tips on ln or are the tips at different indices?
        if all(~n) || all([cinfo(n).DataIndex] ~= ind)
            newdt(length(newdt)+1) = makeDataCursor(ui.dcm, ln{1}, ind, []);
        end
    end
    alldt = [alldt; newdt'];

    % Update lnktips
    curdt = findobj(alldt,'Cursor',curcur);
    lnkdt(length(lnkdt) + 1) = {[curdt newdt]};

    % Add callbacks so that all linked tips are deleted together and moved
    % together etc.
    for t = 1:length(lnkdt{end})
        rt = 1:length(lnkdt{end}) ~= t; % indices of other tips to remove
        % Add delete functions to linked tips
        lnkdt{end}(t).DeleteFcn = {@rmLinkedCursors, ui.dcm, lnkdt{end}(rt), fig};
        % Add property listener to move tips together - listeners are not
        % implimented for Datatips yet.
%         plistenP = addlistener(dt(t),'Position','PostSet',{@rmCursorsNoLink, ui.dcm, dt(rt)});
%         plistenS = addlistener(dt(t),'DataSource','PostSet',{@rmCursorsNoLink, ui.dcm, dt(rt)});
    end

    % Reset current cursor
    ui.dcm.CurrentCursor = curcur;

    % Update callback args
    setLinkedTips(fig, lnkdt)

elseif slctopt.SelectionLinkAxes
    % Ensure datatips move together
    curdt = findobj(alldt,'Cursor',curcur);
    mvSrcLinkedTips(lnkdt, curdt, {slct.xplobj})
    % We don't need to move the tips here because this is handled in the
    % button up function. We only ensure data sources are correct.
end

% Remove extra datatips
if false
    % NOT YET IMPLIMENTED
end

end

function [] = rmLinkedCursors(srcdt, ~, dcm, dt, fig)
% Delete function callback when cursors are linked. Removes linked tips and
% updates lnkdt callback argument for fig button down, button up, and key 
% press functions.

    rmCursorsNoLink([], [], dcm, dt)

    % To avoid needing to update every delete callback for every datatip we
    % will get the lnkdt variable from the mode callback arguments.
    lnkdt = getLinkedTips(fig);
    cind = cellfun(@(C) any(ismember(C, srcdt)), lnkdt);
    lnkdt = lnkdt(~cind);

    setLinkedTips(fig, lnkdt)
end

function [] = mvLinkedTipsButtonUp(src, event, fcn, dcm, lnkdt)
% Callback to call mvLinkedTips with appropriate arguments.

% Run the original callback
if isa(fcn,'function_handle')
   fcn(src, event);
end

% % Execute the specified callback function
% hgfeval(newButtonUpFcn,hFig,evd);

% If a tip was hit and selection type is normal continue.
% Note: undocumented event property "HitObject"
% isTipHit = strcmp(event.HitObject.Tag, 'PointTipLocator');
if ~strcmp(src.SelectionType,'normal')
    return;
end

% Get the current datatip
curdt = findobj([lnkdt{:}],'Cursor',dcm.CurrentCursor);

% Ensure linked tips move with current tip
mvLinkedTips(lnkdt, curdt)

end

function [] = mvLinkedTipsKeyPress(src, event, fcn, dcm, lnkdt)
% Callback to call mvLinkedTips with appropriate arguments.

% Run the original callback
if isa(fcn,'function_handle')
   fcn(src, event);
end

% Exit early if invalid event data
if ~isobject(event) || ~isvalid(event)
    return;
end

% Do nothing if an arrow key wasn't pressed.
if ~strcmp(event.Key, {'leftarrow','rightarrow','uparrow','downarrow'})
    return;
end

% Get the current datatip
curdt = findobj([lnkdt{:}],'Cursor',dcm.CurrentCursor);

% Ensure linked tips move with current tip
mvLinkedTips(lnkdt, curdt)

end

function str = dcmUpdate(pdtobj, eobj, ui, slct)

% Get index of selected object
s = [slct.xplobj] == eobj.Target;
ind = pdtobj.Cursor.DataIndex;

% Get data point
pnt = ind2pnt(eobj.Target, slct(s).data(:,2), ind);

% Make string of data to display
str = cell(1,length(slct(s).data));
try % >= 20XXx
    ui.dcm.Interpreter;
    for d = 1:length(slct(s).data)
        str{d} = [slct(s).data{d,1} ' {\bf\color{DarkGreen}{' num2str(pnt(d),4) '}}'];
    end
catch % <= 20XXx
    pdtobj.TextColor = [0.3 0.6 0.3];
    pdtobj.FontWeight = 'bold';
    for d = 1:length(slct(s).data)
        str{d} = [slct(s).data{d,1} ' ' num2str(pnt(d),4)];
    end
end

end

function [] = updateArgsOnStop(fcn, dcmode)
% Updates lnkdt arguments in uimode functions before disabling the mode so
% that they are correct when the mode is enabled again.

setLinkedTipsMode(dcmode, getLinkedTips(dcmode.FigureHandle));

fcn{1}(fcn{2:end});

end

%% --- Old Functions --- %%
function [txt, edt] = makeValueDispBar(ax, tedat, temarg, varargin)
% Makes a bar of labeled edit fields under the given axes

% Handle optional arguments
%             { txtopt, editopt }
default_opt = { {}    , {}      };
% Overwrite defaults
default_opt(1:length(varargin)) = varargin;
% Assign to pretty variable names
[txtopt, editopt] = default_opt{:};

% Make text/edit pairs for each axes
txt = gobjects(1,length(tedat)/2); edt = txt; % preinitialize
for n = 1:length(tedat)/2

    % Make edit box label
    if n > 1 % place next to last text box
        txt(n) = makePairedText(txt(n-1), temarg, [txtopt, {'String', tedat{2*n-1}}]);
    else % else place under axes
        txt(n) = makePairedText(ax, temarg, [txtopt, {'String', tedat{2*n-1}}]);
    end

    % Make the paired edit box
    edt(n) = makePairedEdit( txt(n), [editopt, {'Userdata', tedat{2*n}}] );

    % Set units to normalized
    txt(n).Units = 'normalized'; edt(n).Units = 'normalized';

end

end

function [txt] = makePairedText(anchor, marg, opt)
% Places and sizes initial text box label for edit box to go underneath
fig = ancestor(anchor,'figure');

% Get max TightInset for this row so all txtedt pairs are aligned
if isgraphics(anchor,'axes')
    figax = findobj(fig.Children, 'Type', 'axes');

    ax_pos = getAsMat(figax, 'Position');
    ax_sub = getSubPlotInd(ax_pos); % subplot-like indices for each axes
    ax_row = cellfun(@max, ax_sub(:,1)); % for our purposes we only want the lowest row each plot is part of

    rowax = figax( ax_row == ax_row(figax==anchor) );
    tins = max(getAsMat(rowax, 'TightInset'), [], 1);
end

% Make text object
txtset = [{fig, 'Style', 'text', 'Units', 'pixels'}, opt];
txt = uicontrol(txtset{:});

% Set label size and position relative to associated axes and previous
% labels

chkExtent(txt); % ensure string doesn't leave box
txt.Position(4) = txt.Extent(4); % remove vertical padding

txt.Units = 'pixels'; 
tpos = txt.Position;

% If an axes is passed as the anchor set below left corner, if a txt object
% set beside.
unit = anchor.Units; anchor.Units = 'pixels';
if isgraphics(anchor,'axes')
    vert = anchor.Position(2) - tins(2) - marg(4) - tpos(4); % set pair below axes
    horz = anchor.Position(1) + marg(1); % set pair on left plot edge
else
    vert = anchor.Position(2); % set pair at same height as last pair
    horz = sum(anchor.Position([1 3])) + sum(marg([1 3])); % set pair next to last pair
end
anchor.Units = unit;

% Finally set the position
txt.Position = [horz vert tpos(3:4)];

end

function [edt] = makePairedEdit(txt, edtopt)
% Places an edit box underneath a text field
fig = ancestor(txt,'figure');

% Make edit box
editset = [{fig, 'Style', 'edit', 'Units', 'pixels'}, edtopt];
edt = uicontrol(editset{:});

% Set position of edit box relative to label
edt.Units = txt.Units;
tpos = txt.Position;
epos = edt.Position;
dif3 = tpos(3) - epos(3); % to keep label and box centred
edt.Position = [tpos(1)+dif3/2 tpos(2)-epos(4) epos(3:4)];

end

function [ax_sub_new] = assignCBSubPlotInd(children, ax_sub, cbmap)
% Assign colorbars the same row and column as their associated axes.

n = 0; % offset
ax_sub_new = cell(length(children),2);
for c = 1:length(children)
    if cbmap(c)
        % Hidden colorbar property 'Axes' used
        ax_sub_new(c,:) = ax_sub( children(c).Axes==children(~cbmap),: );
        n = n + 1;
    else
        ax_sub_new(c,:) = ax_sub(c-n,:);
    end
end

end

function [linkedtips] = rmMovedLinkedTips(dcm, linkedtips, curtip)
% Checks if current tip was moved and removes its linked tips if so

cind = cellfun(@(C) any(ismember(C, curtip)), linkedtips); % cell index of tips linked to curtip
% Is DataIndex equal for the tips?
if any(cind)
    DataIndices = arrayfun(@(A) A.DataIndex, [linkedtips{cind}.Cursor]);
    if ~all( DataIndices == DataIndices(1) )
        rtind = linkedtips{cind} ~= curtip;
        rmCursorsNoLink([], [], dcm, linkedtips{cind}(rtind)); % remove tips linked to curtip
        linkedtips = linkedtips(~cind); % update linkedtips
    end
end

end

function [] = lnChoosePnt(src, event, n, slct, ui, lnksel)
% Assigned as ButtonDownFcn to a line to select nearest point when line is
% clicked
%
% Input
%   n: the index of the selected object in slct
%   slct: the select struct defined in the main function. This struct will
%         hold the index of the selected points
%   ui: the ui struct defined in the main function
%   lnksel: the option to link selections on all selectable lines

ax = ancestor(slct(n).xplobj, 'axes'); % Get the parent axes
fig = ancestor(ax, 'figure'); % Get the parent figure

% Check if the axes holds 3D information
is3D = false;
if ~isempty(slct(n).xplobj.ZData) && ~strcmp(slct(n).xplobj.Type, 'contour')
    is3D = true;
end

if strcmp(fig.SelectionType,'normal') % left click

    % Get plot-frame intersects of mouse
    pos = ax.CurrentPoint;

    % Normalize by figure limits before search
    xn = abs( diff(ax.XLim) );
    x = slct(n).xplobj.XData/xn; 

    yn = abs( diff(ax.YLim) );
    y = slct(n).xplobj.YData/yn;

    zn = abs( diff(ax.ZLim) );
    z = slct(n).xplobj.ZData/zn; 

    pos = pos./[xn yn zn; xn yn zn];
    posln = (pos(2,:)-pos(1,:))./norm(pos(2,:)-pos(1,:)); % line through both intersects

    % 'contour' and 'surface' plots may have gridded data
    if any( strcmp(slct(n).xplobj.Type, {'contour','surface'}) )
        % In this case 'z' will always be gridded, 'x' and 'y' may not
        if isvector(x)
            [x, y] = meshgrid(x, y);
        end
    end

    x = x(:);
    y = y(:);
    z = z(:);

    % Find nearest data point - certain plot types are treated differently
    if ~is3D
        % 2D: just find the point closest to the first intersection
        ind = dsearchn( [ x, y ], pos(1,1:2) );

    else
        % 3D: find the point closest to the intersecting line
        tol = 1e-3; % distance tolerance
        nz = length(z);
        %     If pos1 = pnt on line, posln = direction vector of line, xyz = point to find distance for
        % dst = ||(pos1 - xyz) - ( (pos1 - xyz).posln )*posln||
        xyz_d = sqrt(sum( ((pos(1,:) - [x, y, z]) - sum( (pos(1,:) - [x, y, z]).*repmat(posln,nz,1), 2 ).*repmat(posln,nz,1)).^2, 2));
        ind = 1:nz;
        ind = ind(ismembertol(xyz_d, min(xyz_d), tol));
        % Of points that fit tol, choose the point closest to the screen
        [~, ii] = min(norm(pos(1,:)' - [x(ind), y(ind), z(ind)]'));
        ind = ind(ii);

    end

    % Update all points when link all selections is on, otherwise just the
    % selected point.
    if lnksel
        slist = 1:length(slct);
    else
        slist = n;
    end

    % Updating points indicated above
    for s = slist
        % Set selected index
        slct(s).ind = ind;

        % Snap to nearest data point
        if any( strcmp(slct(s).xplobj.Type, {'contour','surface'}) )
            if isvector(slct(s).xplobj.XData)
                [x,y] = meshgrid(slct(s).xplobj.XData, slct(s).xplobj.YData);
            else
                x = slct(s).xplobj.XData;
                y = slct(s).xplobj.YData;
            end
            p(1) = x(ind);
            p(2) = y(ind);
        else
            p(1) = slct(s).xplobj.XData(ind);
            p(2) = slct(s).xplobj.YData(ind);
        end

        if is3D
            p(3) = slct(s).xplobj.ZData(ind);
        else
            p(3) = 0;
        end

        % Update the graphical line object for the selected point
        if ~isempty(slct(s).pnt)
            delete(slct(s).pnt)
        end
        sax = ancestor(slct(s).xplobj, 'axes');
        slct(s).pnt = line(sax, p(1), p(2), p(3), 'DisplayName', 'Selection', 'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'm'); % TODO: allow specifying point options

        % Set edit box values
        for edt = ui.ax(s).edt
            edt.String = num2str( edt.UserData(ind), 5 );
        end
    end

end

% Update structure by updating the ButtonDownFcn arguments
for s = 1:length(slct)
    slct(s).xplobj.ButtonDownFcn(2:end) = {s, slct, ui, lnksel};
    if ~isempty(slct(s).pnt)
        slct(s).pnt.ButtonDownFcn = {@ lnChoosePnt, s, slct, ui, lnksel}; % let user click on the point too
    end
end

% Update slct structure in 'View Details' callback function
for pb  = ui.fig.pbtn'
    pb.Callback{3} = slct;
end

end
