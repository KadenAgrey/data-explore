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

% Validation functions
checkCharts = @(A) chkCharts(A,fig);

% Required
in.addRequired('fig', @isgraphics); % figure handle to build ui on
in.addRequired('pbtnfcn'); % function handle (with arguments) to call when button is pressed
% Optional Positional
in.addOptional('lines', getAllExplorableLines(fig), checkCharts); % charts to select data points from
% Optional Name/Value Pairs
in.addParameter('DataFromAxes',true); % display data from axes in boxes
in.addParameter('DataFromUser',{}); % display boxes will be added with user data
in.addParameter('SelectionLinkAxes',true); % link selected points accross all selectable objects
% in.addParameter('SelectionPerChart', inf); % number of data points that can be selected per axes
% in.addParameter('SelectionPerAxes', inf); % number of data points that can be selected per axes
% in.addParameter('SelectionProperties', []); % properties of selection marker

in.parse(fig, pbtnfcn, varargin{:})

% Pull fig out of the object for easier access
fig = in.Results.fig;

%% --- Initialize --- %%
% These variables define the size and spacing of the ui objects. The
% pushbuttons are placed at the bottom of the figure, so the figmargins are
% used and pushbutton bottom margin is set to 0.
pbsize = [100 25]; % [pixels] | [width height] PushButton size
pbmargins = [3 0 0 3]; % [pixels] | [left bottom right top] PushButton margins
figmargins = [0 10 0 10]; % [pixels] | [left bottom right top] Figure window inside margins

% Fields to pass to callbacks as options
opt_fields = {'SelectionLinkAxes'};
% opt_fields = {'SelectionLinkAxes','SelectionPerChart',...
%               'SelectionPerAxes'};

% Initialize the struct to reference ui objects, this is passed to 
% callbacks.
%   ui.pbtn is an array holding all pbtn objects
%   ui.dcm holds the datacursormode object
% Calling datacursormode here creates the mode object, whose default
% functions we will change to add our own functionality.
ui = struct('pbtn', [], 'dcm', datacursormode(fig), 'xpl', struct('chart', [], 'data', []));

% Get objects selected points will be associated with. This is passed to
% callbacks.
if isempty(in.Results.DataFromUser)
    user_data = cell(length(in.Results.lines),1);
else
    user_data = in.Results.DataFromUser;
end

for p = 1:length(in.Results.lines)
    ui.xpl(p).chart = in.Results.lines(p);

    % Get labels. If a label is empty use x,y,z as default.
    lab = {'x', 'y', 'z'};
    if ~isempty(ui.xpl(p).chart.Parent.XLabel.String)
        lab{1} = ui.xpl(p).chart.Parent.XLabel.String;
    end
    if ~isempty(ui.xpl(p).chart.Parent.YLabel.String)
        lab{2} = ui.xpl(p).chart.Parent.YLabel.String;
    end
    if ~isempty(ui.xpl(p).chart.Parent.ZLabel.String)
        lab{3} = ui.xpl(p).chart.Parent.ZLabel.String;
    end

    % Construct data cell
    if isempty(ui.xpl(p).chart.ZData)
        ui.xpl(p).data = {lab{1}, ui.xpl(p).chart.XData; ...
                        lab{2}, ui.xpl(p).chart.YData};
    else
        ui.xpl(p).data = {lab{1}, ui.xpl(p).chart.XData; ...
                        lab{2}, ui.xpl(p).chart.YData; ...
                        lab{3}, ui.xpl(p).chart.ZData};
    end
    ui.xpl(p).data = [ui.xpl(p).data; user_data{p}];
end

% Define a struct containing selection options to pass to the callbacks
opt = struct();
for f = opt_fields
    opt.(f{:}) = in.Results.(f{:});
end

%% --- Prepair Main Figure --- %%
% --- Figure size and plot positioning --- %
fig.Units = 'pixels';
set(fig.Children, 'Units', 'pixels');

% Move the children
children = findobj(fig.Children, 'flat', '-not', 'AxisLocationMode', 'auto');
y = getBottomChild(fig, 'pixels', 'OuterPosition');
for ch = children'
    ch.Position(2) = ch.Position(2)-y + pbsize(2) + pbmargins(4) + figmargins(2);
end
% Resize the window to tightly wrap the children on the top and bottom
y = getTopChild(fig, 'pixels', 'OuterPosition');
fig.Position(4) = y + figmargins(4);

% --- Make non-explorable charts 'unpickable' --- %
% Set the 'PickableParts' property for non-explorable axes children to
% 'none'.
axs = findobj(fig.Children, 'flat', 'Type', 'axes', '-or', ...
                                   'Type', 'polaraxes', '-or', ...
                                   'Type', 'geoaxes');
for ax = axs'
    for ch = 1:length(ax.Children)
        if ~ismember(ax.Children(ch), in.Results.lines)
            ax.Children(ch).PickableParts = 'none';
        end
    end
end

%% --- Setup UI --- %%
% Make button to view details of selected point
horz = getLeftChild(fig, 'pixels', 'Position'); % place in line with farthest left plot
ui.pbtn = gobjects(1, size(in.Results.pbtnfcn,1));
for pb = 1:size(in.Results.pbtnfcn,1)
    % Pushbutton position is in line with leftmost figure or next to the 
    % last pushbutton with a size defined in pbsize.
    pbpos = [horz + (pbsize(1) + pbmargins(1))*(pb-1), ...
             pbmargins(2) + figmargins(2), pbsize];
    % Make the button
    ui.pbtn(pb) = uicontrol(fig, 'Style', 'pushbutton',...
                                 'String', in.Results.pbtnfcn{pb,1},...
                                 'Units', 'pixels',...
                                 'Position', pbpos, ...
                                 'Callback', {@pbtnCallback, in.Results.pbtnfcn(pb,2:end), ui, {}, opt});

    ui.pbtn(pb).Units = 'normalized'; % set units
end

% Assign default callbacks to dcmode
dcmode = getuimode(fig, 'Exploration.Datacursor');
% Set button down function
fcn = {@lnSelectPnt, dcmode.WindowButtonDownFcn, ui, {}, opt};
dcmode.WindowButtonDownFcn = fcn;
% Set button up function
fcn = {@mvLinkedTipsButtonUp, dcmode.WindowButtonUpFcn, ui.dcm, {}};
dcmode.WindowButtonUpFcn = fcn;
% Set key press function
fcn = {@mvLinkedTipsKeyPress, dcmode.KeyPressFcn, ui.dcm, {}};
dcmode.KeyPressFcn = fcn;

% Assign UpdateFcn to data cursor mode object
ui.dcm.UpdateFcn = {@dcmUpdate, ui};

%% --- Finalize Figure Properties --- %%
ui.dcm.Enable = 'on'; % turn data cursor mode on

% Reset units on figure and all figure children
fig.Units = 'normalized';
set(fig.Children, 'Units', 'normalized');

end

%% --- Functions --- %%
% --- Getters --- %
function [lines] = getAllExplorableLines(fig)
% If lns is empty returns all axes children of explorable types as the
% default behavior

ch = findobj(fig.Children, 'flat', 'Type', 'axes', '-or', ...
                                   'Type', 'polaraxes', '-or', ...
                                   'Type', 'geoaxes');
lines = findobj([ch.Children], 'flat', '-not', {'Type', 'text', '-or', 'Type', 'light'});

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

function [pos, ind] = getLeftChild(parent, unit, prop)
% Gets leftmost child object
children = [parent.Children]; % all children

oldunit = get(children, 'Units'); % to reset units
set(children, 'Units', unit); % set units

position = [children.(prop)]; % all positions
[pos,ind] = min(position(1:4:end)); % minimum x-position and child index

set(children, {'Units'}, oldunit); % reset units

end

function [pos, ind] = getBottomChild(parent, unit, prop)
% Gets lowest child object
children = [parent.Children]; % all children

oldunit = get(children, 'Units'); % to reset units
set(children, 'Units', unit); % set units

position = [children.(prop)]; % all positions
[pos,ind] = min(position(2:4:end)); % minimum x-position and child index

set(children, {'Units'}, oldunit); % reset units

end

function [pos, ind] = getTopChild(parent, unit, prop)
% Gets lowest child object
children = [parent.Children]; % all children

oldunit = get(children, 'Units'); % to reset units
set(children, 'Units', unit); % set units

position = [children.(prop)]; % all positions
[pos,ind] = max(position(2:4:end) + position(4:4:end)); % maximum y-position and child index

set(children, {'Units'}, oldunit); % reset units
end

function [lnkdt] = getLinkedTips(fig)
% Gets linked tips cell array from fig window button down function
% arguments.

dcmode = getuimode(fig, 'Exploration.Datacursor');
if ~isempty(dcmode)
    lnkdt = dcmode.WindowButtonDownFcn{4};
else % if the figure is being closed return an empty cell
    lnkdt = {};
end

end

% --- Setters --- %
function [] = setTipCallbacks(src, event, fig, ui, opt, lnkdt)
% Assign button down callbacks for figure window
disableFigFcnListener(fig)

% Set button down function
wndwfcn = {@lnSelectPnt, fig.WindowButtonDownFcn, lnkdt, ui, opt};
fig.WindowButtonDownFcn = wndwfcn;
% Set button up function
wndwfcn = {@mvLinkedTipsButtonUp, fig.WindowButtonUpFcn, ui.dcm, lnkdt};
fig.WindowButtonUpFcn = wndwfcn;
% Set key press function
wndwfcn = {@mvLinkedTipsKeyPress, fig.KeyPressFcn, ui.dcm, lnkdt};
fig.KeyPressFcn = wndwfcn;

end

function [] = setLinkedTips(fig, pbtn, lnkdt)
% Update lnkdt argument in figure and mode callbacks

setLinkedTipsFig(fig, lnkdt);
dcmode = getuimode(fig, 'Exploration.Datacursor');
setLinkedTipsMode(dcmode, lnkdt);
setLinkedTipsButton(pbtn, lnkdt);

end

function [] = setLinkedTipsFig(fig, lnkdt)

disableFigFcnListener(fig)
fig.WindowButtonDownFcn{end}{4} = lnkdt;
fig.WindowButtonUpFcn{end}{4} = lnkdt;
fig.KeyPressFcn{end}{4} = lnkdt;

end

function [] = setLinkedTipsMode(mode, lnkdt)

mode.WindowButtonDownFcn{4} = lnkdt;
mode.WindowButtonUpFcn{4} = lnkdt;
mode.KeyPressFcn{4} = lnkdt;

end

function [] = setLinkedTipsButton(pbtn, lnkdt)

for p = pbtn
    if isgraphics(p)
        p.Callback{4} = lnkdt;
    end
end

end

% --- Makers --- %
function [dt] = makeDataCursor(dcm, target, ind, properties)

% Get the cursor position in data units
if isempty(target.ZData)
    % Reconsider doing things this way, it requires constructing a cell 
    % array out of the plot data. Alternative: write a function to do the
    % same thing as ind2pnt but using target data instead of a cell array.
    pnt = ind2pnt(target, {target.XData, target.YData}, ind);
else
    pnt = ind2pnt(target, {target.XData, target.YData, target.ZData}, ind);
end

% Get cursor position in pixels
units = get( target.Parent, 'Units' ); % store fig units
set( target.Parent, 'Units', 'pixels' ); % set to pixels
figpnt = data2fig(target.Parent, pnt); % get position
set( target.Parent, 'Units', units ); % % reset units

dt = dcm.createDatatip(target, figpnt);

% Sometimes the wrong point is selected (probably due rounding in 
% createDatatip), so we need to correct this using increment functions. 
% This should never need to increment very many steps, unless the grid is 
% extremely fine compared to the figure size.
incrementCursorToIndex(dt.Cursor, ind);

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
function [] = chkCharts(charts, fig)
% Checks that all charts belong to fig.
for c = charts(:)'
    if ~isequal(ancestor(c, 'figure'), fig)
        error(['The gobject ' c.DisplayName 'is not part of the ' ...
               'figure. All exploreable charts must be on the same figure.'])
    end
end

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
% accounts for gridded data. If the data are from a meshgrid object
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
        [sub(2), sub(1)] = ind2sub(size(data{levi}), ind);
    else
        [sub(2), sub(1), sub(3)] = ind2sub(size(data{levi}), ind);
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

if any(strcmp(line.Type, {'contour', 'surface'}))
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
    ind = sub2ind(size(line.ZData),y,x);
else
% Data is a 1D sequence.
    ind = find(line.XData == pnt(1), 1);
end

end

function [pntf] = data2fig(ax, pnt)
% Converts data coordinates (<x> and <y>) on an axes to equivalent
% coordinates in the figure with the same units as <ax>.
%
% Using third party functions - as a comparison for bug fixing
%   This works, but not perfectly. 
% [pntf(1), pntf(2)] = ds2fig(ax, pnt(1), pnt(2), pnt(3));
% pntf = pntf.*ax.Parent.Position(3:4);
% return;

% If pnt is 2D we set the z coordinate to 0.
if length(pnt) == 2
    pnt = [pnt 0];
end

% Get bounding box of projection from axis limits. We will project this
% into our 2D data space to get the 2D axis limits.
box3D = [ax.XLim(1), ax.YLim(1), ax.ZLim(1), 1; ...
         ax.XLim(1), ax.YLim(1), ax.ZLim(2), 1; ...
         ax.XLim(1), ax.YLim(2), ax.ZLim(1), 1; ...
         ax.XLim(2), ax.YLim(1), ax.ZLim(1), 1; ...
         ax.XLim(2), ax.YLim(2), ax.ZLim(2), 1; ...
         ax.XLim(1), ax.YLim(2), ax.ZLim(2), 1; ...
         ax.XLim(2), ax.YLim(1), ax.ZLim(2), 1; ...
         ax.XLim(2), ax.YLim(2), ax.ZLim(1), 1;]';

% Data stretching is accounted for by DataAspectRatio.
% pScl = [ax.PlotBoxAspectRatio, 1]';
% pScl = [1 1 1 1]';
dScl = [ax.DataAspectRatio, 1]';
% dScl = [1 1 1 1]';

% Get view transformation matrix. This will project our 3D data into 2D
% space.
if strcmp(ax.Projection, 'orthographic')
% Orthographic project is used
    d2f = viewmtx(ax.View(1),ax.View(2));
else
% Perspective projection is used
    dnorm = [diff(ax.XLim), diff(ax.YLim), diff(ax.ZLim)];
    d2f = viewmtx(ax.View(1),ax.View(2),ax.CameraViewAngle,ax.CameraTarget./dnorm);
end

% Transform bounding box to 2D space
% We scale the 3D box by the PlotBoxAspectRatio then transform it to the 2D
% space. The same procedure will be followed for the pnt.
box3DTrans = d2f*( box3D./dScl );
box3DTrans(1:3,:) = box3DTrans(1:3,:)./box3DTrans(4,:); % scale by the homogenous vector for perspective projection
% Each column of box2D is a pair of 2D axis limits so that 
% box2D = [xlim_lo, ylim_lo; xlim_hi, ylim_hi]
box2D = [min(box3DTrans(1:2,:),[],2), max(box3DTrans(1:2,:),[],2)]';

% Get pnt in 2D space
% We add the 1 to account for the homogenous vector, then scale as with
% box3D and transform the point to 2D space. The third element of pnt2D is
% effectively a measure of "depth" into the screen.
pnt2D = d2f*( [pnt 1]'./dScl );

% Convert to figure space
pos = ax.Position;

% Account for the 2D box aspect ratio not filling the position rectangle.
ARmodes = {ax.PlotBoxAspectRatioMode, ax.DataAspectRatioMode};
boxAR = diff(box2D(:,1))/diff(box2D(:,2));
posAR = pos(3)/pos(4);
if any(strcmp(ARmodes, 'manual')) && posAR < boxAR
% Adjust y-pos
    d = pos(4) - pos(3)/boxAR;
    pos(2) = pos(2) + d/2;
    pos(4) = pos(4) - d;
elseif any(strcmp(ARmodes, 'manual')) && posAR > boxAR
% Adjust x-pos
    d = pos(3) - pos(4)*boxAR;
    pos(1) = pos(1) + d/2;
    pos(3) = pos(3) - d;
end

pntax = ( (pnt2D(1:2)'/pnt2D(4) - box2D(1,:))./diff(box2D) ); % position of point normalized to 2D plot box
pntf = pos(1:2) + pos(3:4).*pntax;

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
    if c.DataIndex == ind
    % If the index hasn't changed then we don't need to move this tip.
        continue;
    end

    if isempty(c.DataSource.ZData)
        % Reconsider doing things this way, it requires constructing a cell 
        % array out of the plot data. Alternative: write a function to do the
        % same thing as ind2pnt but using target data instead of a cell array.
        pnt = ind2pnt(c.DataSource, {c.DataSource.XData, c.DataSource.YData}, ind);
    else
        pnt = ind2pnt(c.DataSource, {c.DataSource.XData, c.DataSource.YData, c.DataSource.ZData}, ind);
    end

    units = c.DataSource.Parent.Units;
    c.DataSource.Parent.Units = 'pixels';
    unitsfig = get(ancestor(c.DataSource, 'figure'), 'Units');
    set(ancestor(c.DataSource, 'figure'), 'Units', 'pixels');

    figpnt = data2fig(c.DataSource.Parent, pnt);
    c.moveTo(figpnt);

    c.DataSource.Parent.Units = units;
    set(ancestor(c.DataSource, 'figure'), 'Units', unitsfig);

    % Sometimes the wrong point is selected (probably due rounding in 
    % createDatatip), so we need to correct this using increment functions. 
    % This should never need to increment very many steps, unless the grid is 
    % extremely fine compared to the figure size.
    incrementCursorToIndex(c, ind);
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
        X = xplobj(m);

        for S = Sources
            if X == S{1}
                xplUsed(m) = true;
                break;
            end
        end
    end

    curxpl = xplobj == cursrc; % current xplobj
    if ~all(xplUsed) % if all xplobj aren't used
        xplUsed = curxpl; % reset xplUsed
        curlnk = linkedtips{cind} == curtip; % find current tip in linkedtips
        % Assign cursors to remaining data sources
        for tip = linkedtips{cind}(~curlnk)
            x = find(~xplUsed,1); % location of first unused src
            tip.DataSource = xplobj(x);
            xplUsed(x) = true;
        end
    end

end

end

function [] = incrementCursorToIndex(cursor, ind)
% Increments data cursor to <ind>

% If indices are already the same, return
if cursor.DataIndex == ind
    return;
end

target = cursor.DataSource;

% Get the x and y indices
if any(strcmp(target.Type, {'contour', 'surface'}))
    [ytrue, xtrue] = ind2sub(size(target.ZData),ind);
    [y, x] = ind2sub(size(target.ZData),cursor.DataIndex);
else
    xtrue = ind; 
    ytrue = ind;
    y = cursor.DataIndex;
    x = ind; % it's not a grid, we only need to increment one dimension
end

% Get distance to increment
dx = xtrue - x;
dy = ytrue - y;

% Increment tip set distance
for inc = x:sign(dx):xtrue-sign(dx) % number of increments one less than abs(dx)
    if dx < 0
        direction = 'left';
    else
        direction = 'right';
    end
    cursor.increment(direction);
end
for inc = y:sign(dy):ytrue-sign(dy) % number of increments one less than abs(dy)
    if dy < 0
        direction = 'down';
    else
        direction = 'up';
    end
    cursor.increment(direction);
end

end

% --- Callbacks --- %
function [] = pbtnCallback(src, event, fcn, ui, lnkdt, opt)
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

% Contains information associated with selected points data.
slct = struct('chart', [], 'chartnum', [], 'links', [], 'index', [], 'point', []);

cinfo = ui.dcm.getCursorInfo;
for p = 1:length(cinfo)
    chartnum = find( [ui.xpl.chart]==cinfo(p).Target, 1 ); % index of line in list of explorable charts

    % Get Target .(chart)
    slct(p).chart = cinfo(p).Target;
    slct(p).chartnum = chartnum;

    % DataIndex .(index) and from cursor info if available. Else find it
    % from cursor info .(Position).
    if isfield(cinfo, 'DataIndex')
        slct(p).index = cinfo(p).DataIndex;
    else
        slct(p).index = pnt2ind(cinfo(p).Target, cinfo(p).Position);
    end

    % Get .(point) manually
    slct(p).point = ind2pnt(cinfo(p).Target, ui.xpl(chartnum).data(:,2), slct(p).index);
end

% Get links
if opt.SelectionLinkAxes
    % For each set of linked tips, get their indices in the slct struct
    for lnk = lnkdt
        ind = find([slct.index] == lnk{1}(1).Cursor.DataIndex);

        % For each ind, assign the others as links in each slct.links
        nn = 1:length(ind);
        for n = nn
            slct(ind(n)).links = ind(n~=nn);
        end
    end
end

% Call the user's function and pass arguments through
fcn{1}( src, event, ui, slct, fcn{2:end} );

end

function [] = lnSelectPnt(fig, event, fcn, ui, lnkdt, opt)
% Callback assigned to each selectable object (line, contour, surface
% etc.). This is interjected before the standard matlab datatip mode
% callback is executed and allows us to manage the datatips directly.

% Run the original callback
if isa(fcn,'function_handle')
   fcn(fig, event);
end

% Get index of selected object
% Note: undocumented event property "HitObject"
xlns = [ui.xpl.chart];
s = xlns == event.HitObject;

% If an explorable object was hit we need to continue with "normal" or 
% "extend" if the modifier is "shift" or "alt" as this means a new data 
% cursor was created.
% Note: undocumented event property "HitObject"
isxplhit = any(s);
% Copied from %matlabroot%/toolbox/matlab/graphics/
%               datacursormanager.m@localWindowButtonDownFcn()
mod = get(fig,'CurrentModifier');
% hMode = getuimode(hFig,'Exploration.Datacursor'); % don't need to check this if using change in alldt as the check

isAddSelType = strcmp(fig.SelectionType, 'normal'); 
isAddMod = numel(mod)==1  && strcmp(fig.SelectionType, 'extend') ...
    && any( strcmp(mod{1}, {'shift','alt'})  );
% isAddOption = hMode.ModeStateData.newCursor;
if ~isxplhit || ~( isAddSelType || isAddMod )
    return;
end

% Store the current cursor to reset after
curcur = ui.dcm.CurrentCursor;
alldt = findall(fig.Children, 'Type', 'hggroup');

% Determine if a new cursor was added
isAddRequest = false;
if opt.SelectionLinkAxes
    if numel(alldt) < numel(ui.xpl) || ( numel([lnkdt{:}]) < numel(alldt) )
        isAddRequest = true;
    end
end

% Link Axes
if opt.SelectionLinkAxes && isAddRequest
    % Get index of point
    ind = ui.dcm.CurrentCursor.DataIndex;

    % Make Linked data tips
    clns = [alldt.DataSource]; % cursor lines
    newdt = gobjects(0);
    for ln = xlns(~s)
        n = clns == ln;
        cindices = arrayfun(@(A) A.DataIndex, [alldt(n).Cursor]); % cursor indices
        % Are there no tips on ln or are the tips at different indices?
        if all(~n) || all(cindices ~= ind)
            newdt(length(newdt)+1) = makeDataCursor(ui.dcm, ln, ind, []);
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
        lnkdt{end}(t).DeleteFcn = {@rmLinkedCursors, fig, ui.dcm, ui.pbtn, lnkdt{end}(rt)};
        % Add property listener to move tips together - listeners are not
        % implimented for Datatips yet.
%         plistenP = addlistener(dt(t),'Position','PostSet',{@rmCursorsNoLink, ui.dcm, dt(rt)});
%         plistenS = addlistener(dt(t),'DataSource','PostSet',{@rmCursorsNoLink, ui.dcm, dt(rt)});
    end

    % Reset current cursor
    ui.dcm.CurrentCursor = curcur;

    % Update callback args
    setLinkedTips(fig, ui.pbtn, lnkdt)

elseif opt.SelectionLinkAxes
    % Ensure datatips move together
    curdt = findobj(alldt,'Cursor',curcur);
    mvSrcLinkedTips(lnkdt, curdt, xlns)
    % We don't need to move the tips here because this is handled in the
    % button up function. We only ensure data sources are correct.
end

% Remove extra datatips
if false
    % NOT YET IMPLIMENTED
end

end

function [] = rmLinkedCursors(srcdt, ~, fig, dcm, pbtns, dt)
% Delete function callback when cursors are linked. Removes linked tips and
% updates lnkdt callback argument for fig button down, button up, and key 
% press functions.

    rmCursorsNoLink([], [], dcm, dt)

    % To avoid needing to update every delete callback for every datatip we
    % will get the lnkdt variable from the mode callback arguments.
    lnkdt = getLinkedTips(fig);
    cind = cellfun(@(C) any(ismember(C, srcdt)), lnkdt);
    lnkdt = lnkdt(~cind);

    setLinkedTips(fig, pbtns, lnkdt)
end

function [] = mvLinkedTipsButtonUp(src, event, fcn, dcm, lnkdt)
% Callback to call mvLinkedTips with appropriate arguments.

% Run the original callback
if isa(fcn,'function_handle')
   fcn(src, event);
end

% % Execute the specified callback function
% hgfeval(newButtonUpFcn,hFig,evd);

% If selection type is normal continue.
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
direction = event.Key(1:end-length('arrow'));

% Get the current datatip
curdt = findobj([lnkdt{:}],'Cursor',dcm.CurrentCursor);

% Increment linked tips
cind = cellfun(@(C) any(ismember(C, curdt)), lnkdt); % cell index of tips linked to curtip
if any(cind)
    rtind = lnkdt{cind} ~= curdt;
    for c = [lnkdt{cind}(rtind).Cursor]
        c.increment(direction); % move cursor with curtip
    end
end

end

function str = dcmUpdate(pdtobj, eobj, ui)

% Get index of selected object
s = [ui.xpl.chart] == eobj.Target;
ind = pdtobj.Cursor.DataIndex;

% Get data point
pnt = ind2pnt(eobj.Target, ui.xpl(s).data(:,2), ind);

% Make string of data to display
str = cell(1,length(ui.xpl(s).data));
try % >= 20XXx
    ui.dcm.Interpreter;
    for d = 1:length(ui.xpl(s).data)
        str{d} = [ui.xpl(s).data{d,1} ' {\bf\color{DarkGreen}{' num2str(pnt(d),4) '}}'];
    end
catch % <= 20XXx
    pdtobj.TextColor = [0.3 0.6 0.3];
    pdtobj.FontWeight = 'bold';
    for d = 1:length(ui.xpl(s).data)
        str{d} = [ui.xpl(s).data{d,1} ' ' num2str(pnt(d),4)];
    end
end

end
