function [ fig ] = exploreResults( fig, pbtnfcn, varargin )
% exploreResults launches an interactive plot from which detailed results 
% are passed to user defined functions from selected data points.
% By: Kaden Agrey
% v2.0 2019.07.12
% 
% Input
% fig: 
%   Figure handle on which to place ui elements. The script will 
%   automatically add space for the for the ui elements.
% 
% 
% pbtnfcn: 
%   This argument is used to initialize the push buttons. It should be a 
%   cell array with the button name as the first column, a function handle 
%   as the second column and arguments to the function as later columns. 
%   Each row corresponds to a different pushbutton, which will be added to 
%   the figure as needed.
% 
%   The first 4 arguments to the function MUST be reserved for
%   exploreResults. The latter arguments should correspond to those in the
%   cell array. Eg:
%       pbtnfcn = {'My Button', fcn, arg1, arg2};
%   will correspond to the user function,
%       function argout = myFunc(src, event, ui, slct, arg1, arg2, ... )
% 
%   The reserved arguments are:
%       src: matlab variable pointing to the src of the callback, which is
%       the pushbutton object.
%       event: matlab variable giving event information.
%       ui: exploreResults variable with handles for all ui elements and
%       explorable chart objects.
%           ui.dcm: data cursor manager object
%           ui.pbtn: push button object
%           ui.xpl(i).chart: explorable chart object (a line or contour)
%           ui.xpl(i).data: cell array containing axes and/or user data
%           associated with the chart object.
%       slct: exploreResults variable with information on the selected
%       points. 
%           slct(j).chart: chart object point is on
%           slct(j).chartnum: index of chart object in ui.xpl(i).chart
%           slct(j).links: if SelectionLinkCharts is true, index of other 
%           selected points, slct(j), linked to this one.
%           slct(j).index: linear index of selected point in associated 
%           chart data, ui.xpl(chartnum).data.
%           slct(j).point: value of associated data at index.
% 
% 
% Optional
% charts: 
%   A list of the graphics objects that will be selectable. This is not a 
%   list of axes, but of the axes children that will be selectable. 
%   Supported types are: lines (2D, 3D), scatters (2D, 3D), contours and
%   surfaces. Other types may work but have not been tested. These object 
%   options are set correctly by default in MATLAB, ensure you haven't 
%   changed them.
%       chart.HitTest = 'on'
%       chart.PickableParts = 'visible'
% 
% TODO: Find out how to force Matlab to wait for the renderer before
% executing the mainplot formatting. - drawnow; doesn't seem to work
% 
% 
% Name/Value Pairs
% DataFromAxes:
%   Boolean, true allows the function to automatically get the data from 
%   axes and chart objects.
% 
% 
% DataFromUser: 
%   A cell array with name/data pairs that will be displayed with the 
%   datatips. The cell array should be nested so that the first index 
%   corresponds to the chart list, and the next index to name/data pairs.
%       data = { {name1, data1; name2, data2}, ... % chart 1
%                {name1, data1}, ...               % chart 2
%                {name1, data1;, name2, data2} }   % chart 3
%   Data arrays must have the same number of elements and dimension as the 
%   chart from which it will be selected.
% 
% 
% SelectionLinkCharts: 
%   Boolean, if true will force all explorable charts to have the same
%   indices selected. Selecting a point on one chart will select the same 
%   point on all charts. 
% 
% TODO: Let user specify which plots to link
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
%   2019.07.12: Switched data display to use MATLABs datatips

%% --- Parse Inputs --- %%
in = inputParser(); % initialize parser object

% Validation functions
chkCharts = @(charts) checkCharts(charts,fig);

% Required
in.addRequired('fig', @isgraphics); % figure handle to build ui on
in.addRequired('pbtnfcn'); % function handle (with arguments) to call when button is pressed
% Optional Positional
in.addOptional('charts', getAllExplorableCharts(fig), chkCharts); % charts to select data points from
% Optional Name/Value Pairs
in.addParameter('DataFromAxes',true); % display data from axes in boxes
in.addParameter('DataFromUser',{}); % display boxes will be added with user data
in.addParameter('SelectionLinkCharts',true); % link selected points accross all selectable objects
% in.addParameter('SelectionPerChart', inf); % number of data points that can be selected per chart
% in.addParameter('SelectionProperties', []); % properties of selection marker

in.parse(fig, pbtnfcn, varargin{:})

% Pull fig out of the object for easier access
fig = in.Results.fig;

%% --- Initialize --- %%
% --- ui --- %
% Initialize the struct to reference ui objects, this is passed to 
% callbacks.
% Calling datacursormode here creates the mode object, whose default
% functions we will change to add our own functionality.
ui = struct('pbtn', [], 'dcm', datacursormode(fig), 'xpl', struct('chart', [], 'data', []));

% --- ui.xpl --- %
% Get chart objects and data that selected points will be associated with. 
% This is passed to callbacks.
for p = 1:length(in.Results.charts)
    ui.xpl(p).chart = in.Results.charts(p); % chart object handle

    % Get data from axes
    if in.Results.DataFromAxes
        % Get labels. If a label is empty use x,y,z as default.
        lab = {'x', 'y', 'z'};
        if ~isempty( ui.xpl(p).chart.Parent.XLabel.String )
            lab{1} = ui.xpl(p).chart.Parent.XLabel.String;
        end
        if ~isempty( ui.xpl(p).chart.Parent.YLabel.String )
            lab{2} = ui.xpl(p).chart.Parent.YLabel.String;
        end
        if ~isempty( ui.xpl(p).chart.Parent.ZLabel.String )
            lab{3} = ui.xpl(p).chart.Parent.ZLabel.String;
        end

        % Construct data cell
        if isempty( ui.xpl(p).chart.ZData )
            ui.xpl(p).data = {lab{1}, ui.xpl(p).chart.XData; ...
                              lab{2}, ui.xpl(p).chart.YData};
        else
            ui.xpl(p).data = {lab{1}, ui.xpl(p).chart.XData; ...
                              lab{2}, ui.xpl(p).chart.YData; ...
                              lab{3}, ui.xpl(p).chart.ZData};
        end
    end

    % Combine with user data if provided
    if ~isempty(in.Results.DataFromUser)
        ui.xpl(p).data = [ui.xpl(p).data; in.Results.DataFromUser{p}];
    end
end

% --- Window size and space --- %
% These variables define the size and spacing of the ui objects. The
% pushbuttons are placed at the bottom of the figure, so the figmargins are
% used and pushbutton bottom margin is set to 0.
pbsize = [100 25]; % [pixels] | [width height] PushButton size
pbmargins = [3 0 0 3]; % [pixels] | [left bottom right top] PushButton margins
figmargins = [0 10 0 10]; % [pixels] | [left bottom right top] Figure window inside margins

% --- Collection options to pass to callbacks --- %
% Fields to pass to callbacks as options
opt_fields = {'SelectionLinkCharts'};
% opt_fields = {'SelectionLinkCharts','SelectionPerChart'};

% Define struct containing options to pass to the callbacks
opt = struct();
for field = opt_fields
    opt.(field{:}) = in.Results.(field{:});
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
        if ~ismember(ax.Children(ch), in.Results.charts)
            ax.Children(ch).PickableParts = 'none';
        end
    end
end

%% --- Setup UI --- %%
% --- Make pushbuttons specified by user input --- %
horz = getLeftChild(fig, 'pixels', 'Position'); % place in line with farthest left plot
ui.pbtn = gobjects(1, size(in.Results.pbtnfcn,1));
for pb = 1:size(in.Results.pbtnfcn,1)
    % Pushbutton position is in line with leftmost figure or next to the 
    % last pushbutton with a size defined in pbsize.
    pbpos = [horz + (pbsize(1) + pbmargins(1))*(pb-1), ...
             pbmargins(2) + figmargins(2), pbsize]; % [x y w h];
    % Make the button
    ui.pbtn(pb) = uicontrol(fig, 'Style', 'pushbutton',...
                                 'String', in.Results.pbtnfcn{pb,1},...
                                 'Units', 'pixels',...
                                 'Position', pbpos, ...
                                 'Callback', {@pbtnCallback, in.Results.pbtnfcn(pb,2:end), ui, {}, opt});

    ui.pbtn(pb).Units = 'normalized'; % set units
end

% --- Assign new callbacks to dcmode --- %
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

% --- Assign UpdateFcn to data cursor mode object --- %
ui.dcm.UpdateFcn = {@dcmUpdate, ui};

%% --- Finalize Figure Properties --- %%
ui.dcm.Enable = 'on'; % turn data cursor mode on

% Reset units on figure and all figure children
fig.Units = 'normalized';
set(fig.Children, 'Units', 'normalized');

end

%% --- Functions --- %%
% --- Getters --- %
function [charts] = getAllExplorableCharts(fig)
% Returns all axes children of explorable types

ch = findobj(fig.Children, 'flat', 'Type', 'axes', '-or', ...
                                   'Type', 'polaraxes', '-or', ...
                                   'Type', 'geoaxes');
charts = findobj([ch.Children], 'flat', '-not', {'Type', 'text', '-or', 'Type', 'light'});

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
function [] = setLinkedTips(fig, pbtn, lnkdt)
% Update <lnkdt> argument in figure and mode callbacks

% --- Update figure callback arguments --- %
disableFigFcnListener(fig)
fig.WindowButtonDownFcn{end}{4} = lnkdt;
fig.WindowButtonUpFcn{end}{4} = lnkdt;
fig.KeyPressFcn{end}{4} = lnkdt;

% --- Update uimode callback arguments --- %
dcmode = getuimode(fig, 'Exploration.Datacursor');
dcmode.WindowButtonDownFcn{4} = lnkdt;
dcmode.WindowButtonUpFcn{4} = lnkdt;
dcmode.KeyPressFcn{4} = lnkdt;

% --- Update pushbutton callback arguments --- %
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

% Sometimes data2fig selects the wrong point. This can be due to rounding
% in createDatatip, or data2fig not accounting for all
% transformations for 3D plots (a surprisingly difficult topic). Too fix 
% this we need to correct the point using increment functions. This should 
% never need to increment very many steps, unless the grid is extremely 
% fine compared to the figure size.
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
function [] = checkCharts(charts, fig)
% Checks that all <charts> belong to <fig>.
for c = charts(:)'
    if ~isequal(ancestor(c, 'figure'), fig)
        error(['The gobject ' c.DisplayName 'is not part of the ' ...
               'figure. All exploreable charts must be on the same figure.'])
    end
end

end

function [] = disableFigFcnListener(fig)
% This uses undocumented functionality, see link below if it breaks. We
% need to disable some listeners so that we can change the button down
% function of the figure and interject with our own code.
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
isvec = cellfun(@(C) isvector(C), data);
if any(isvec)
% The spatial data is in vector form. Note that this assumes the user
% provides only data with the same dimensions as the level data, which 
% they should.
    % Get subscript indices of point to return.
    levi = find(~isvec,1); % index of 'level' data
    sub = zeros(1,length(data));
    if levi < 4
        [sub(2), sub(1)] = ind2sub(size(data{levi}), ind);
    else
        [sub(2), sub(1), sub(3)] = ind2sub(size(data{levi}), ind);
    end
    sub(levi:end) = ind; % the level data is gridded, use linear index

else
% The spatial data is in matrix form. Note that this assumes the user
% provides only data with the same dimensions as the level data, which 
% they should.
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
    else % XData is a nd grid ERROR: This won't work for 3D data ----------------------------------------------------------------------
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
% Converts data coordinates in <pnt> on an axes to equivalent
% coordinates in the figure window with the same units as <ax>. This code 
% is influenced by (but not copied from) the excellent script ds2fig() by 
% MinLong Kwong on file exchange. It gets a reasonable approximation, but
% isn't perfect. The results are corrected later.
%
% Using third party functions - as a comparison for bug fixing
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
dScl = [ax.DataAspectRatio, 1]';

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
% 
% dcm:
%   Datacursormanager object.
% 
% dc:
%   Array of data cursors to be removed.

for c = dc(:)'
    dcm.removeDataCursor(c)
end

end

function [] = rmCursorsNoLink(~, ~, dcm, dt)
% Invokes the datacursormanager.removeCursor function to remove an array of
% datatips after disabling the delete callback function
%
% Input
% ~: 
%   First two arguments are placeholders for when this function is set as a
%   callback.
% 
% dcm:
%   Datacursormanager object.
% 
% dt:
%   Array of datatips to be removed.

set(dt, 'DeleteFcn', []);
rmCursors([], [], dcm, [dt.Cursor])

end

function [] = mvCursors(~,~, dc, cur)
% Invokes the datacursor.moveTo function to move an array of cursors to the
% same index on their sources as the cursor <cur>
%
% Input
% ~: 
%   First two arguments are placeholders for when this function is set as a
%   callback.
% 
% dc:
%   Array of data cursors to be moved.
% 
% cur:
%   Cursor whose data index to match.

ind = cur.DataIndex;
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

    figpnt = data2fig(c.DataSource.Parent, pnt);
    c.moveTo(figpnt);

    c.DataSource.Parent.Units = units;

    % Sometimes data2fig selects the wrong point. This can be due to rounding
    % in createDatatip, or data2fig not accounting for all
    % transformations for 3D plots (a surprisingly difficult topic). Too fix 
    % this we need to correct the point using increment functions. This should 
    % never need to increment very many steps, unless the grid is extremely 
    % fine compared to the figure size.
    incrementCursorToIndex(c, ind);
end

end

function [] = mvLinkedTips(linkedtips, curtip)
% Moves all tips linked to the current datatip to the same index on their
% source.
% 
% Input
% linkedtips: 
%   Cell array where each element is a list of linked datatips.
% 
% curtip:
%   Current tip in the figure window.

cind = cellfun(@(C) any(ismember(C, curtip)), linkedtips); % cell index of tips linked to curtip

if any(cind)
    % Is DataIndex equal for the tips?
    DataIndices = arrayfun(@(A) A.DataIndex, [linkedtips{cind}.Cursor]);
    if ~all( DataIndices == DataIndices(1) )
        rtind = linkedtips{cind} ~= curtip; % remaining tips
        mvCursors([], [], [linkedtips{cind}(rtind).Cursor], curtip.Cursor); % move tips linked to curtip
    end
end

end

function [] = mvSrcLinkedTips(linkedtips, curtip, xplcharts)
% If the source of the current cursor has changed then its linked tips are 
% moved to the remaining sources. This algorithm is fast enough for
% reassigning a few tips, but can be made much more efficient if required.
% 
% Input
% linkedtips: 
%   Cell array where each element is a list of linked datatips.
% 
% curtip:
%   Current tip in the figure window.
% 
% xplcharts:
%   List of selectable charts

cind = cellfun(@(C) any(ismember(C, curtip)), linkedtips); % cell index of tips linked to curtip

if any(cind)
    % Are all explorable sources covered?
    sources = arrayfun(@(A) A.DataSource, [linkedtips{cind}.Cursor], 'UniformOutput', false); % sources of linked tips
    cursrc = curtip.Cursor.DataSource; % source of current tip

    % Find which xplcharts have been used
    xplused = false(length(xplcharts),1);
    for m = 1:length(xplcharts)
        X = xplcharts(m);

        for S = sources
            if X == S{1}
                xplused(m) = true;
                break;
            end
        end
    end

    curxpl = xplcharts == cursrc; % current xplchart
    if ~all(xplused) % if all xplcharts aren't used
        xplused = curxpl; % reset xplused
        curlnk = linkedtips{cind} == curtip; % find current tip in linkedtips
        % Assign cursors to remaining data sources
        for tip = linkedtips{cind}(~curlnk)
            x = find(~xplused,1); % location of first unused source
            tip.DataSource = xplcharts(x);
            xplused(x) = true;
        end
    end

end

end

function [] = incrementCursorToIndex(cursor, ind)
% Increments data cursor to <ind>
% 
% Input
% cursor:
%   Cursor object to increment.
% 
% ind:
%   Index to increment to.

% If indices are already the same, return
if cursor.DataIndex == ind
    return;
end

source = cursor.DataSource;

% Get the x and y indices
if any(strcmp(source.Type, {'contour', 'surface'}))
    [ytrue, xtrue] = ind2sub(size(source.ZData),ind);
    [y, x] = ind2sub(size(source.ZData),cursor.DataIndex);
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
function str = dcmUpdate(pdtobj, eobj, ui)
% Builds string for data cursor display.

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

% If an explorable object was hit we continue with "normal" selection type 
% or "extend" if the modifier is "shift" or "alt", as this means a new data 
% cursor was created.
isxplhit = any(s);
% Adapted from %matlabroot%/toolbox/matlab/graphics/
%               datacursormanager.m@localWindowButtonDownFcn()
mod = get(fig,'CurrentModifier');
isAddSelType = strcmp(fig.SelectionType, 'normal');  % selection type is 'normal'?
isAddMod = numel(mod)==1  && strcmp(fig.SelectionType, 'extend') ...
    && any( strcmp(mod{1}, {'shift','alt'})  ); % selection type is 'extend' with 'shift' or 'alt'?
if ~isxplhit || ~( isAddSelType || isAddMod )
    return; % return if no tip was made
end

% Store the current cursor to reset after
curcur = ui.dcm.CurrentCursor;
alldt = findall(fig.Children, 'Type', 'hggroup');

% Even if the correct key combinations were pressed it may not have 
% resulted in a new data tip. Determine if a new cursor was added.
isAdded = false;
if opt.SelectionLinkAxes
    if numel(alldt) < numel(ui.xpl) || ( numel([lnkdt{:}]) < numel(alldt) )
        isAdded = true;
    end
end

% Link Axes
if opt.SelectionLinkAxes && isAdded
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
    % NOT YET REQUIRED
end

end

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

    % DataIndex .(index) from cursor info if available. Else find it from 
    % cursor info .(Position).
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

function [] = mvLinkedTipsButtonUp(src, event, fcn, dcm, lnkdt)
% If click and drag is used to move a tip then we need to update tip
% locations on the button up action.

% Run the original callback
if isa(fcn,'function_handle')
   fcn(src, event);
end

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
% If arrow keys are used to move a tip then we need to update linked tip 
% locations the same way.

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

function [] = rmLinkedCursors(srcdt, ~, fig, dcm, pbtns, dt)
% Delete function callback when cursors are linked. Removes linked tips and
% updates <lnkdt> callback argument for callback arguments.

    rmCursorsNoLink([], [], dcm, dt)

    % To avoid needing to update every delete callback for every datatip we
    % will get the lnkdt variable from the mode callback arguments.
    lnkdt = getLinkedTips(fig);

    % Remove deleted tips from lnkdt
    cind = cellfun(@(C) any(ismember(C, srcdt)), lnkdt);
    lnkdt = lnkdt(~cind); 

    setLinkedTips(fig, pbtns, lnkdt)
end
