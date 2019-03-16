function [ fig ] = exploreResults( mainfig, pbtnfcn, lines, txtedtdat, usefigdat, linkselect )
% exploreResults launches an interactive plot from which data points can be
% selected and detailed results for that case viewed.
% By: Kaden Agrey
% v1.1 2018.07.11
% 
% Input
%   mainplot: This can be a figure handle or anonymous function that
%   returns a figure (with no input arguments). The script will 
%   automatically add space for the for the ui elements.
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
%   instead of aligns them.
%   TODO: Find out how to force Matlab to wait for the renderer before
%   executing the mainplot formatting.
% 
%   pbtnfcn: This argument is used to initialize the push button called
%   'View Details'. It should be a cell array with an anonymous function
%   handle as the first element and arguments to the function as later
%   elements.
%   The first 4 arguments to this function MUST be reserved for
%   exploreResults. The latter arguments should correspond to those in the
%   cell array. 
%       function argout = myFunc(src, event, slct, ui, arg1, arg2, ... )
%       src: matlab variable pointing to the src of the callback.
%       event: matlab variable giving event information
%       slct: exploreResults variable with information on the selected
%       points. 
%           slct(i).xplobj: graphics object with explorable data
%           slct(i).pnt: list of graphics objects for selected points
%           slct(i).ind: list of selected linear indices from explorable 
%           data
%       ui: exploreResults variable with pointers to all ui elements
%       created
% 
%   TODO: Currently only one push button can be set - allow multiple
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
%   usefigdat: Additionally the user can allow the function to
%   automatically get the data from the axes and line objects by setting 
%   this option to true. This will create two text/edit pairs for each 
%   "Explorable" axes.
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

%% --- Initialize --- %%
% These variables define the size and spacing of the ui objects
pbtn_h = 30; % [pixels] PushButton height
pbtn_blw_marg = 5; % [pixels] Margin below PushButton
txtedt_h = 20*2; % [pixels] Text & Edit box height
txtedtmarg = [3 10 0 10]; % [pixels] Text-Edit box margins

% Generate or pass on the figure
if isa(mainfig, 'function_handle')
    fig = mainfig();
elseif isgraphics(mainfig,'figure')
    fig = mainfig;
else
    error('The main plot variable must be a figure or an anonymous function that will generate a figure')
end

% Initialize xplr struct with figure objects and slct struct with graphics 
% objects to select points from. ui struct is also initialized, each field
% holds the ui objects by associated graphics object. (push buttons for the
% figure, text/edit boxes for each axes they are placed under etc.)

% Get all axes parenting explorable lines
xplr = struct('ax', [], 'ln', lines);
xplr.ax = getExplorableAx(fig, lines);

% Initialize the struct to reference all ui objects. Each field will be an
% array with numel(ui.field) equal to the number of explorable objects in
% the corresponding field in xplr. (Eg, one figure, possibly several axes
% and several lines for each axes). Stored in each field will be 
ui = struct('fig', [], 'ax', []);

% Get objects that selected points will be associated with
slct = struct('pnt', [], 'xplobj', [], 'ind', []);
for p = 1:length(xplr.ln)
    slct(p).xplobj = xplr.ln(p);
end

%% --- Prepair Main Figure --- %%
% --- Figure size and plot positioning --- %
% TODO: clean up algorithm
% Resize the figure and position the plots to fit the ui
children = findobj(fig.Children, 'flat', '-not', 'AxisLocationMode', 'auto');
nch = length(children); % number of children in figure
if nch > 1
    cbmap = cellfun(@(C) strcmp(C,'colorbar'), get(children, 'Type'));
else
    cbmap = strcmp(get(children, 'Type'),'colorbar');
end
fig.Units = 'pixels';
set(fig.Children, 'Units', 'pixels');

% Get size of all the figure children
ax_tins = zeros(nch,4);
for ch = 1:nch
    % Some objects don't have a TightInset
    if isprop(children(ch), 'TightInset')
        ax_tins(ch,:) = children(ch).TightInset;
    else
        ax_tins(ch,:) = [0 0 0 0];
    end
end

ax_pos = getAsMat(children, 'Position');
ax_siz = ax_pos(:,3:4); % save plot sizes now - they can automagically change
ax_sub = getSubPlotInd(ax_pos(~cbmap,:)); % subplot-like indices for children (cell array)
ax_sub = assignCBSubPlotInd(children, ax_sub, cbmap);
ax_row = cellfun(@max, ax_sub(:,1)); % for our purposes we only want the lowest row each plot is part of
nrow = max(ax_row);

% Get rows that require padding (columns won't require padding)
pad_row = []; % list of rows requiring padding
for row = 1:nrow
    rowch = find(row==ax_row)';
    for ch = rowch

        % If the child is Explorable then add the row to the pad list
        ind = children(ch)==xplr.ax;
        if any(ind) && ( usefigdat || ~isempty(txtedtdat(ind)) )
            pad_row(length(pad_row)+1,1) = row;
            break
        end

    end
end

% Reposition children
% TODO: manually positioned colorbars shuoldn't always be aligned with the
% row they belong to - switch to an algorithm that offsets row members
% instead of aligns them.
ax_pos_new = ax_pos; % initially nothing has moved
for row = nrow:-1:1 % start from the bottom
    % Add text/edit padding if this row requires it
    if ismember(row, pad_row)
        axpad = txtedt_h + sum(txtedtmarg([2 4]));
    else
        axpad = txtedtmarg(2);
    end

    % Find the axes from which to get the height this row should start at
    rowch = find(row==ax_row); % children in this row
    rowcol = unique( cell2mat(ax_sub(rowch,2)) );

    rowt = max(ax_tins(rowch,:),[],1); % the largest TightInsets from this row (we only use the bottom and top)
    rowy = ax_pos(rowch(1),2) - rowt(2); % lowest y-position of this row

    % Note: here a 10 pixel overlap is allowed when checking if a figure is
    % below the row
    allsharecol = cellfun(@(C) any(ismember(C, rowcol)), ax_sub(:,2));
    allblw = ( sum(ax_pos(:,[2 4]), 2) + ax_tins(:,4) - 10 < rowy );
    chblw = find( allsharecol.*allblw ); % children below and sharing a column with this row

    if ~isempty(chblw)
        rowy = max( sum(ax_pos_new(chblw,[2 4]), 2) + ax_tins(chblw,4) ); % new y-position of this row
    else
        rowy = pbtn_h - txtedtmarg(2) + 5; % if this is the bottom row place it 5 pixels above the PushButtons
    end

    % Move all the axes of this row accordingly
    for ch = rowch'
        children(ch).Position(2) = rowy + rowt(2) + axpad;
    end

    % Update the axes positions recorded
    ax_pos_new = getAsMat(children, 'Position'); % ################## CHNAGE LATER
    ax_pos_new(:,3:4) = ax_siz; % maintain axes sizes

end

% Resize figure window
fig.OuterPosition(2) = 50; % move to the bottom of the screen
fig.Position(4) = max( sum(ax_pos_new(:,[2 4]), 2) + ax_tins(:,4) ) + 10;

% Reset the position sizes incase they have magically changed
for ch = 1:length(children)
    children(ch).Position(3:4) = ax_pos_new(ch, 3:4);
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
if isempty(txtedtdat)
    txtedtdat = cell(1,length(xplr.ax));
end

% Make button to view details of selected point
horz = getLeftChild(fig, 'pixels'); % place in line with farthest left plot

ui.fig.pbtn = uicontrol(fig, 'Style', 'pushbutton',... % make the button
                             'String', 'View Details',...
                             'Units', 'pixels',...
                             'Position', [0 0 150 pbtn_h]);

ui.fig.pbtn.Position = ui.fig.pbtn.Position + [horz pbtn_blw_marg 0 0];
% ui.fig.pbtn.Units = 'normalized'; % reset units

% Make text and edit objects for each explorable axes
nax = length(xplr.ax);
for a = 1:nax

    % Initial text and edit box size and options
    txtopt = {'Position', [0 0 75 txtedt_h/2]};
    editopt = {'Position', [0 0 75 txtedt_h/2], 'String', 'Empty', 'Enable', 'inactive'};

    if usefigdat % user setting to use figure axes data
        axdat = getAxesData(xplr.ax(a), xplr.ln);
        if isempty(txtedtdat{a})
            txtedtdat{a} = axdat;
        else
            txtedtdat{a} = [axdat, txtedtdat{a}];
        end
    end

    % Make the text/edit pair
    [ui.ax(a).txt, ui.ax(a).edt] = makeValueDispBar(xplr.ax(a), txtedtdat{a}, txtedtmarg, txtopt, editopt);

end

% Assign user callback function to push button
ui.fig.pbtn.Callback = {@ viewDetailsCallback, slct, ui, pbtnfcn};

% Assign ButtonDownFcn to selectable objects in figures
for ob = 1:length(slct)
    % Pass structs with selected point information (slct) and updatable ui 
    % objects (ui). Pass index of current selectable object.
    slct(ob).xplobj.ButtonDownFcn = {@ lnChoosePnt, ob, slct, ui, linkselect};
end

%% Finalize Figure Properties
% Reset units on figure and all figure children
fig.Units = 'normalized';
set(fig.Children, 'Units', 'normalized');

end

%% --- Functions --- %%
% --- Getters --- %
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

function [axs, lns] = getExplorable(fig)
% Gets pointers to axes and lines from parent figure marked 'explorable'

ca = 1; cl = 1; % counters
% Loop over all children of parent figure
for a = length(fig.Children):-1:1
% We loop backwards because figure children are in reverse order of
% creation, and it's likely that the users data is organized in order of
% creation.

    % If user set 'Tag' field to 'Explorable' add child to axes list
    if strcmp(fig.Children(a).Tag, 'Explorable')

        axs(ca) = fig.Children(a); % Get axes pointer

        % Loop over all children of parent axes
        for lnn = length(axs(ca).Children):-1:1 % looping backwards for the same reason
            % If user set 'Tag' field to 'Explorable' add child to line list
            if strcmp(axs(ca).Children(lnn).Tag, 'Explorable')
                lns(cl) = axs(ca).Children(lnn); % Get line pointer
                cl = cl + 1; % increment line counter
            end
        end

        ca = ca + 1; % increment axes counter

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

% --- Makers --- %
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

% --- Callbacks --- %
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
ui.fig.pbtn.Callback{2} = slct;

end

function [] = viewDetailsCallback(src, event, slct, ui, usrcall)
% Callback function for the view details push button. Takes the currently
% selected point information and the user supplied anonymous function with
% arguments. After checking that points have been properly selected will
% launch the user function with the first three arguments as matlab defined
% <src>, <event>, and explore_results defined <slct>.

% If points haven't all been selected do nothing.
if ~isfield(slct, 'ind')
    return
end
for p = 1:length(slct)
    if isempty(slct(p).ind)
        return
    end
end

% Call the user's function and pass options through
usrcall{1}( src, event, slct, ui, usrcall{2:end} );

end
