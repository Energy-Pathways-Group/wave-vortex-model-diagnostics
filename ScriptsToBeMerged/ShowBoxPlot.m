% % % function ShowBoxPlot
% function ShowBoxPlot(fluxes and reservoirs)

% inputs from ShowEnergyFluxesFromDiagnosticsFile.m
% could assemble these inputs directly from DiagnosticsFile, but I've
% already done a bunch of the work in ShowEnergyFluxesFromDiagnosticsFile
% so might as well call ShowBoxPlot from there.
meanToGeostrophic = (...
        sum(EnergyForcing_jkR{'te_g_geostrophic_mean_flow'}(:))...
        +sum(EnergyForcing_jkR{'te_mda_geostrophic_mean_flow'}(:))...
        )*flux_conversion;
nioToWave = (...
        sum(EnergyForcing_jkR{'te_igw_inertial_forcing'}(:))...
        +sum(EnergyForcing_jkR{'te_io_inertial_forcing'}(:))...
        )*flux_conversion;
tideToWave = (...
        sum(EnergyForcing_jkR{'te_igw_M2_tidal_forcing'}(:))...
        +sum(EnergyForcing_jkR{'te_io_M2_tidal_forcing'}(:))...
        )*flux_conversion;
geostrophicToFriction = (...
        sum(EnergyForcing_jkR{'te_g_quadratic_bottom_friction'}(:))...
        +sum(EnergyForcing_jkR{'te_mda_quadratic_bottom_friction'}(:))...
        )*flux_conversion;
geostrophicToDamping = (...
        sum(EnergyForcing_jkR{'te_g_adaptive_damping'}(:))...
        +sum(EnergyForcing_jkR{'te_mda_adaptive_damping'}(:))...
        )*flux_conversion;
waveToFriction = (...
        sum(EnergyForcing_jkR{'te_igw_quadratic_bottom_friction'}(:))...
        +sum(EnergyForcing_jkR{'te_io_quadratic_bottom_friction'}(:))...
        )*flux_conversion;
waveToDamping = (...
        sum(EnergyForcing_jkR{'te_igw_adaptive_damping'}(:))...
        +sum(EnergyForcing_jkR{'te_io_adaptive_damping'}(:))...
        )*flux_conversion;
geostrophicReservoir = (filter(E_g)+filter(E_mda))/E_gm;
waveReservoir = (filter(E_w)+filter(E_io))/E_gm;
% geostrophicToWave = % don't know what to include here.

% % dummy input
% meanToGeostrophic = 1;
% nioToWave = 1;
% tideToWave = 1;
% geostrophicToFriction = 1;
% geostrophicToDamping = 1;
% waveToFriction = 1;
% waveToDamping = 1;
% geostrophicReservoir = 1;
% waveReservoir = 1;

%==========================================================================
% 0. MAKE TITLE
%==========================================================================
titlePos = {
    'geostrophic-mean-flow', 'mean flow';
    'inertial-forcing', 'NIO';
    'M2-tidal-forcing', 'M2 tide'
};
if isa(wvt,'WVTransformBoussinesq')
    titleText = "Boussinesq";
elseif isa(wvt,'WVTransformHydrostatic')
    titleText = "Hydrostatic";
elseif isa(wvt,'WVTransformStratifiedQG')
    titleText = "Stratified QG";
elseif isa(wvt,'WVTransformBarotropicQG')
    titleText = "Barotropic QG";
elseif isa(wvt,'WVTransformConstantStratification')
    titleText = "Constant Stratification";
else
    error("Unrecognized WVTransform type.")
end
for i = 1:size(titlePos,1)
    if any(contains(wvt.forcingNames, titlePos{i,1}))
        titleText = titleText+", "+titlePos{i,2};  % add to titleText
    end
end

%==========================================================================
% 1. MASTER COLOUR TABLE
%==========================================================================
col.sources     = [191 191 250]/255;   % blue
col.geo         = [205 253 254]/255;   % green
col.wave        = [205 253 197]/255;   % yellow-orange
col.sinks       = [245 194 193]/255;   % red

%==========================================================================
% 2. DEFINE BOXES BY ROW  (label, colour, rounded?, boxText)
%==========================================================================

% static row1
% row1 = { ...
%     {'mean flow',       col.sources, 0, 'mean flow'}, ...
%     {'M2 tide',         col.sources, 0, 'M2 tide'}, ...
%     {'NIO',             col.sources, 0, 'NIO'}...
%     };

% build row1 of boxes depending on forcing applied
% List of possible row1 entries: {test_string, label, color, rounded?, boxText}
row1pos = {
    'geostrophic-mean-flow', 'mean flow',   col.sources, 0, 'mean flow';
    'inertial-forcing',      'NIO',         col.sources, 0, 'NIO';
    'M2-tidal-forcing',      'M2 tide',     col.sources, 0, 'M2 tide'
};
row1 = {};
for i = 1:size(row1pos,1)
    if any(contains(wvt.forcingNames, row1pos{i,1}))
        row1{end+1} = row1pos(i,2:5);  % add {label, color, rounded?, boxText}
    end
end
% row2 reservoirs
row2 = { ...
    {'geostrophic',     col.geo,    1, sprintf('geostrophic\n %0.2f GM',geostrophicReservoir)}, ...
    {'wave',            col.wave,   1, sprintf('wave\n %0.2f GM',waveReservoir)}...
    };
% row3 sinks
row3 = { ...
    {'bottom friction', col.sinks,  0, 'bottom friction'}, ...
    {'damping',         col.sinks,  0, 'damping'}...
    };
rowDef = {row1; row2; row3};          % 3×1 cell array (one element per row)

rowLbl = {'Sources','Reservoirs','Sinks'};


%==========================================================================
% 3. DEFINE ARROWS  – ONLY BY BOX LABELS
%    {fromLabel, toLabel, value, text, labPos}
%==========================================================================

% static arrowDef
% arrowDef = { ...
%     {'mean flow','geostrophic', meanToGeostrophic,sprintf('%0.2f GM/yr',meanToGeostrophic),0.5}, ...
%     {'NIO',      'wave',        nioToWave,sprintf('%0.2f GM/yr',nioToWave),0.5}, ...
%     {'geostrophic','wave',        1.1,'1.1',0.5}, ...
%     {'M2 tide',  'wave',        tideToWave,sprintf('%0.2f GM/yr',tideToWave),0.5}, ...
%     {'geostrophic','bottom friction',abs(geostrophicToFriction),sprintf('%0.2f GM/yr',abs(geostrophicToFriction)),0.5}, ...
%     {'geostrophic','damping',       abs(geostrophicToDamping),sprintf('%0.2f GM/yr',abs(geostrophicToDamping)),0.25}, ...
%     {'wave','bottom friction',      abs(waveToFriction),sprintf('%0.2f GM/yr',abs(waveToFriction)),0.25}, ...
%     {'wave','damping',              abs(waveToDamping),sprintf('%0.2f GM/yr',abs(waveToDamping)),0.5} };

% build ArrowDef depending on forcing applied
% start with required fluxes
arrowDef = { ...
    {'geostrophic','wave',        1.1,'1.1',0.5}, ...
    {'geostrophic','bottom friction',abs(geostrophicToFriction),sprintf('%0.2f GM/yr',abs(geostrophicToFriction)),0.5}, ...
    {'geostrophic','damping',       abs(geostrophicToDamping),sprintf('%0.2f GM/yr',abs(geostrophicToDamping)),0.25}, ...
    {'wave','bottom friction',      abs(waveToFriction),sprintf('%0.2f GM/yr',abs(waveToFriction)),0.25}, ...
    {'wave','damping',              abs(waveToDamping),sprintf('%0.2f GM/yr',abs(waveToDamping)),0.5} };
% add additional fluxes depending on forcing
arrowPos = {
    'geostrophic-mean-flow', 'mean flow','geostrophic', meanToGeostrophic,sprintf('%0.2f GM/yr',meanToGeostrophic),0.5;
    'inertial-forcing',      'NIO',      'wave',        nioToWave,sprintf('%0.2f GM/yr',nioToWave),0.5;
    'M2-tidal-forcing',      'M2 tide',  'wave',        tideToWave,sprintf('%0.2f GM/yr',tideToWave),0.5
};
for i = 1:size(arrowPos,1)
    if any(contains(wvt.forcingNames, arrowPos{i,1}))
        arrowDef{end+1} = arrowPos(i,2:6);  % add {label, color, rounded?, boxText}
    end
end

%==========================================================================
% 4. LAYOUT CONSTANTS (easy to tweak)
%==========================================================================
w = 2.2;                    % box width
h = [1,1.5,1];              % box height for each row
if length(row1)>2; dx=.5; else; dx=1.5; end % inter-box gap
dy = 2;                     % inter-row gap
x0 = 0;                     % left anchor
y0 = sum(h)+dy*(numel(rowDef)-1)+2;    % top anchor.  Sum of heights in h and dy, plus 2.
font.box = 16;  font.row = 18;  font.arrow = 14;  font.title = 20;
arrow.baseWidth  = 1.5;    % (already present)
arrow.headBasePt = 6;     % point size when shaft is baseWidth
arrow.headCutFrac = 1.0;   % 1 = stop at the head's base, 0.5 = middle

%==========================================================================
% 5. BUILD BOX STRUCTURE WITH POSITIONS
%==========================================================================
boxes = struct('label',[],'boxText',[],'pos',[],'faceColor',[],'round',[]);
rows  = cell(numel(rowDef),1);
boxPtr = 0;

for r = 1:numel(rowDef)
    row = rowDef{r};                     % the cell array for this row
    n   = numel(row);
    xs  = x0 + ((max(cellfun(@numel,rowDef))-n)*(w+dx))/2 + (0:n-1)*(w+dx);
    y   = y0 - (r-1)*dy - sum(h(1:r));
    rows{r} = [];
    for k = 1:n
        spec = row{k};                   % {label, color, rounded?, boxText}
        boxPtr = boxPtr + 1;
        boxes(boxPtr).label     = spec{1};
        boxes(boxPtr).boxText   = spec{4};
        boxes(boxPtr).faceColor = spec{2};
        boxes(boxPtr).round     = spec{3};
        boxes(boxPtr).pos       = [xs(k), y, w, h(r)];
        rows{r}(end+1) = boxPtr;
    end
end

% manually put "geostrophic" and "bottom friction" boxes under "mean flow"
if length(row1)==3
    boxes(4).pos(1) = boxes(1).pos(1);
    boxes(6).pos(1) = boxes(1).pos(1);
end

% Helper maps: label → index
label2idx = containers.Map( ...
    {boxes.label}, num2cell(1:numel(boxes)));

%==========================================================================
% 6. TRANSLATE ARROW DEFINITIONS INTO INDEXED STRUCT
%==========================================================================
arrows = struct('from',[],'to',[],'value',[],'txt',[],'labPos',[]);
for n = 1:numel(arrowDef)
    arrows(n).from  = label2idx(arrowDef{n}{1});
    arrows(n).to    = label2idx(arrowDef{n}{2});
    arrows(n).value = arrowDef{n}{3};
    arrows(n).txt   = arrowDef{n}{4};
    arrows(n).labPos = arrowDef{n}{5};
end

%==========================================================================
% 7. DRAW EVERYTHING
%==========================================================================
figure('Color','w'); clf
ax = axes('Position',[0.05 0.05 0.9 0.9]); hold(ax,'on'); axis(ax,'equal','off');
title(titleText,'FontSize',font.title,'FontWeight','bold')

% Compute box centres
ctr = @(p)[p(1)+w/2, p(2)+p(4)/2];

% Arrow line-widths
lw = arrow.baseWidth + 6*[arrows.value]/max([arrows.value]);
% lw = arrow.baseWidth + sqrt([arrows.value]/min([arrows.value]));

% Draw arrows
for n = 1:numel(arrows)

    % from and to box indices
    i = arrows(n).from;  j = arrows(n).to;

    % box centers
    p0 = ctr(boxes(i).pos);
    p1 = ctr(boxes(j).pos);
    v  = p1 - p0;
    p0 = p0 + 0.9*edgeOffset(v,[w/2,boxes(i).pos(4)/2]); % reduce edgeOffset at p0 and draw boxes over arrows because sometimes arrows didn't touch rounded box corners with full edgeOffset.
    p1 = p1 - edgeOffset(v,[w/2,boxes(j).pos(4)/2]);
    % plot(ax,[p0(1) p1(1)],[p0(2) p1(2)],'k','LineWidth',lw(n));
    % 
    % % scale head length ≈ proportional to the shaft's actual width
    % headSizePt = arrow.headBasePt * (lw(n) / arrow.baseWidth);
    % head = arrowHead(p0, p1, headSizePt);
    % patch(ax,head(:,1),head(:,2),'k','EdgeColor','none');

    % --- scale head to this arrow's width --------------------------------
    headSizePt = arrow.headBasePt * (lw(n) / arrow.baseWidth);
    headLen    = headSizePt / 72;              % pt → axis units
    
    % --- shorten shaft so it ends *before* the head ----------------------
    vHat   = (p1 - p0) / norm(p1 - p0);        % unit direction vector
    p1Base = p1 - arrow.headCutFrac * headLen * vHat;
    
    plot(ax,[p0(1) p1Base(1)],[p0(2) p1Base(2)],'k','LineWidth',lw(n));
    
    % --- draw arrowhead (tip at p1) --------------------------------------
    headTri = arrowHead(p1, vHat, headLen);    % helper updated below
    patch(ax,headTri(:,1),headTri(:,2),'k','EdgeColor','none');
    % label with white background
    pmid = p0 + arrows(n).labPos * (p1 - p0);
    text(ax,pmid(1),pmid(2),arrows(n).txt,'FontSize',font.arrow, ...
        'BackgroundColor','w','Margin',0.5, ...
        'HorizontalAlignment','center','VerticalAlignment','middle');
end

% Boxes
for k = 1:numel(boxes)
    r  = boxes(k).pos;
    hBox = boxes(k).pos(4);
    curv = boxes(k).round * [0.2 0.2*w/hBox];
    rectangle(ax,'Position',r,'Curvature',curv,'FaceColor',boxes(k).faceColor, ...
              'EdgeColor','k','LineWidth',1.2);
    text(ax,r(1)+w/2,r(2)+hBox/2,boxes(k).boxText,'FontSize',font.box, ...
        'HorizontalAlignment','center','VerticalAlignment','middle');
end

% row labels
xLabelOffset = 0.5;        % how far left of the left-most box (in axis units)
for r = 1:numel(rowDef)
    idx  = rows{r};                                % indices of boxes in row r
    yMid = boxes(idx(1)).pos(2) + h(r)/2;             % vertical centre of the row
    xLeft= 0*min(arrayfun(@(i) boxes(i).pos(1), idx)) - xLabelOffset;

    text(xLeft, yMid, rowLbl{r}, ...
         'FontSize', font.row, 'FontWeight', 'bold', ...
         'Rotation', 90, ...                       % 90° CCW
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment',   'middle');
end

%==========================================================================
% Helper functions
%==========================================================================
function d = edgeOffset(v,wh)                 % shortens arrow to meet box edge
    vu = v/norm(v);
    d  = vu * min([wh(1)/abs(vu(1)), wh(2)/abs(vu(2))]);
end
function tri = arrowHead(tip, dirHat, len)
    q = [-dirHat(2) dirHat(1)];            % perpendicular
    tri = [ tip; ...
            tip - len*(dirHat + 0.4*q); ...
            tip - len*(dirHat - 0.4*q) ];
end
% function tri = arrowHead(p0,p1,sizePt)       % small triangle
%     v = p1 - p0; v = v/norm(v);  q = [-v(2) v(1)];
%     L = sizePt/72; tri = [p1; p1-L*(v+0.4*q); p1-L*(v-0.4*q)];
% end

%==========================================================================
% % % end
