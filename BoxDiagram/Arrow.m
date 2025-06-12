classdef Arrow < handle
    %ARROW Connect two Box objects.  Magnitude is a data value which the
    %layout engine can translate into a line‑width; LineWidth is filled in
    %later.  Arrowhead geometry matches the original ShowBoxPlot script.

    properties
        Source Box
        Target Box
        Label string = ""
        Magnitude double {mustBeNonnegative} = 1   % data value
        LineWidth double = NaN                     % data-space width, set later
        Color (1,3) double = [0 0 0]
        LabelOffset double = 0.5
        FontSize double {mustBePositive} = 10
    end

    properties (Constant, Access=private)
        BaseWidthPt = 6;    % head size in printer points when LineWidth = 1 pt
        WidthRatio  = 0.8;  % head base width : head length
    end

    methods
        function obj = Arrow(src,dst,options)
            arguments
                src Box
                dst Box
                options.Label string = ""
                options.Magnitude double = 1
                options.LineWidth double = NaN
                options.Color (1,3) double = [0 0 0]
                options.LabelOffset double = 0.5
                options.FontSize double = 10
            end
            obj.Source      = src;
            obj.Target      = dst;
            obj.Label       = options.Label;
            obj.Magnitude   = options.Magnitude;
            obj.LineWidth   = options.LineWidth;
            obj.Color       = options.Color;
            obj.LabelOffset = options.LabelOffset;
            obj.FontSize    = options.FontSize;
        end

        function draw(obj,ax)
            if isnan(obj.LineWidth) || obj.LineWidth<=0
                error("Arrow:LineWidthUnset","LineWidth must be assigned before drawing");
            end
            if nargin<2 || isempty(ax), ax = gca; end

            %%  Head dimensions in data units (proportional to LineWidth)
            headLen = (Arrow.BaseWidthPt/72) * obj.LineWidth; % approx
            headWidth = Arrow.WidthRatio * headLen;

            %%  Geometry
            [xc1,yc1] = obj.Source.getGeometry();
            [xc2,yc2] = obj.Target.getGeometry();
            ang = atan2(yc2-yc1, xc2-xc1);
            [x1,y1] = obj.Source.edgeIntersect(ang);
            [x2,y2] = obj.Target.edgeIntersect(ang+pi);

            %%  Shaft (trimmed before head)
            x2b = x2 - headLen*cos(ang);
            y2b = y2 - headLen*sin(ang);
            plot(ax,[x1 x2b],[y1 y2b],'Color',obj.Color,'LineWidth',obj.LineWidth);

            %%  Head triangle
            left  = [x2b + headWidth/2*sin(ang), y2b - headWidth/2*cos(ang)];
            right = [x2b - headWidth/2*sin(ang), y2b + headWidth/2*cos(ang)];
            patch(ax,[x2 left(1) right(1)],[y2 left(2) right(2)],obj.Color,'EdgeColor',obj.Color);

            %%  Optional label
            if strlength(obj.Label)>0
                tx = x1 + obj.LabelOffset*(x2 - x1);
                ty = y1 + obj.LabelOffset*(y2 - y1);
                text(ax,tx,ty,obj.Label,'HorizontalAlignment','center', ...
                    'VerticalAlignment','bottom','FontSize',obj.FontSize, ...
                    'BackgroundColor',[1 1 1 0.7],'Margin',0.5,'Interpreter','none');
            end
        end
    end
end