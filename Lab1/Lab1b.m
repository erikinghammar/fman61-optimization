classdef Lab1b < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure           matlab.ui.Figure
        GridLayout         matlab.ui.container.GridLayout
        LeftPanel          matlab.ui.container.Panel
        GridLayout2        matlab.ui.container.GridLayout
        SDLinesearchmathbfx_n1mathbfx_nlambda_nmathbfd_nLabel  matlab.ui.control.Label
        Label              matlab.ui.control.Label
        FValue             matlab.ui.control.TextArea
        NewtonXnLabel_2    matlab.ui.control.Label
        Iternum            matlab.ui.control.TextArea
        NewtonXnLabel      matlab.ui.control.Label
        FunDropDown        matlab.ui.control.DropDown
        fxyLabel           matlab.ui.control.Label
        MethodButtonGroup  matlab.ui.container.ButtonGroup
        Epsilon            matlab.ui.control.Slider
        varepsilonLabel    matlab.ui.control.Label
        SD                 matlab.ui.control.RadioButton
        NM                 matlab.ui.control.RadioButton
        Ymax               matlab.ui.control.NumericEditField
        YmaxLabel          matlab.ui.control.Label
        Ymin               matlab.ui.control.NumericEditField
        YminLabel          matlab.ui.control.Label
        Xmax               matlab.ui.control.NumericEditField
        XmaxLabel          matlab.ui.control.Label
        Xmin               matlab.ui.control.NumericEditField
        XminLabel          matlab.ui.control.Label
        RightPanel         matlab.ui.container.Panel
        Instructions       matlab.ui.control.Label
        LevelsSlider       matlab.ui.control.Slider
        LevelsSliderLabel  matlab.ui.control.Label
        Axes3D             matlab.ui.control.UIAxes
        Axes2D             matlab.ui.control.UIAxes
        ContextMenu        matlab.ui.container.ContextMenu
        RestartMenu        matlab.ui.container.Menu
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end


    properties (Access = private)
        F = @(x,y) exp(x) + x.^2 + 2*x.*y + 4*y.^2 + 10*(x + 2*y - 6).^2;
        DF = matlabFunction(gradient(sym(@(x,y) exp(x) + x.^2 + 2*x.*y + 4*y.^2 + 10*(x + 2*y - 6).^2)),'Vars',{'x','y'})
        DDF = matlabFunction(hessian(sym(@(x,y) exp(x) + x.^2 + 2*x.*y + 4*y.^2 + 10*(x + 2*y - 6).^2)),'Vars',{'x','y'})

        x0
        y0

        FunLine2D
        FunLine3D

        Fun2D
        Fun3D

        Niter = 0;
        restart = true;
        
        draw3 = true;
    end

    methods (Access = private)

        function drawFunction(app)
            delete(app.Fun2D)
            [xx,yy] = meshgrid(linspace(app.Xmin.Value,app.Xmax.Value),linspace(app.Ymin.Value,app.Ymax.Value));
            zz = app.F(xx,yy);

            shift = -1.1*min(min(zz(:)),0) + 1e-6;
            zz = zz + shift;
            tt = logspace(log10(min(zz(:))),log10(max(zz(:))),app.LevelsSlider.Value);

            tt = tt - shift;
            zz = zz - shift;
            [~,app.Fun2D] = contour(app.Axes2D,xx,yy,zz,tt);
            app.Fun2D.HitTest = 'off';
            colorbar(app.Axes2D); colormap(app.Axes2D,'winter');

            if app.draw3
                delete(app.Fun3D)
                [xx,yy] = meshgrid(linspace(app.Xmin.Value,app.Xmax.Value,30),linspace(app.Ymin.Value,app.Ymax.Value,30));
                zz = app.F(xx,yy);
                app.Fun3D = mesh(app.Axes3D,xx,yy,zz);
                app.Fun3D.FaceAlpha = .2;
                shading(app.Axes3D,'interp')
                colorbar(app.Axes3D); colormap(app.Axes3D,'winter');
                view(app.Axes3D,3);
            end
        end

        function NewtonIteration(app)
            A = app.Epsilon.Value*eye(2) + app.DDF(app.x0,app.y0);
            d = -A\app.DF(app.x0,app.y0);
            app.x0 = app.x0 + d(1); app.y0 = app.y0 + d(2);
            UpdatePlots(app)
        end

        function SteepestDescentIteration(app)
            d = -app.DF(app.x0,app.y0);
            d = d/norm(d);
            lambda = 1;
            alpha = 2;
            epsilon = 0.2;

            T = @(L) app.F(app.x0,app.y0) + L*epsilon*d'*app.DF(app.x0,app.y0);
            F1 = @(L) app.F(app.x0 + L*d(1),app.y0 + L*d(2));

            while lambda < 100 && lambda > 1e-4
                if F1(lambda) > T(lambda)
                    lambda = lambda/alpha;
                    continue
                end

                if F1(alpha*lambda) < T(alpha*lambda)
                    lambda = lambda*alpha;
                    continue
                end
                break;
            end

            app.x0 = app.x0 + lambda*d(1);
            app.y0 = app.y0 + lambda*d(2);
            UpdatePlots(app)
        end

        function UpdatePlots(app)
            app.FunLine2D.XData = [app.FunLine2D.XData app.x0];
            app.FunLine2D.YData = [app.FunLine2D.YData app.y0];

            app.FunLine3D.XData = app.x0;
            app.FunLine3D.YData = app.y0;
            app.FunLine3D.ZData = app.F(app.x0,app.y0);

            app.Niter = app.Niter + 1;
            app.Iternum.Value = num2str(app.Niter);

            app.FValue.Value = num2str(app.F(app.x0,app.y0));
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            warning('off','MATLAB:RankDeficientMatrix')
            warning('off','MATLAB:ui:Slider:fixedHeight')
            app.restart = true;
            delete([app.FunLine2D app.FunLine3D])
            app.Iternum.Value = '';
            app.FValue.Value = '';
            drawFunction(app)
        end

        % Value changed function: FunDropDown
        function FunDropDownValueChanged(app, event)
            value = app.FunDropDown.Value;
            switch value
                case "Function 1"
                    app.F = @(x,y) exp(x) + x.^2 + 2*x.*y + 4*y.^2 + 10*(x + 2*y - 6).^2;
                    app.LevelsSlider.Value = 10;
                case "Function 2"
                    app.F = @(x,y) sqrt(x.^2 + (y - 2).^2) + sqrt(x.^2 + y.^2) + sqrt((x - 6).^2 + (y + 2).^2);
                    app.LevelsSlider.Value = 10;
                case "Function 3"
                    app.F = @(x,y) 100*(y - x.^2).^2 + (1-x).^2 + 2;
                    app.LevelsSlider.Value = 2;
                case "Function 4"
                    app.F = @(x,y) (0.35*(11*x + 6)./12).^3.*(0.2 + 0.35*(8*y - 24)/12).*exp((-(0.35*(11*x + 6)/12).^2 - 2*(0.35*(8*y - 24)/12).^2)/2);
                    app.LevelsSlider.Value = 2;
            end
            startupFcn(app);
            app.DF = matlabFunction(gradient(sym(app.F)),'Vars',{'x','y'});
            app.DDF = matlabFunction(hessian(sym(app.F)),'Vars',{'x','y'});
        end

        % Button down function: Axes2D
        function Axes2DButtonDown(app, event)
            if event.Button == 1
                if app.restart
                    app.x0 = event.IntersectionPoint(1);
                    app.y0 = event.IntersectionPoint(2);
                    app.FunLine2D = plot(app.Axes2D,app.x0,app.y0,'rx--');
                    app.FunLine2D.HitTest = 'off';
                    app.restart = false;

                    app.FValue.Value = num2str(app.F(app.x0,app.y0));
                    app.FunLine3D = plot3(app.Axes3D,app.x0,app.y0,app.F(app.x0,app.y0),'r.','MarkerSize',25);
                else
                    switch app.MethodButtonGroup.SelectedObject.Tag
                        case "Newton"
                            NewtonIteration(app)
                        case "SD"
                            SteepestDescentIteration(app)
                    end
                end
            end
        end

        % Value changed function: LevelsSlider
        function LevelsSliderValueChanged(app, event)
            app.draw3 = false;
            drawFunction(app);
            app.draw3 = true;
        end

        % Menu selected function: RestartMenu
        function RestartMenuSelected(app, event)
            delete([app.FunLine2D,app.FunLine3D])
            app.Niter = 0;
            app.Iternum.Value = '';
            app.FValue.Value = '';
            startupFcn(app)
        end

        % Value changed function: Xmax, Xmin, Ymax, Ymin
        function XminValueChanged(app, event)
            drawFunction(app)
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {637, 637};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {276, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)
            warning('off','MATLAB:RankDeficientMatrix');
            warning('off','MATLAB:ui:Slider:fixedHeight');

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 854 637];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {276, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.LeftPanel);
            app.GridLayout2.ColumnWidth = {31, 67, '3x', '1x'};
            app.GridLayout2.RowHeight = {'100x', '50x', '50x', 25, 25, '1.46x', 24, 24, 24, 24, 'fit', 22, 'fit'};
            app.GridLayout2.ColumnSpacing = 4.8;
            app.GridLayout2.Padding = [4.8 10 4.8 10];

            % Create XminLabel
            app.XminLabel = uilabel(app.GridLayout2);
            app.XminLabel.Interpreter = 'latex';
            app.XminLabel.HorizontalAlignment = 'right';
            app.XminLabel.FontName = 'Courier New';
            app.XminLabel.FontSize = 16;
            app.XminLabel.Layout.Row = 7;
            app.XminLabel.Layout.Column = 2;
            app.XminLabel.Text = '$x_\text{min}$';

            % Create Xmin
            app.Xmin = uieditfield(app.GridLayout2, 'numeric');
            app.Xmin.ValueChangedFcn = createCallbackFcn(app, @XminValueChanged, true);
            app.Xmin.FontName = 'Courier New';
            app.Xmin.Layout.Row = 7;
            app.Xmin.Layout.Column = 3;
            app.Xmin.Value = -6;

            % Create XmaxLabel
            app.XmaxLabel = uilabel(app.GridLayout2);
            app.XmaxLabel.Interpreter = 'latex';
            app.XmaxLabel.HorizontalAlignment = 'right';
            app.XmaxLabel.FontSize = 16;
            app.XmaxLabel.Layout.Row = 8;
            app.XmaxLabel.Layout.Column = 2;
            app.XmaxLabel.Text = '$x_\text{max}$';

            % Create Xmax
            app.Xmax = uieditfield(app.GridLayout2, 'numeric');
            app.Xmax.ValueChangedFcn = createCallbackFcn(app, @XminValueChanged, true);
            app.Xmax.FontName = 'Courier New';
            app.Xmax.Layout.Row = 8;
            app.Xmax.Layout.Column = 3;
            app.Xmax.Value = 6;

            % Create YminLabel
            app.YminLabel = uilabel(app.GridLayout2);
            app.YminLabel.Interpreter = 'latex';
            app.YminLabel.HorizontalAlignment = 'right';
            app.YminLabel.FontSize = 16;
            app.YminLabel.Layout.Row = 9;
            app.YminLabel.Layout.Column = 2;
            app.YminLabel.Text = '$y_\text{min}$';

            % Create Ymin
            app.Ymin = uieditfield(app.GridLayout2, 'numeric');
            app.Ymin.ValueChangedFcn = createCallbackFcn(app, @XminValueChanged, true);
            app.Ymin.FontName = 'Courier New';
            app.Ymin.Layout.Row = 9;
            app.Ymin.Layout.Column = 3;
            app.Ymin.Value = -6;

            % Create YmaxLabel
            app.YmaxLabel = uilabel(app.GridLayout2);
            app.YmaxLabel.Interpreter = 'latex';
            app.YmaxLabel.HorizontalAlignment = 'right';
            app.YmaxLabel.FontSize = 16;
            app.YmaxLabel.Layout.Row = 10;
            app.YmaxLabel.Layout.Column = 2;
            app.YmaxLabel.Text = '$y_\text{max}$';

            % Create Ymax
            app.Ymax = uieditfield(app.GridLayout2, 'numeric');
            app.Ymax.ValueChangedFcn = createCallbackFcn(app, @XminValueChanged, true);
            app.Ymax.FontName = 'Courier New';
            app.Ymax.Layout.Row = 10;
            app.Ymax.Layout.Column = 3;
            app.Ymax.Value = 6;

            % Create MethodButtonGroup
            app.MethodButtonGroup = uibuttongroup(app.GridLayout2);
            app.MethodButtonGroup.AutoResizeChildren = 'off';
            app.MethodButtonGroup.Title = 'Choose a method';
            app.MethodButtonGroup.Layout.Row = 1;
            app.MethodButtonGroup.Layout.Column = [1 4];

            % Create NM
            app.NM = uiradiobutton(app.MethodButtonGroup);
            app.NM.Tag = 'Newton';
            app.NM.Text = 'Modified Newton method (without l.s.)';
            app.NM.FontName = 'Verdana';
            app.NM.Position = [11 108 258 22];
            app.NM.Value = true;

            % Create SD
            app.SD = uiradiobutton(app.MethodButtonGroup);
            app.SD.Tag = 'SD';
            app.SD.Text = 'Steepest descent (with Armijo''s l.s.)';
            app.SD.FontName = 'Verdana';
            app.SD.Position = [9 8 244 22];

            % Create varepsilonLabel
            app.varepsilonLabel = uilabel(app.MethodButtonGroup);
            app.varepsilonLabel.Interpreter = 'latex';
            app.varepsilonLabel.HorizontalAlignment = 'right';
            app.varepsilonLabel.FontSize = 16;
            app.varepsilonLabel.Position = [24 68 25 22];
            app.varepsilonLabel.Text = '$\varepsilon$';

            % Create Epsilon
            app.Epsilon = uislider(app.MethodButtonGroup);
            app.Epsilon.Limits = [0 30];
            app.Epsilon.Position = [70 77 133 7];

            % Create fxyLabel
            app.fxyLabel = uilabel(app.GridLayout2);
            app.fxyLabel.Interpreter = 'latex';
            app.fxyLabel.HorizontalAlignment = 'right';
            app.fxyLabel.FontSize = 16;
            app.fxyLabel.Layout.Row = 12;
            app.fxyLabel.Layout.Column = [1 2];
            app.fxyLabel.Text = '$f(x,y) =$';

            % Create FunDropDown
            app.FunDropDown = uidropdown(app.GridLayout2);
            app.FunDropDown.Items = {'Function 1', 'Function 2', 'Function 3', 'Function 4'};
            app.FunDropDown.ValueChangedFcn = createCallbackFcn(app, @FunDropDownValueChanged, true);
            app.FunDropDown.FontName = 'Verdana';
            app.FunDropDown.Layout.Row = [11 13];
            app.FunDropDown.Layout.Column = 3;
            app.FunDropDown.Value = 'Function 1';

            % Create NewtonXnLabel
            app.NewtonXnLabel = uilabel(app.GridLayout2);
            app.NewtonXnLabel.Interpreter = 'latex';
            app.NewtonXnLabel.HorizontalAlignment = 'right';
            app.NewtonXnLabel.FontName = 'Courier New';
            app.NewtonXnLabel.FontSize = 16;
            app.NewtonXnLabel.Layout.Row = 4;
            app.NewtonXnLabel.Layout.Column = [1 2];
            app.NewtonXnLabel.Text = 'No. iterations';

            % Create Iternum
            app.Iternum = uitextarea(app.GridLayout2);
            app.Iternum.FontName = 'Courier New';
            app.Iternum.FontSize = 16;
            app.Iternum.Layout.Row = 4;
            app.Iternum.Layout.Column = 3;

            % Create NewtonXnLabel_2
            app.NewtonXnLabel_2 = uilabel(app.GridLayout2);
            app.NewtonXnLabel_2.Interpreter = 'latex';
            app.NewtonXnLabel_2.HorizontalAlignment = 'right';
            app.NewtonXnLabel_2.FontSize = 16;
            app.NewtonXnLabel_2.Layout.Row = 5;
            app.NewtonXnLabel_2.Layout.Column = [1 2];
            app.NewtonXnLabel_2.Text = '$f(x_n,y_n) = $';

            % Create FValue
            app.FValue = uitextarea(app.GridLayout2);
            app.FValue.FontName = 'Courier New';
            app.FValue.FontSize = 16;
            app.FValue.Layout.Row = 5;
            app.FValue.Layout.Column = 3;

            % Create Label
            app.Label = uilabel(app.GridLayout2);
            app.Label.Interpreter = 'latex';
            app.Label.HorizontalAlignment = 'center';
            app.Label.Layout.Row = 2;
            app.Label.Layout.Column = [1 4];
            app.Label.Text = {'Modified Newton '; '$\mathbf{x}_{n+1} = \mathbf{x}_n - (\varepsilon_n \mathbb{I} + \nabla^2f(\mathbf{x}_n))^{-1}\nabla f(\mathbf{x}_n)$'};

            % Create SDLinesearchmathbfx_n1mathbfx_nlambda_nmathbfd_nLabel
            app.SDLinesearchmathbfx_n1mathbfx_nlambda_nmathbfd_nLabel = uilabel(app.GridLayout2);
            app.SDLinesearchmathbfx_n1mathbfx_nlambda_nmathbfd_nLabel.Interpreter = 'latex';
            app.SDLinesearchmathbfx_n1mathbfx_nlambda_nmathbfd_nLabel.HorizontalAlignment = 'center';
            app.SDLinesearchmathbfx_n1mathbfx_nlambda_nmathbfd_nLabel.Layout.Row = 3;
            app.SDLinesearchmathbfx_n1mathbfx_nlambda_nmathbfd_nLabel.Layout.Column = [1 4];
            app.SDLinesearchmathbfx_n1mathbfx_nlambda_nmathbfd_nLabel.Text = {'SD + Line search'; '$\mathbf{d}_n = -\nabla f(\mathbf{x}_n)$'; '$\lambda_n = {\rm{argmin}\atop \scriptstyle \lambda} f(\mathbf{x}_n + \lambda\mathbf{d}_n)$'; '$\mathbf{x}_{n+1} = \mathbf{x}_n + \lambda_n \mathbf{d}_n$'};

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Create Axes2D
            app.Axes2D = uiaxes(app.RightPanel);
            app.Axes2D.NextPlot = 'add';
            app.Axes2D.ButtonDownFcn = createCallbackFcn(app, @Axes2DButtonDown, true);
            app.Axes2D.Position = [7 319 566 263];

            % Create Axes3D
            app.Axes3D = uiaxes(app.RightPanel);
            app.Axes3D.NextPlot = 'add';
            app.Axes3D.Position = [7 1 566 274];

            % Create LevelsSliderLabel
            app.LevelsSliderLabel = uilabel(app.RightPanel);
            app.LevelsSliderLabel.HorizontalAlignment = 'right';
            app.LevelsSliderLabel.Position = [3 605 40 22];
            app.LevelsSliderLabel.Text = 'Levels';

            % Create LevelsSlider
            app.LevelsSlider = uislider(app.RightPanel);
            app.LevelsSlider.Limits = [2 60];
            app.LevelsSlider.ValueChangedFcn = createCallbackFcn(app, @LevelsSliderValueChanged, true);
            app.LevelsSlider.Position = [64 615 497 3];
            app.LevelsSlider.Value = 10;

            % Create Instructions
            app.Instructions = uilabel(app.RightPanel);
            app.Instructions.HorizontalAlignment = 'center';
            app.Instructions.FontSize = 16;
            app.Instructions.Position = [19 283 533 19];
            app.Instructions.Text = 'Right-click and choose "Restart" to reset the program!';

            % Create ContextMenu
            app.ContextMenu = uicontextmenu(app.UIFigure);

            % Create RestartMenu
            app.RestartMenu = uimenu(app.ContextMenu);
            app.RestartMenu.MenuSelectedFcn = createCallbackFcn(app, @RestartMenuSelected, true);
            app.RestartMenu.Text = 'Restart';
            
            % Assign app.ContextMenu
            app.RightPanel.ContextMenu = app.ContextMenu;
            app.Axes2D.ContextMenu = app.ContextMenu;
            app.Axes3D.ContextMenu = app.ContextMenu;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Lab1b

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end