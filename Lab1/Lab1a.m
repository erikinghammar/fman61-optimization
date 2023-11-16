classdef Lab1a < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure               matlab.ui.Figure
        GridLayout             matlab.ui.container.GridLayout
        LeftPanel              matlab.ui.container.Panel
        TabGroup               matlab.ui.container.TabGroup
        GoldenTab              matlab.ui.container.Tab
        GoldenGrid             matlab.ui.container.GridLayout
        GoldenIternum          matlab.ui.control.TextArea
        NewtonXnLabel_4        matlab.ui.control.Label
        GoldenFL               matlab.ui.control.TextArea
        NewtonYnLabel_5        matlab.ui.control.Label
        GoldenFM               matlab.ui.control.TextArea
        NewtonXnLabel_3        matlab.ui.control.Label
        GoldenL                matlab.ui.control.TextArea
        NewtonY0Label_2        matlab.ui.control.Label
        GoldenTheory           matlab.ui.control.Label
        NewtonTab              matlab.ui.container.Tab
        GridLayout2            matlab.ui.container.GridLayout
        NewtonIternum          matlab.ui.control.TextArea
        NewtonXnLabel_2        matlab.ui.control.Label
        NewtonXn               matlab.ui.control.TextArea
        NewtonDn               matlab.ui.control.TextArea
        NewtonYn               matlab.ui.control.TextArea
        NewtonYnLabel          matlab.ui.control.Label
        NewtonXnLabel          matlab.ui.control.Label
        NewtonY0Label          matlab.ui.control.Label
        NewtonTheory           matlab.ui.control.Label
        ArmijoTab              matlab.ui.container.Tab
        GridLayout3            matlab.ui.container.GridLayout
        alphaLabel             matlab.ui.control.Label
        ArmijoAlpha            matlab.ui.control.Slider
        ArmijoEpsilon          matlab.ui.control.Slider
        varepsilonSliderLabel  matlab.ui.control.Label
        ArmijoTheory1          matlab.ui.control.Label
        ArmijoTheory2          matlab.ui.control.Label
        Xmax                   matlab.ui.control.NumericEditField
        XmaxLabel              matlab.ui.control.Label
        Xmin                   matlab.ui.control.NumericEditField
        XminLabel              matlab.ui.control.Label
        fxDropDown             matlab.ui.control.DropDown
        fxDropDownLabel        matlab.ui.control.Label
        RightPanel             matlab.ui.container.Panel
        Instructions           matlab.ui.control.Label
        UIAxes                 matlab.ui.control.UIAxes
        ContextMenu            matlab.ui.container.ContextMenu
        RestartMenu            matlab.ui.container.Menu
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    
    properties (Access = private)
        
        F = @(x) (-x.^4+6*x.^2)./(1+(x-0.5).^4)
        DF = matlabFunction(diff(sym(@(x) (-x.^4+6*x.^2)./(1+(x-0.5).^4))))
        DDF = matlabFunction(diff(sym(@(x) (-x.^4+6*x.^2)./(1+(x-0.5).^4)),2))
        iter = 0
        
        phi = .5*(sqrt(5) - 1);
        xa = -Inf
        xb = Inf
        GoldenLine1
        GoldenLine2
        GoldenLine3

        x0 = -Inf
        NewtonRestart = true
        NewtonLine

        lambda0 = 0.02
        ArmijoStartPoint
        ArmijoRedLine
        ArmijoGreenLine
        ArmijoBlueLine

        
        Flabels = ["\frac{-x^4+6x^2}{1 + (x - 0.5)^4}";
                   "\max\left\{\frac{(x+7)^2}{50}+4,10-\frac{(x+10)^2}{50})\right\}"]
        FunLine
    end
    
    methods (Access = private)

        function UpdatePlot(app)
            delete(app.FunLine)
            xx = linspace(app.UIAxes.XLim(1),app.UIAxes.XLim(2),500);
            app.FunLine = plot(app.UIAxes,xx,app.F(xx),'b-');
            app.FunLine.HitTest = 'off';

            YL = max(app.F(xx)) - min(app.F(xx));
            Y0 = (max(app.F(xx)) + min(app.F(xx)))/2;
            R = 1.1;

            app.UIAxes.YLim = [-R/2*YL + Y0, R/2*YL + Y0];
        end
        
        function UpdateXLim(app)
            X = [app.Xmin.Value, app.Xmax.Value];
            if X(2) > X(1)
                app.UIAxes.XLim = X;
                UpdatePlot(app)
            else
                return
            end
        end
        
        function GoldenStep(app)
            delete([app.GoldenLine1,app.GoldenLine2,app.GoldenLine3]);
            L = app.xb - app.xa;
            lk = app.xa + (1 - app.phi)*L;
            mk = app.xa + app.phi*L;

            app.GoldenLine1 = plot(app.UIAxes,[app.xa app.xa],[app.UIAxes.YLim(1) app.F(app.xa)],'ro-');
            app.GoldenLine3 = plot(app.UIAxes,[app.xb app.xb],[app.UIAxes.YLim(1) app.F(app.xb)],'ro-');

            app.GoldenL.Value = num2str(app.xb - app.xa);
            
            app.GoldenFL.Value = num2str(app.F(lk));
            app.GoldenFM.Value = num2str(app.F(mk));
               
            if ~isinf(app.xa*app.xb)
                if app.F(lk) < app.F(mk)
                    app.xb = mk;
                    app.GoldenLine2 = plot(app.UIAxes,[app.xb app.xb],[app.UIAxes.YLim(1) app.F(app.xb)],'b*-');
                else
                    app.xa = lk;
                    app.GoldenLine2 = plot(app.UIAxes,[app.xa app.xa],[app.UIAxes.YLim(1) app.F(app.xa)],'b*-');
                end
                app.iter = app.iter + 1;
                app.GoldenIternum.Value = num2str(app.iter);
            end
        end

        function NewtonStep(app)
            if app.NewtonRestart
                app.NewtonRestart = false;
                app.NewtonLine = plot(app.UIAxes,app.x0,app.F(app.x0),'r.--','MarkerSize',14);
                app.NewtonXn.Value = num2str(app.x0);
            else
                app.iter = app.iter + 1;
                app.NewtonXn.Value = num2str(app.x0);
                app.x0 = app.x0 - app.DF(app.x0)/app.DDF(app.x0);
                app.NewtonLine.XData = [app.NewtonLine.XData app.x0];
                app.NewtonLine.YData = [app.NewtonLine.YData app.F(app.x0)];
                app.NewtonDn.Value = num2str(- app.DF(app.x0)/app.DDF(app.x0));
                app.NewtonYn.Value = num2str(app.F(app.x0));
                app.NewtonIternum.Value = num2str(app.iter);
            end
        end

        function ArmijoStep(app)
            if app.ArmijoEpsilon.Value >= 1 || app.ArmijoAlpha.Value <= 1
                return
            end

            delete([app.ArmijoBlueLine;app.ArmijoRedLine;app.ArmijoGreenLine;app.ArmijoStartPoint])
            d = -sign(app.DF(app.x0));
            T = @(L) app.F(app.x0) + app.ArmijoEpsilon.Value*L*app.DF(app.x0);

            xx = linspace(app.x0,d*max(d*app.UIAxes.XLim));
            app.ArmijoStartPoint = plot(app.UIAxes,app.x0,app.F(app.x0),'go');
            app.ArmijoGreenLine = plot(app.UIAxes,xx,T(xx - app.x0),'g-');

            lambda = d*app.lambda0;
            while true
                if app.F(app.x0 + lambda) > T(lambda)
                    %if (1) does not hold
                    lambda = lambda/app.ArmijoAlpha.Value;
                    continue
                end

                if app.F(app.x0 + app.ArmijoAlpha.Value*lambda) < T(app.ArmijoAlpha.Value*lambda)
                    %if (2) does not hold
                    lambda = lambda*app.ArmijoAlpha.Value;
                    continue
                end

                %if both (1) and (2) are satisfied
                break;
            end
            ArmijoIter = round(log(abs(lambda/app.lambda0))/log(app.ArmijoAlpha.Value));
            % if lambda was increasing, add extra redstar to plot
            if ArmijoIter > 0
                ArmijoIter = ArmijoIter + 1;
            end
            kk = sign(ArmijoIter)*(0:abs(ArmijoIter));
            LL = app.x0 + d*app.lambda0*(app.ArmijoAlpha.Value.^kk);

            app.ArmijoRedLine = plot(app.UIAxes,LL,app.F(LL),'r*');
            app.ArmijoBlueLine = plot(app.UIAxes,app.x0 + lambda,app.F(app.x0 + lambda),'b*');
            
            app.ArmijoBlueLine.HitTest = 'off';
            app.ArmijoRedLine.HitTest = 'off';
            app.ArmijoGreenLine.HitTest = 'off';
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.NewtonRestart = true;
            app.iter = 0;
            app.x0 = -Inf;
            app.xa = -Inf;
            app.xb = Inf;
            cla(app.UIAxes);
            %app.UIAxes.Title.String = "N. Iterations: 0";
            
            %fxDropDownValueChanged(app)
            app.fxDropDownLabel.Text = sprintf("$f(x) = %s$",app.Flabels(app.fxDropDown.Value(end) - 48));
            
            app.NewtonXn.Value = '';
            app.NewtonDn.Value = '';
            app.NewtonYn.Value = '';
            app.NewtonIternum.Value = '';

            app.ArmijoEpsilon.Value = .2;
            app.ArmijoAlpha.Value = 2;

            app.GoldenL.Value = '';
            app.GoldenFM.Value = '';
            app.GoldenFL.Value = '';
            app.GoldenIternum.Value = '';

            UpdatePlot(app)
        end

        % Button down function: UIAxes
        function UIAxesButtonDown(app, event)
            if event.Button == 1
                switch app.TabGroup.SelectedTab.Title
                    case "Newton"
                        if app.NewtonRestart
                            app.x0 = event.IntersectionPoint(1);
                        end
                        NewtonStep(app)
                        app.NewtonRestart = false;
                    case "Armijo"
                        app.x0 = event.IntersectionPoint(1);
                        ArmijoStep(app)
                    case "Golden"
                        if app.xa == -Inf
                            app.xa = event.IntersectionPoint(1);
                        elseif app.xb == Inf
                            if app.xa < event.IntersectionPoint(1)
                                app.xb = event.IntersectionPoint(1);
                            elseif app.xa > event.IntersectionPoint(1)
                                app.xb = app.xa;
                                app.xa = event.IntersectionPoint(1);
                            elseif app.xa == app.xb
                                return
                            end
                        end
                        GoldenStep(app)
                end
            end
        end

        % Value changed function: Xmin
        function XminValueChanged(app, event)
            UpdateXLim(app)
        end

        % Value changed function: Xmax
        function XmaxValueChanged(app, event)
            UpdateXLim(app)
        end

        % Callback function
        function MethodButtonGroupSelectionChanged(app, event)
            selectedButton = app.MethodButtonGroup.SelectedObject;
            newtonGroup = [app.NewtonX0;app.NewtonX0Label;
                                   app.NewtonXn;app.NewtonY0Label;
                                   app.NewtonDn;app.NewtonXnLabel;
                                   app.NewtonYn;app.NewtonYnLabel];
            armijoGroup = [];
            goldenGroup = [];

            switch selectedButton.Text
                case "Newton"
                    set(newtonGroup,'Visible','on');
                    set(armijoGroup,'Visible','off');
                    set(goldenGroup,'Visible','off');
                    startupFcn(app)
                case "Armijo"
                    set(newtonGroup,'Visible','off');
                    set(armijoGroup,'Visible','on');
                    set(goldenGroup,'Visible','off');
                    startupFcn(app)
                case "Golden section"
                    startupFcn(app)
                    set(newtonGroup,'Visible','off');
                    set(armijoGroup,'Visible','off');
                    set(goldenGroup,'Visible','on');
                    startupFcn(app)
            end
        end

        % Callback function: RestartMenu, TabGroup
        function RestartMenuSelected(app, event)
            app.NewtonRestart = true;
            startupFcn(app)
        end

        % Value changed function: fxDropDown
        function fxDropDownValueChanged(app, event)
            value = app.fxDropDown.Value;
            switch value
                % there must be an easier way to deal with symbolics
                % and then convert at the end...
                case "Function 1"
                    app.F = @(x) (-x.^4+6*x.^2)./(1+(x-0.5).^4);
                    app.Xmin.Value = -0.5;
                    app.Xmax.Value = 0.5;
                case "Function 2"
                    app.F = @(x) abs(x.^2/25 + (17*x)/25 - 151/50)/2 - (3*x)/50 + 649/100;
                    app.Xmin.Value = -10;
                    app.Xmax.Value = 10;
            end
            UpdateXLim(app)
            app.DF = matlabFunction(diff(sym(app.F)));
            app.DDF = matlabFunction(diff(sym(app.F),2));
            startupFcn(app)
        end

        % Value changed function: ArmijoAlpha, ArmijoEpsilon
        function ArmijoEpsilonValueChanged(app, event)
            if app.x0 > -Inf
                ArmijoStep(app)
            end
        end

        % Value changing function: ArmijoEpsilon
        function ArmijoEpsilonValueChanging(app, event)
            app.ArmijoEpsilon.Value = event.Value;
            if app.x0 > -Inf
                ArmijoStep(app)
                drawnow
            end
        end

        % Value changing function: ArmijoAlpha
        function ArmijoAlphaValueChanging(app, event)
            app.ArmijoAlpha.Value = event.Value;
            if app.x0 > -Inf
                ArmijoStep(app)
                drawnow
            end
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {571, 571};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {262, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 949 571];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {262, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create fxDropDownLabel
            app.fxDropDownLabel = uilabel(app.LeftPanel);
            app.fxDropDownLabel.Interpreter = 'latex';
            app.fxDropDownLabel.HorizontalAlignment = 'center';
            app.fxDropDownLabel.Position = [7 164 249 56];
            app.fxDropDownLabel.Text = '$f(x)$';

            % Create fxDropDown
            app.fxDropDown = uidropdown(app.LeftPanel);
            app.fxDropDown.Items = {'Function 1', 'Function 2'};
            app.fxDropDown.ValueChangedFcn = createCallbackFcn(app, @fxDropDownValueChanged, true);
            app.fxDropDown.FontName = 'Verdana';
            app.fxDropDown.Position = [57 115 161 28];
            app.fxDropDown.Value = 'Function 1';

            % Create XminLabel
            app.XminLabel = uilabel(app.LeftPanel);
            app.XminLabel.Interpreter = 'latex';
            app.XminLabel.Tag = 'xminlabel';
            app.XminLabel.HorizontalAlignment = 'right';
            app.XminLabel.Position = [86 53 32 22];
            app.XminLabel.Text = '$x_{\rm min}$';

            % Create Xmin
            app.Xmin = uieditfield(app.LeftPanel, 'numeric');
            app.Xmin.ValueChangedFcn = createCallbackFcn(app, @XminValueChanged, true);
            app.Xmin.Tag = 'xmin';
            app.Xmin.FontName = 'Courier New';
            app.Xmin.Position = [133 54 56 22];
            app.Xmin.Value = -0.5;

            % Create XmaxLabel
            app.XmaxLabel = uilabel(app.LeftPanel);
            app.XmaxLabel.Interpreter = 'latex';
            app.XmaxLabel.HorizontalAlignment = 'right';
            app.XmaxLabel.Position = [86 24 32 22];
            app.XmaxLabel.Text = '$x_{\rm max}$';

            % Create Xmax
            app.Xmax = uieditfield(app.LeftPanel, 'numeric');
            app.Xmax.ValueChangedFcn = createCallbackFcn(app, @XmaxValueChanged, true);
            app.Xmax.FontName = 'Courier New';
            app.Xmax.Position = [133 25 56 22];
            app.Xmax.Value = 0.5;

            % Create TabGroup
            app.TabGroup = uitabgroup(app.LeftPanel);
            app.TabGroup.SelectionChangedFcn = createCallbackFcn(app, @RestartMenuSelected, true);
            app.TabGroup.Position = [1 224 260 346];

            % Create GoldenTab
            app.GoldenTab = uitab(app.TabGroup);
            app.GoldenTab.AutoResizeChildren = 'off';
            app.GoldenTab.Title = 'Golden';

            % Create GoldenGrid
            app.GoldenGrid = uigridlayout(app.GoldenTab);
            app.GoldenGrid.RowHeight = {'100x', 25, 25, 25, 25};

            % Create GoldenTheory
            app.GoldenTheory = uilabel(app.GoldenGrid);
            app.GoldenTheory.Interpreter = 'latex';
            app.GoldenTheory.HorizontalAlignment = 'center';
            app.GoldenTheory.FontSize = 18;
            app.GoldenTheory.Layout.Row = 1;
            app.GoldenTheory.Layout.Column = [1 2];
            app.GoldenTheory.Text = {'$\varphi := \frac{\sqrt{5} - 1}{2}$'; '$\lambda_k = a_k + (1 - \varphi)(b_k - a_k)$'; '$\mu_k = a_k + \varphi(b_k - a_k)$'};

            % Create NewtonY0Label_2
            app.NewtonY0Label_2 = uilabel(app.GoldenGrid);
            app.NewtonY0Label_2.Interpreter = 'latex';
            app.NewtonY0Label_2.HorizontalAlignment = 'right';
            app.NewtonY0Label_2.FontName = 'Courier New';
            app.NewtonY0Label_2.FontSize = 16;
            app.NewtonY0Label_2.Layout.Row = 2;
            app.NewtonY0Label_2.Layout.Column = 1;
            app.NewtonY0Label_2.Text = '$L$';

            % Create GoldenL
            app.GoldenL = uitextarea(app.GoldenGrid);
            app.GoldenL.FontName = 'Courier New';
            app.GoldenL.FontSize = 14;
            app.GoldenL.Layout.Row = 2;
            app.GoldenL.Layout.Column = 2;

            % Create NewtonXnLabel_3
            app.NewtonXnLabel_3 = uilabel(app.GoldenGrid);
            app.NewtonXnLabel_3.Interpreter = 'latex';
            app.NewtonXnLabel_3.HorizontalAlignment = 'right';
            app.NewtonXnLabel_3.FontName = 'Courier New';
            app.NewtonXnLabel_3.FontSize = 16;
            app.NewtonXnLabel_3.Layout.Row = 3;
            app.NewtonXnLabel_3.Layout.Column = 1;
            app.NewtonXnLabel_3.Text = '$f(\lambda)$';

            % Create GoldenFM
            app.GoldenFM = uitextarea(app.GoldenGrid);
            app.GoldenFM.FontName = 'Courier New';
            app.GoldenFM.FontSize = 14;
            app.GoldenFM.Layout.Row = 3;
            app.GoldenFM.Layout.Column = 2;

            % Create NewtonYnLabel_5
            app.NewtonYnLabel_5 = uilabel(app.GoldenGrid);
            app.NewtonYnLabel_5.Interpreter = 'latex';
            app.NewtonYnLabel_5.HorizontalAlignment = 'right';
            app.NewtonYnLabel_5.FontName = 'Courier New';
            app.NewtonYnLabel_5.FontSize = 16;
            app.NewtonYnLabel_5.Layout.Row = 4;
            app.NewtonYnLabel_5.Layout.Column = 1;
            app.NewtonYnLabel_5.Text = '$f(\mu)$';

            % Create GoldenFL
            app.GoldenFL = uitextarea(app.GoldenGrid);
            app.GoldenFL.FontName = 'Courier New';
            app.GoldenFL.FontSize = 14;
            app.GoldenFL.Layout.Row = 4;
            app.GoldenFL.Layout.Column = 2;

            % Create NewtonXnLabel_4
            app.NewtonXnLabel_4 = uilabel(app.GoldenGrid);
            app.NewtonXnLabel_4.Interpreter = 'latex';
            app.NewtonXnLabel_4.HorizontalAlignment = 'right';
            app.NewtonXnLabel_4.FontName = 'Courier New';
            app.NewtonXnLabel_4.FontSize = 16;
            app.NewtonXnLabel_4.Layout.Row = 5;
            app.NewtonXnLabel_4.Layout.Column = 1;
            app.NewtonXnLabel_4.Text = 'No. iterations';

            % Create GoldenIternum
            app.GoldenIternum = uitextarea(app.GoldenGrid);
            app.GoldenIternum.FontName = 'Courier New';
            app.GoldenIternum.FontSize = 14;
            app.GoldenIternum.Layout.Row = 5;
            app.GoldenIternum.Layout.Column = 2;

            % Create NewtonTab
            app.NewtonTab = uitab(app.TabGroup);
            app.NewtonTab.AutoResizeChildren = 'off';
            app.NewtonTab.Title = 'Newton';

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.NewtonTab);
            app.GridLayout2.RowHeight = {'1.13x', 25, 25, 25, 25};

            % Create NewtonTheory
            app.NewtonTheory = uilabel(app.GridLayout2);
            app.NewtonTheory.Interpreter = 'latex';
            app.NewtonTheory.HorizontalAlignment = 'center';
            app.NewtonTheory.FontSize = 18;
            app.NewtonTheory.Layout.Row = 1;
            app.NewtonTheory.Layout.Column = [1 2];
            app.NewtonTheory.Text = {'$\delta_{n} = - \frac{f''(x_n)}{f''''(x_n)}$'; ''; '$x_{n+1} = x_n + \delta_n$'};

            % Create NewtonY0Label
            app.NewtonY0Label = uilabel(app.GridLayout2);
            app.NewtonY0Label.Interpreter = 'latex';
            app.NewtonY0Label.HorizontalAlignment = 'right';
            app.NewtonY0Label.FontName = 'Courier New';
            app.NewtonY0Label.FontSize = 16;
            app.NewtonY0Label.Layout.Row = 2;
            app.NewtonY0Label.Layout.Column = 1;
            app.NewtonY0Label.Text = '$x_n$';

            % Create NewtonXnLabel
            app.NewtonXnLabel = uilabel(app.GridLayout2);
            app.NewtonXnLabel.Interpreter = 'latex';
            app.NewtonXnLabel.HorizontalAlignment = 'right';
            app.NewtonXnLabel.FontName = 'Courier New';
            app.NewtonXnLabel.FontSize = 16;
            app.NewtonXnLabel.Layout.Row = 3;
            app.NewtonXnLabel.Layout.Column = 1;
            app.NewtonXnLabel.Text = '$\delta_n$';

            % Create NewtonYnLabel
            app.NewtonYnLabel = uilabel(app.GridLayout2);
            app.NewtonYnLabel.Interpreter = 'latex';
            app.NewtonYnLabel.HorizontalAlignment = 'right';
            app.NewtonYnLabel.FontName = 'Courier New';
            app.NewtonYnLabel.FontSize = 16;
            app.NewtonYnLabel.Layout.Row = 4;
            app.NewtonYnLabel.Layout.Column = 1;
            app.NewtonYnLabel.Text = '$f(x_{n+1})$';

            % Create NewtonYn
            app.NewtonYn = uitextarea(app.GridLayout2);
            app.NewtonYn.FontName = 'Courier New';
            app.NewtonYn.FontSize = 14;
            app.NewtonYn.Layout.Row = 4;
            app.NewtonYn.Layout.Column = 2;

            % Create NewtonDn
            app.NewtonDn = uitextarea(app.GridLayout2);
            app.NewtonDn.FontName = 'Courier New';
            app.NewtonDn.FontSize = 14;
            app.NewtonDn.Layout.Row = 3;
            app.NewtonDn.Layout.Column = 2;

            % Create NewtonXn
            app.NewtonXn = uitextarea(app.GridLayout2);
            app.NewtonXn.FontName = 'Courier New';
            app.NewtonXn.FontSize = 14;
            app.NewtonXn.Layout.Row = 2;
            app.NewtonXn.Layout.Column = 2;

            % Create NewtonXnLabel_2
            app.NewtonXnLabel_2 = uilabel(app.GridLayout2);
            app.NewtonXnLabel_2.Interpreter = 'latex';
            app.NewtonXnLabel_2.HorizontalAlignment = 'right';
            app.NewtonXnLabel_2.FontName = 'Courier New';
            app.NewtonXnLabel_2.FontSize = 16;
            app.NewtonXnLabel_2.Layout.Row = 5;
            app.NewtonXnLabel_2.Layout.Column = 1;
            app.NewtonXnLabel_2.Text = 'No. iterations';

            % Create NewtonIternum
            app.NewtonIternum = uitextarea(app.GridLayout2);
            app.NewtonIternum.FontName = 'Courier New';
            app.NewtonIternum.FontSize = 14;
            app.NewtonIternum.Layout.Row = 5;
            app.NewtonIternum.Layout.Column = 2;

            % Create ArmijoTab
            app.ArmijoTab = uitab(app.TabGroup);
            app.ArmijoTab.Title = 'Armijo';

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.ArmijoTab);
            app.GridLayout3.ColumnWidth = {25, '100x'};
            app.GridLayout3.RowHeight = {'2x', '2x', 45, 45};
            app.GridLayout3.ColumnSpacing = 5.13751220703125;
            app.GridLayout3.Padding = [5.13751220703125 10 5.13751220703125 10];

            % Create ArmijoTheory2
            app.ArmijoTheory2 = uilabel(app.GridLayout3);
            app.ArmijoTheory2.Interpreter = 'latex';
            app.ArmijoTheory2.HorizontalAlignment = 'center';
            app.ArmijoTheory2.FontSize = 18;
            app.ArmijoTheory2.Layout.Row = 1;
            app.ArmijoTheory2.Layout.Column = [1 2];
            app.ArmijoTheory2.Text = {'$F(\lambda) := f(\mathbf{x}_0 + \lambda\mathbf{d})$'; ''; '$T(\lambda) = F(0) + \varepsilon \lambda F''(0)$'};

            % Create ArmijoTheory1
            app.ArmijoTheory1 = uilabel(app.GridLayout3);
            app.ArmijoTheory1.Interpreter = 'latex';
            app.ArmijoTheory1.FontSize = 18;
            app.ArmijoTheory1.Layout.Row = 2;
            app.ArmijoTheory1.Layout.Column = [1 2];
            app.ArmijoTheory1.Text = {'Condition 1: $F(\bar{\lambda}) \leq T(\bar{\lambda})$'; ''; 'Condition 2: $F(\alpha \bar{\lambda}) \geq T(\alpha \bar{\lambda})$'};

            % Create varepsilonSliderLabel
            app.varepsilonSliderLabel = uilabel(app.GridLayout3);
            app.varepsilonSliderLabel.Interpreter = 'latex';
            app.varepsilonSliderLabel.HorizontalAlignment = 'right';
            app.varepsilonSliderLabel.FontSize = 16;
            app.varepsilonSliderLabel.Layout.Row = 3;
            app.varepsilonSliderLabel.Layout.Column = 1;
            app.varepsilonSliderLabel.Text = '$\varepsilon$';

            % Create ArmijoEpsilon
            app.ArmijoEpsilon = uislider(app.GridLayout3);
            app.ArmijoEpsilon.Limits = [1e-06 1];
            app.ArmijoEpsilon.ValueChangedFcn = createCallbackFcn(app, @ArmijoEpsilonValueChanged, true);
            app.ArmijoEpsilon.ValueChangingFcn = createCallbackFcn(app, @ArmijoEpsilonValueChanging, true);
            app.ArmijoEpsilon.Layout.Row = 3;
            app.ArmijoEpsilon.Layout.Column = 2;
            app.ArmijoEpsilon.Value = 0.476703069318499;

            % Create ArmijoAlpha
            app.ArmijoAlpha = uislider(app.GridLayout3);
            app.ArmijoAlpha.Limits = [1 3];
            app.ArmijoAlpha.MajorTicks = [1 1.5 2 2.5 3];
            app.ArmijoAlpha.MajorTickLabels = {'1', '1.5', '2', '2.5', '3'};
            app.ArmijoAlpha.ValueChangedFcn = createCallbackFcn(app, @ArmijoEpsilonValueChanged, true);
            app.ArmijoAlpha.ValueChangingFcn = createCallbackFcn(app, @ArmijoAlphaValueChanging, true);
            app.ArmijoAlpha.MinorTicks = [1 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2 2.05 2.1 2.15 2.2 2.25 2.3 2.35 2.4 2.45 2.5 2.55 2.6 2.65 2.7 2.75 2.8 2.85 2.9 2.95 3];
            app.ArmijoAlpha.Layout.Row = 4;
            app.ArmijoAlpha.Layout.Column = 2;
            app.ArmijoAlpha.Value = 1.5;

            % Create alphaLabel
            app.alphaLabel = uilabel(app.GridLayout3);
            app.alphaLabel.Interpreter = 'latex';
            app.alphaLabel.HorizontalAlignment = 'right';
            app.alphaLabel.FontSize = 16;
            app.alphaLabel.Layout.Row = 4;
            app.alphaLabel.Layout.Column = 1;
            app.alphaLabel.Text = '$\alpha$';

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Create UIAxes
            app.UIAxes = uiaxes(app.RightPanel);
            app.UIAxes.XLim = [-0.5 0.5];
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.NextPlot = 'add';
            app.UIAxes.ButtonDownFcn = createCallbackFcn(app, @UIAxesButtonDown, true);
            app.UIAxes.Position = [6 47 675 513];

            % Create Instructions
            app.Instructions = uilabel(app.RightPanel);
            app.Instructions.HorizontalAlignment = 'center';
            app.Instructions.FontSize = 16;
            app.Instructions.Position = [6 16 675 22];
            app.Instructions.Text = 'Right-click and choose "Restart" to reset the program!';

            % Create ContextMenu
            app.ContextMenu = uicontextmenu(app.UIFigure);

            % Create RestartMenu
            app.RestartMenu = uimenu(app.ContextMenu);
            app.RestartMenu.MenuSelectedFcn = createCallbackFcn(app, @RestartMenuSelected, true);
            app.RestartMenu.Text = 'Restart';
            
            % Assign app.ContextMenu
            app.UIAxes.ContextMenu = app.ContextMenu;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Lab1a

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