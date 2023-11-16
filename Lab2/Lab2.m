classdef Lab2 < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure               matlab.ui.Figure
        GridLayout             matlab.ui.container.GridLayout
        LeftPanel              matlab.ui.container.Panel
        LeftGrid               matlab.ui.container.GridLayout
        ObjectiveLabel         matlab.ui.control.Label
        LoadpenaltymButton     matlab.ui.control.Button
        LoadbarriermButton     matlab.ui.control.Button
        TheoryLabel            matlab.ui.control.Label
        FValue                 matlab.ui.control.TextArea
        FValueLabel            matlab.ui.control.Label
        BarrierPanel           matlab.ui.container.Panel
        BarrierGrid            matlab.ui.container.GridLayout
        EpsilonButton          matlab.ui.control.Button
        Epsilon                matlab.ui.control.NumericEditField
        EpsilonLabel           matlab.ui.control.Label
        BarrierIterations      matlab.ui.control.Label
        BarrierCheckBox        matlab.ui.control.CheckBox
        PenaltyPanel           matlab.ui.container.Panel
        PenaltyGrid            matlab.ui.container.GridLayout
        MuButton               matlab.ui.control.Button
        Mu                     matlab.ui.control.NumericEditField
        MuLabel                matlab.ui.control.Label
        PenaltyIterations      matlab.ui.control.Label
        PenaltyCheckBox        matlab.ui.control.CheckBox
        DropDown               matlab.ui.control.DropDown
        RightPanel             matlab.ui.container.Panel
        RightGrid              matlab.ui.container.GridLayout
        IterateButton          matlab.ui.control.Button
        z_0                    matlab.ui.control.NumericEditField
        z_0Label               matlab.ui.control.Label
        y_0                    matlab.ui.control.NumericEditField
        y_0Label               matlab.ui.control.Label
        x_0                    matlab.ui.control.NumericEditField
        x_0Label               matlab.ui.control.Label
        RightLabel1            matlab.ui.control.Label
        RightLabel2            matlab.ui.control.Label
        PersonalCheckBox       matlab.ui.control.CheckBox
        SystemCheckBox         matlab.ui.control.CheckBox
        UIAxes                 matlab.ui.control.UIAxes
        ContextMenu            matlab.ui.container.ContextMenu
        RestartMenu            matlab.ui.container.Menu
        ContextMenu2           matlab.ui.container.ContextMenu
        ShowcontoursonoffMenu  matlab.ui.container.Menu
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end


    properties (Access = private)
        A = [-1, -1;1, -1;1, 2;-1, 0;0, -1]
        b = [-2; 0; 6; 0; 0]
        c = [1;2];
        options = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton')

        X = Inf

        mu = 1
        epsilon = 1

        almostInf = 1e200;
        
        F
        IterLine
        StartLine
        ContourLine
        ContourBar

        NumIter = 0;

        Alpha
        Beta
    end

    methods (Access = private)

        function createPolygon(app)
            comb = nchoosek(1:size(app.A,1),size(app.A,2));
            v = [];
            for k = 1:size(comb,1)
                Ak = app.A(comb(k,:),:);
                bk = app.b(comb(k,:));
                if rank(Ak) == size(app.A,2)
                    xk = Ak\bk;
                    if all(app.A*xk - app.b <= 1e-12)
                        v = [v xk];
                    end
                end
            end

            switch size(app.A,2)
                case 2
                    set(app.UIAxes,'View',[0 90]);
                    [~,I] = unique(angle(v(1,:) + 1i*v(2,:) - mean(v(1,:) + 1i*v(2,:))));
                    fill(app.UIAxes,v(1,I),v(2,I),'y','FaceAlpha',0.5,'HitTest','off')
                case 3
                    v = unique(v.','rows');
                    set(app.UIAxes,'View',[-37.5, 30]);
                    plot(alphaShape(v(:,1),v(:,2),v(:,3)),'Parent',app.UIAxes,'FaceAlpha',0.5,'FaceColor','yellow','HitTest','off');
            end
        end

        function PenaltyStep(app)
            updateIteration(app)
            if app.SystemCheckBox.Value
                app.mu = 2*app.mu;
            end
            app.PenaltyIterations.Text = sprintf('No. Iterations: %i',app.NumIter);
            updateContour(app)
        end

        function BarrierStep(app)
            updateIteration(app)
            if app.SystemCheckBox.Value
                app.epsilon = app.epsilon/2;
            end
            app.BarrierIterations.Text = sprintf('No. Iterations: %i',app.NumIter);
            updateContour(app)
        end

        function chooseFunctions(app)
            if app.SystemCheckBox.Value
                app.Alpha = @(x) AG_penalty(app,x);
                app.Beta = @(x) AG_barrier(app,x);
                app.Mu.Visible = 'off';
                app.MuLabel.Visible = 'off';
                app.MuButton.Visible = 'off';
                app.Epsilon.Visible = 'off';
                app.EpsilonLabel.Visible = 'off';
                app.EpsilonButton.Visible = 'off';
                if app.DropDown.Value == "Problem 3"
                    app.IterateButton.Visible = 'on';
                end
            else
                app.IterateButton.Visible = 'off';
                app.Alpha = @(x) penalty(x,app.A,app.b);
                app.Beta = @(x) barrier(x,app.A,app.b);
                app.mu = 0;
                app.epsilon = 0;
                if app.PenaltyCheckBox.Value
                    app.Mu.Visible = 'on';
                    app.MuLabel.Visible = 'on';
                    app.MuButton.Visible = 'on';
                    app.Epsilon.Visible = 'off';
                    app.EpsilonLabel.Visible = 'off';
                    app.EpsilonButton.Visible = 'off';
                else
                    app.Mu.Visible = 'off';
                    app.MuLabel.Visible = 'off';
                    app.MuButton.Visible = 'off';
                    app.Epsilon.Visible = 'on';
                    app.EpsilonLabel.Visible = 'on';
                    app.EpsilonButton.Visible = 'on';
                end
            end
        end

        function a = AG_penalty(app,x)
            % function by Andrey Gulchek, 2016
            g = app.A*x - app.b;
            a = g'*diag(g > 0)*g;
        end

        function b = AG_barrier(app,x)
            % function by Andrey Gulchek, 2016
            g = app.A*x - app.b;
            if any(g >= 0)
                b = Inf;
            else
                b = -sum(log(-g));
            end
        end

        function loadFunction(~,funcname)
            [filename,pathname] = uigetfile(funcname,'Pick a MATLAB code file');
            while true
                if filename == 0
                    break
                elseif ~strcmp(filename,funcname)
                    uiwait(errordlg(sprintf('Please choose file with name %s',funcname)));
                    [filename,pathname] = uigetfile(funcname,'Pick a MATLAB code file');
                else
                    if ~strcmp([pwd '\'],pathname)
                        fullname = [pathname filename];
                        copyfile(fullname,funcname)
                    end
                    break
                end
            end
        end
        
        function updateIteration(app)
            updateF(app)
            app.X = fminunc(app.F,app.X,app.options);
            app.FValue.Value = num2str(app.c'*app.X);
            if app.DropDown.Value == "Problem 3"
                set(app.IterLine,...
                    'XData',[app.IterLine.XData, app.X(1)],...
                    'YData',[app.IterLine.YData, app.X(2)],...
                    'ZData',[app.IterLine.ZData, app.X(3)]);
            else
                set(app.IterLine,...
                    'XData',[app.IterLine.XData, app.X(1)],...
                    'YData',[app.IterLine.YData, app.X(2)]);
            end
            app.NumIter = app.NumIter + 1;
        end
        
        function [zz,tt,xx,yy] = computeContour(app)
            [xx,yy] = meshgrid(linspace(app.UIAxes.XLim(1),app.UIAxes.XLim(2),100),linspace(app.UIAxes.YLim(1),app.UIAxes.YLim(2)),100);
            zz = 0*xx;
            for i = 1:size(xx,1)
                for j = 1:size(xx,2)
                    zz(i,j) = app.F([xx(i,j);yy(i,j)]);
                end
            end

            shift = -1.1*min(min(zz(:)),0) + 1e-6;
            zz = zz + shift;
            
            zn = unique(zz(:));
            z1 = zn(1);
            if length(zn) > 1
                if app.BarrierCheckBox.Value
                    z2 = zn(end-1);
                else
                    z2 = zn(end);
                end
            else
                z2 = z1;
            end

            tt = logspace(log10(z1),log10(z2),20);

            zz = zz - shift;
            tt = tt - shift;
        end
        
        function updateF(app)
            if app.PenaltyCheckBox.Value
                app.F = @(x) app.c'*x + app.mu*app.Alpha(x);
            else
                app.F = @(x) min(app.c'*x + app.epsilon*app.Beta(x),app.almostInf);
            end
        end
        
        function updateContour(app)
            if any(isgraphics(app.ContourLine))
                [zz,tt,xx,yy] = computeContour(app);
                set(app.ContourLine,'XData',xx,...
                    'YData',yy,'ZData',zz,'LevelList',tt);
                drawnow
            end
        end
        
        function updateTheory(app)
            if app.PenaltyCheckBox.Value
                app.TheoryLabel.Text = "${\min \atop \scriptstyle \mathbf{x} \in \mathbb{R}^d} f(\mathbf{x}) + \mu_k \alpha(\mathbf{x})$";
            else
                app.TheoryLabel.Text = "${\min \atop \scriptstyle \mathbf{x} \in \mathrm{int}(S)} f(\mathbf{x}) + \varepsilon_k \beta(\mathbf{x})$";
            end
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.F = @(x) app.c'*x;
            createPolygon(app)
            chooseFunctions(app)
        end

        % Callback function: RestartMenu, x_0, y_0, z_0
        function RestartMenuSelected(app, event)
            delete([app.IterLine,app.StartLine])
            app.Mu.Value = app.mu;
            app.Epsilon.Value = app.epsilon;
            app.FValue.Value = "";
            app.NumIter = 0;
            app.PenaltyIterations.Text = "No. Iterations: 0";
            app.BarrierIterations.Text = "No. Iterations: 0";
            app.X = Inf;
            if app.DropDown.Value == "Problem 3"
                app.X = [app.x_0.Value;app.y_0.Value;app.z_0.Value];
                app.IterLine = plot3(app.UIAxes,app.X(1),app.X(2),app.X(3),'r.--','MarkerSize',12,'HitTest','off');
                app.FValue.Value = num2str(app.c'*app.X);
            end

            if app.SystemCheckBox.Value
                app.mu = 1;
                app.epsilon = 1;
                app.F = @(x) app.c'*x;
            end
            updateContour(app)
        end

        % Button down function: UIAxes
        function UIAxesButtonDown(app, event)
            if event.Button == 1 && ~strcmp(app.DropDown.Value,"Problem 3")
                if any(isinf(app.X))
                    app.X = [event.IntersectionPoint(1);event.IntersectionPoint(2)];
                    app.IterLine = plot(app.UIAxes,app.X(1),app.X(2),'rx--','HitTest','Off');
                    app.StartLine = plot(app.UIAxes,app.X(1),app.X(2),'ro','HitTest','Off');
                else
                    if app.PersonalCheckBox.Value
                        return
                    end
                    if app.BarrierCheckBox.Value
                        BarrierStep(app)
                    else
                        PenaltyStep(app)
                    end
                end                
                app.FValue.Value = num2str(app.c'*app.X);
            end
        end

        % Value changed function: PenaltyCheckBox
        function PenaltyCheckBoxValueChanged(app, event)
            app.BarrierCheckBox.Value = ~app.BarrierCheckBox.Value;
            updateTheory(app)
            chooseFunctions(app)
            RestartMenuSelected(app, event)
        end

        % Value changed function: BarrierCheckBox
        function BarrierCheckBoxValueChanged(app, event)
            app.PenaltyCheckBox.Value = ~app.PenaltyCheckBox.Value;
            updateTheory(app)
            chooseFunctions(app)
            RestartMenuSelected(app, event)
        end

        % Value changed function: SystemCheckBox
        function SystemCheckBoxValueChanged(app, event)
            app.PersonalCheckBox.Value = ~app.PersonalCheckBox.Value;
            chooseFunctions(app)
            RestartMenuSelected(app, event)
        end

        % Value changed function: PersonalCheckBox
        function PersonalCheckBoxValueChanged(app, event)
            app.SystemCheckBox.Value = ~app.SystemCheckBox.Value;
            chooseFunctions(app)
            RestartMenuSelected(app, event)
        end

        % Value changed function: DropDown
        function DropDownValueChanged(app, event)
            cla(app.UIAxes)
            switch app.DropDown.Value
                case "Problem 1"
                    rotate3d(app.UIAxes,'off');
                    app.A = [-1, -1;1, -1;1, 2;-1, 0;0, -1];
                    app.b = [-2; 0; 6; 0; 0];
                    app.c = [1;2];

                    app.x_0.Visible = 'off';
                    app.y_0.Visible = 'off';
                    app.z_0.Visible = 'off';
                    app.x_0Label.Visible = 'off';
                    app.y_0Label.Visible = 'off';
                    app.z_0Label.Visible = 'off';
                    app.IterateButton.Visible = 'off';
                    app.ObjectiveLabel.Text = "$f(x_1,x_2) = x_1 + 2x_2$";

                    app.UIAxes.XLim = [-2, 4];
                    app.UIAxes.YLim = [-1, 5];
                case "Problem 2"
                    rotate3d(app.UIAxes,'off');
                    app.A = [-1 3;1 -3;-1 -2;3 5;-1 0;0 -1];
                    app.b = [12;2;-2;34;0;0];
                    app.c = [1;2];

                    app.x_0.Visible = 'off';
                    app.y_0.Visible = 'off';
                    app.z_0.Visible = 'off';
                    app.x_0Label.Visible = 'off';
                    app.y_0Label.Visible = 'off';
                    app.z_0Label.Visible = 'off';
                    app.IterateButton.Visible = 'off';
                    app.ObjectiveLabel.Text = "$f(x_1,x_2) = x_1 + 2x_2$";

                    app.UIAxes.XLim = [-4, 12];
                    app.UIAxes.YLim = [-6, 10];
                case "Problem 3"
                    delete(app.ContourBar)
                    rotate3d(app.UIAxes,'on');
                    app.A = [1 3 2;1 1 1;3 5 3;-1, 0, 0;0, -1, 0;0, 0,-1];
                    app.b = [30; 24; 60; 0; 0; 0];
                    app.c = [-2;-4;-3];

                    app.x_0.Visible = 'on';
                    app.y_0.Visible = 'on';
                    app.z_0.Visible = 'on';
                    app.x_0Label.Visible = 'on';
                    app.y_0Label.Visible = 'on';
                    app.z_0Label.Visible = 'on';
                    if app.SystemCheckBox.Value
                        app.IterateButton.Visible = 'on';
                    end
                    app.ObjectiveLabel.Text = "$f(x_1,x_2,x_3) = -2x_1 -4x_2 -3x_3$";
                    app.UIAxes.XLim = [-10, 20];
                    app.UIAxes.YLim = [-10, 20];
                    app.UIAxes.ZLim = [-10, 20];
                    RestartMenuSelected(app, event)
            end
            RestartMenuSelected(app, event)
            createPolygon(app)
        end

        % Button pushed function: LoadbarriermButton
        function LoadBarrierButtonPushed(app, event)
            loadFunction(app,'barrier.m')
        end

        % Button pushed function: LoadpenaltymButton
        function LoadPenaltyButtonPushed(app, event)
            loadFunction(app,'penalty.m')
        end

        % Button pushed function: MuButton
        function MuButtonPushed(app, event)
            if ~any(isinf(app.X))
                PenaltyStep(app)
            end
            updateF(app)
            updateContour(app)
        end

        % Button pushed function: EpsilonButton
        function EpsilonButtonPushed(app, event)
            if ~any(isinf(app.X))
                BarrierStep(app)
            end
            updateF(app)
            updateContour(app)
        end

        % Button pushed function: IterateButton
        function IterateButtonPushed(app, event)
            if app.BarrierCheckBox.Value
                BarrierStep(app)
            else
                PenaltyStep(app)
            end
        end

        % Value changed function: Mu
        function MuValueChanged(app, event)
            app.mu = app.Mu.Value;
            %if app.PenaltyCheckBox.Value && ~any(isinf(app.X))
            %    PenaltyStep(app)
            %end
            %updateF(app)
            %updateContour(app)
        end

        % Value changed function: Epsilon
        function EpsilonValueChanged(app, event)
            app.epsilon = app.Epsilon.Value;
            %if app.BarrierCheckBox.Value && ~any(isinf(app.X))
            %    BarrierStep(app)
            %end
            %updateF(app)
            %updateContour(app)
        end

        % Menu selected function: ShowcontoursonoffMenu
        function ShowContoursMenuSelected(app, event)
            if app.DropDown.Value ~= "Problem 3"
                if ~any(isgraphics(app.ContourLine))
                    [zz,tt,xx,yy] = computeContour(app);
                    [~,app.ContourLine] = contour(app.UIAxes,xx,yy,zz,tt);
                    app.ContourLine.HitTest = 'off';
                    app.ContourBar = colorbar(app.UIAxes); colormap(app.UIAxes,'winter');
                else
                    delete([app.ContourLine,app.ContourBar])
                end
            end
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {569, 569};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {271, '1x'};
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
            app.UIFigure.Position = [100 100 763 569];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {271, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create LeftGrid
            app.LeftGrid = uigridlayout(app.LeftPanel);
            app.LeftGrid.RowHeight = {'2x', '2x', '2x', '2x', '1x', 20, 20, 25, 25};

            % Create DropDown
            app.DropDown = uidropdown(app.LeftGrid);
            app.DropDown.Items = {'Problem 1', 'Problem 2', 'Problem 3'};
            app.DropDown.ValueChangedFcn = createCallbackFcn(app, @DropDownValueChanged, true);
            app.DropDown.Layout.Row = 8;
            app.DropDown.Layout.Column = [1 2];
            app.DropDown.Value = 'Problem 1';

            % Create PenaltyPanel
            app.PenaltyPanel = uipanel(app.LeftGrid);
            app.PenaltyPanel.Title = 'Penalty method';
            app.PenaltyPanel.Layout.Row = [1 2];
            app.PenaltyPanel.Layout.Column = [1 2];

            % Create PenaltyGrid
            app.PenaltyGrid = uigridlayout(app.PenaltyPanel);
            app.PenaltyGrid.ColumnWidth = {'50x', 25, '50x', '50x'};
            app.PenaltyGrid.RowHeight = {'1x', 25, 25};

            % Create PenaltyCheckBox
            app.PenaltyCheckBox = uicheckbox(app.PenaltyGrid);
            app.PenaltyCheckBox.ValueChangedFcn = createCallbackFcn(app, @PenaltyCheckBoxValueChanged, true);
            app.PenaltyCheckBox.Text = 'Run penalty method';
            app.PenaltyCheckBox.Layout.Row = 1;
            app.PenaltyCheckBox.Layout.Column = [1 4];
            app.PenaltyCheckBox.Value = true;

            % Create PenaltyIterations
            app.PenaltyIterations = uilabel(app.PenaltyGrid);
            app.PenaltyIterations.Interpreter = 'latex';
            app.PenaltyIterations.Layout.Row = 3;
            app.PenaltyIterations.Layout.Column = [1 3];
            app.PenaltyIterations.Text = 'No. iterations: 0';

            % Create MuLabel
            app.MuLabel = uilabel(app.PenaltyGrid);
            app.MuLabel.Interpreter = 'latex';
            app.MuLabel.HorizontalAlignment = 'right';
            app.MuLabel.FontSize = 16;
            app.MuLabel.Visible = 'off';
            app.MuLabel.Layout.Row = 2;
            app.MuLabel.Layout.Column = 2;
            app.MuLabel.Text = '$\mu_k$';

            % Create Mu
            app.Mu = uieditfield(app.PenaltyGrid, 'numeric');
            app.Mu.ValueChangedFcn = createCallbackFcn(app, @MuValueChanged, true);
            app.Mu.FontName = 'Courier New';
            app.Mu.Visible = 'off';
            app.Mu.Layout.Row = 2;
            app.Mu.Layout.Column = [3 4];

            % Create MuButton
            app.MuButton = uibutton(app.PenaltyGrid, 'push');
            app.MuButton.ButtonPushedFcn = createCallbackFcn(app, @MuButtonPushed, true);
            app.MuButton.Visible = 'off';
            app.MuButton.Layout.Row = 3;
            app.MuButton.Layout.Column = [3 4];
            app.MuButton.Text = 'Iterate';

            % Create BarrierPanel
            app.BarrierPanel = uipanel(app.LeftGrid);
            app.BarrierPanel.Title = 'Barrier method';
            app.BarrierPanel.Layout.Row = [3 4];
            app.BarrierPanel.Layout.Column = [1 2];

            % Create BarrierGrid
            app.BarrierGrid = uigridlayout(app.BarrierPanel);
            app.BarrierGrid.ColumnWidth = {'50x', 25, '50x', '50x'};
            app.BarrierGrid.RowHeight = {'1x', 25, 25};

            % Create BarrierCheckBox
            app.BarrierCheckBox = uicheckbox(app.BarrierGrid);
            app.BarrierCheckBox.ValueChangedFcn = createCallbackFcn(app, @BarrierCheckBoxValueChanged, true);
            app.BarrierCheckBox.Text = 'Run barrier method';
            app.BarrierCheckBox.Layout.Row = 1;
            app.BarrierCheckBox.Layout.Column = [1 4];

            % Create BarrierIterations
            app.BarrierIterations = uilabel(app.BarrierGrid);
            app.BarrierIterations.Interpreter = 'latex';
            app.BarrierIterations.Layout.Row = 3;
            app.BarrierIterations.Layout.Column = [1 3];
            app.BarrierIterations.Text = 'No. iterations: 0';

            % Create EpsilonLabel
            app.EpsilonLabel = uilabel(app.BarrierGrid);
            app.EpsilonLabel.Interpreter = 'latex';
            app.EpsilonLabel.HorizontalAlignment = 'right';
            app.EpsilonLabel.FontSize = 16;
            app.EpsilonLabel.Visible = 'off';
            app.EpsilonLabel.Layout.Row = 2;
            app.EpsilonLabel.Layout.Column = 2;
            app.EpsilonLabel.Text = '$\varepsilon_k$';

            % Create Epsilon
            app.Epsilon = uieditfield(app.BarrierGrid, 'numeric');
            app.Epsilon.ValueChangedFcn = createCallbackFcn(app, @EpsilonValueChanged, true);
            app.Epsilon.FontName = 'Courier New';
            app.Epsilon.Visible = 'off';
            app.Epsilon.Layout.Row = 2;
            app.Epsilon.Layout.Column = [3 4];

            % Create EpsilonButton
            app.EpsilonButton = uibutton(app.BarrierGrid, 'push');
            app.EpsilonButton.ButtonPushedFcn = createCallbackFcn(app, @EpsilonButtonPushed, true);
            app.EpsilonButton.Visible = 'off';
            app.EpsilonButton.Layout.Row = 3;
            app.EpsilonButton.Layout.Column = [3 4];
            app.EpsilonButton.Text = 'Iterate';

            % Create FValueLabel
            app.FValueLabel = uilabel(app.LeftGrid);
            app.FValueLabel.Interpreter = 'latex';
            app.FValueLabel.HorizontalAlignment = 'right';
            app.FValueLabel.FontSize = 14;
            app.FValueLabel.Layout.Row = 7;
            app.FValueLabel.Layout.Column = 1;
            app.FValueLabel.Text = '$f(\mathbf{x}_k) =$';

            % Create FValue
            app.FValue = uitextarea(app.LeftGrid);
            app.FValue.Editable = 'off';
            app.FValue.FontName = 'Courier New';
            app.FValue.Layout.Row = 7;
            app.FValue.Layout.Column = 2;

            % Create TheoryLabel
            app.TheoryLabel = uilabel(app.LeftGrid);
            app.TheoryLabel.Interpreter = 'latex';
            app.TheoryLabel.HorizontalAlignment = 'center';
            app.TheoryLabel.FontSize = 15;
            app.TheoryLabel.Layout.Row = 5;
            app.TheoryLabel.Layout.Column = [1 2];
            app.TheoryLabel.Text = '${\min \atop \scriptstyle \mathbf{x} \in \mathbb{R}^d} f(\mathbf{x}) + \mu_k \alpha(\mathbf{x})$';

            % Create LoadbarriermButton
            app.LoadbarriermButton = uibutton(app.LeftGrid, 'push');
            app.LoadbarriermButton.ButtonPushedFcn = createCallbackFcn(app, @LoadBarrierButtonPushed, true);
            app.LoadbarriermButton.Layout.Row = 9;
            app.LoadbarriermButton.Layout.Column = 1;
            app.LoadbarriermButton.Text = 'Load barrier.m';

            % Create LoadpenaltymButton
            app.LoadpenaltymButton = uibutton(app.LeftGrid, 'push');
            app.LoadpenaltymButton.ButtonPushedFcn = createCallbackFcn(app, @LoadPenaltyButtonPushed, true);
            app.LoadpenaltymButton.Layout.Row = 9;
            app.LoadpenaltymButton.Layout.Column = 2;
            app.LoadpenaltymButton.Text = 'Load penalty.m';

            % Create ObjectiveLabel
            app.ObjectiveLabel = uilabel(app.LeftGrid);
            app.ObjectiveLabel.Interpreter = 'latex';
            app.ObjectiveLabel.HorizontalAlignment = 'center';
            app.ObjectiveLabel.FontSize = 16;
            app.ObjectiveLabel.Layout.Row = 6;
            app.ObjectiveLabel.Layout.Column = [1 2];
            app.ObjectiveLabel.Text = '$f(x_1,x_2) = x_1 + 2x_2$';

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Create RightGrid
            app.RightGrid = uigridlayout(app.RightPanel);
            app.RightGrid.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.RightGrid.RowHeight = {25, '3.92x', 30, 25, 25};
            app.RightGrid.Padding = [15 10 15 10];

            % Create UIAxes
            app.UIAxes = uiaxes(app.RightGrid);
            xlabel(app.UIAxes, 'X_1')
            ylabel(app.UIAxes, 'X_2')
            zlabel(app.UIAxes, 'X_3')
            app.UIAxes.XLim = [-2 4];
            app.UIAxes.YLim = [-1 5];
            app.UIAxes.ZLim = [-1 1];
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.ZGrid = 'on';
            app.UIAxes.NextPlot = 'add';
            app.UIAxes.Layout.Row = 2;
            app.UIAxes.Layout.Column = [1 7];
            app.UIAxes.ButtonDownFcn = createCallbackFcn(app, @UIAxesButtonDown, true);

            % Create SystemCheckBox
            app.SystemCheckBox = uicheckbox(app.RightGrid);
            app.SystemCheckBox.ValueChangedFcn = createCallbackFcn(app, @SystemCheckBoxValueChanged, true);
            app.SystemCheckBox.Text = 'System';
            app.SystemCheckBox.Layout.Row = 5;
            app.SystemCheckBox.Layout.Column = [3 4];
            app.SystemCheckBox.Value = true;

            % Create PersonalCheckBox
            app.PersonalCheckBox = uicheckbox(app.RightGrid);
            app.PersonalCheckBox.ValueChangedFcn = createCallbackFcn(app, @PersonalCheckBoxValueChanged, true);
            app.PersonalCheckBox.Text = 'Personal';
            app.PersonalCheckBox.Layout.Row = 5;
            app.PersonalCheckBox.Layout.Column = [5 6];

            % Create RightLabel2
            app.RightLabel2 = uilabel(app.RightGrid);
            app.RightLabel2.HorizontalAlignment = 'center';
            app.RightLabel2.Layout.Row = 4;
            app.RightLabel2.Layout.Column = [1 7];
            app.RightLabel2.Text = 'Choose which penalty/barrier functions to use';

            % Create RightLabel1
            app.RightLabel1 = uilabel(app.RightGrid);
            app.RightLabel1.HorizontalAlignment = 'center';
            app.RightLabel1.Layout.Row = 1;
            app.RightLabel1.Layout.Column = [1 7];
            app.RightLabel1.Text = 'Right-click and choose "Restart" to restart the program!';

            % Create x_0Label
            app.x_0Label = uilabel(app.RightGrid);
            app.x_0Label.Interpreter = 'latex';
            app.x_0Label.HorizontalAlignment = 'right';
            app.x_0Label.Visible = 'off';
            app.x_0Label.Layout.Row = 3;
            app.x_0Label.Layout.Column = 1;
            app.x_0Label.Text = '$x_1^{(0)} =$';

            % Create x_0
            app.x_0 = uieditfield(app.RightGrid, 'numeric');
            app.x_0.ValueChangedFcn = createCallbackFcn(app, @RestartMenuSelected, true);
            app.x_0.Visible = 'off';
            app.x_0.Layout.Row = 3;
            app.x_0.Layout.Column = 2;

            % Create y_0Label
            app.y_0Label = uilabel(app.RightGrid);
            app.y_0Label.Interpreter = 'latex';
            app.y_0Label.HorizontalAlignment = 'right';
            app.y_0Label.Visible = 'off';
            app.y_0Label.Layout.Row = 3;
            app.y_0Label.Layout.Column = 3;
            app.y_0Label.Text = '$x_2^{(0)} =$';

            % Create y_0
            app.y_0 = uieditfield(app.RightGrid, 'numeric');
            app.y_0.ValueChangedFcn = createCallbackFcn(app, @RestartMenuSelected, true);
            app.y_0.Visible = 'off';
            app.y_0.Layout.Row = 3;
            app.y_0.Layout.Column = 4;

            % Create z_0Label
            app.z_0Label = uilabel(app.RightGrid);
            app.z_0Label.Interpreter = 'latex';
            app.z_0Label.HorizontalAlignment = 'right';
            app.z_0Label.Visible = 'off';
            app.z_0Label.Layout.Row = 3;
            app.z_0Label.Layout.Column = 5;
            app.z_0Label.Text = '$x_3^{(0)} = $';

            % Create z_0
            app.z_0 = uieditfield(app.RightGrid, 'numeric');
            app.z_0.ValueChangedFcn = createCallbackFcn(app, @RestartMenuSelected, true);
            app.z_0.Visible = 'off';
            app.z_0.Layout.Row = 3;
            app.z_0.Layout.Column = 6;

            % Create IterateButton
            app.IterateButton = uibutton(app.RightGrid, 'push');
            app.IterateButton.ButtonPushedFcn = createCallbackFcn(app, @IterateButtonPushed, true);
            app.IterateButton.Visible = 'off';
            app.IterateButton.Layout.Row = 3;
            app.IterateButton.Layout.Column = 7;
            app.IterateButton.Text = 'Iterate';

            % Create ContextMenu
            app.ContextMenu = uicontextmenu(app.UIFigure);

            % Create RestartMenu
            app.RestartMenu = uimenu(app.ContextMenu);
            app.RestartMenu.MenuSelectedFcn = createCallbackFcn(app, @RestartMenuSelected, true);
            app.RestartMenu.Text = 'Restart';
            
            % Assign app.ContextMenu
            app.MuButton.ContextMenu = app.ContextMenu;
            app.EpsilonButton.ContextMenu = app.ContextMenu;
            app.RightGrid.ContextMenu = app.ContextMenu;
            app.UIAxes.ContextMenu = app.ContextMenu;

            % Create ContextMenu2
            app.ContextMenu2 = uicontextmenu(app.UIFigure);

            % Create ShowcontoursonoffMenu
            app.ShowcontoursonoffMenu = uimenu(app.ContextMenu2);
            app.ShowcontoursonoffMenu.MenuSelectedFcn = createCallbackFcn(app, @ShowContoursMenuSelected, true);
            app.ShowcontoursonoffMenu.Text = 'Show contours (on/off)';
            
            % Assign app.ContextMenu2
            app.TheoryLabel.ContextMenu = app.ContextMenu2;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Lab2

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