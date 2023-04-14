classdef main < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        SignalparametersPanel           matlab.ui.container.Panel
        Lamp_signal                     matlab.ui.control.Lamp
        samplingfrenquencyfactorEditField  matlab.ui.control.NumericEditField
        samplingfrenquencyfactorEditFieldLabel  matlab.ui.control.Label
        SNR_textnormaloutSlider         matlab.ui.control.Slider
        SNR_textnormaloutSliderLabel    matlab.ui.control.Label
        HistoryTextArea                 matlab.ui.control.TextArea
        Signal0directLOSPanel           matlab.ui.container.Panel
        degLabel_2                      matlab.ui.control.Label
        HzLabel_2                       matlab.ui.control.Label
        sampleLabel_2                   matlab.ui.control.Label
        phasephi_0EditField             matlab.ui.control.NumericEditField
        phasephi_0EditField_2Label      matlab.ui.control.Label
        amplituderho_0EditField         matlab.ui.control.NumericEditField
        amplituderho_0EditField_2Label  matlab.ui.control.Label
        Dopplerfrequencyf_d0EditField   matlab.ui.control.NumericEditField
        Dopplerfrequencyf_d0EditField_3Label  matlab.ui.control.Label
        timedelaytau_0EditField         matlab.ui.control.NumericEditField
        timedelaytau_0EditField_2Label  matlab.ui.control.Label
        Signal1reflectedNLOSPanel       matlab.ui.container.Panel
        degLabel                        matlab.ui.control.Label
        HzLabel                         matlab.ui.control.Label
        sampleLabel                     matlab.ui.control.Label
        phasephi_1EditField             matlab.ui.control.NumericEditField
        phasephi_1EditField_2Label      matlab.ui.control.Label
        amplituderho_1EditField         matlab.ui.control.NumericEditField
        amplituderho_1EditField_2Label  matlab.ui.control.Label
        Dopplerfrequencyf_d1EditField   matlab.ui.control.NumericEditField
        Dopplerfrequencyf_d1EditField_2Label  matlab.ui.control.Label
        timedelaytau_1EditField         matlab.ui.control.NumericEditField
        timedelaytau_1EditField_2Label  matlab.ui.control.Label
        LoadButton                      matlab.ui.control.Button
        basebandsignalsamplesEditField  matlab.ui.control.EditField
        basebandsignalsamplesEditFieldLabel  matlab.ui.control.Label
        DualSourceCramrRaoBound2SCRBPanel  matlab.ui.container.Panel
        Lamp_2SCRB                      matlab.ui.control.Lamp
        CRB2Sphi_1EditField             matlab.ui.control.NumericEditField
        sqrttextnormalCRBphi_1Label     matlab.ui.control.Label
        CRB2Srho_1EditField             matlab.ui.control.NumericEditField
        sqrttextnormalCRBrho_1Label     matlab.ui.control.Label
        CRB2Sf_d1EditField              matlab.ui.control.NumericEditField
        sqrttextnormalCRBf_d1Label      matlab.ui.control.Label
        CRB2Stau_1EditField             matlab.ui.control.NumericEditField
        sqrttextnormalCRBtau_1Label     matlab.ui.control.Label
        degLabel_7                      matlab.ui.control.Label
        HzLabel_7                       matlab.ui.control.Label
        sampleLabel_7                   matlab.ui.control.Label
        CRB2Sphi_0EditField             matlab.ui.control.NumericEditField
        sqrttextnormalCRBphi_0Label     matlab.ui.control.Label
        CRB2Srho_0EditField             matlab.ui.control.NumericEditField
        sqrttextnormalCRBrho_0Label     matlab.ui.control.Label
        CRB2Sf_d0EditField              matlab.ui.control.NumericEditField
        sqrttextnormalCRBf_d0Label      matlab.ui.control.Label
        CRB2Stau_0EditField             matlab.ui.control.NumericEditField
        sqrttextnormalCRBtau_0Label     matlab.ui.control.Label
        degLabel_6                      matlab.ui.control.Label
        HzLabel_6                       matlab.ui.control.Label
        sampleLabel_6                   matlab.ui.control.Label
        ProcessButton_3                 matlab.ui.control.Button
        CleantoCompositeBoundRatioCCBRPanel  matlab.ui.container.Panel
        Lamp_CCBR                       matlab.ui.control.Lamp
        CCBRf_dEditField                matlab.ui.control.NumericEditField
        CCBRf_dEditFieldLabel           matlab.ui.control.Label
        CCBRtauEditField                matlab.ui.control.NumericEditField
        CCBRtauEditFieldLabel           matlab.ui.control.Label
        ProcessButton_4                 matlab.ui.control.Button
        MisspecifiedCramrRaoBoundMCRBPanel  matlab.ui.container.Panel
        Lamp_MCRB                       matlab.ui.control.Lamp
        sqrttextnormalMCRBphi_0EditField  matlab.ui.control.NumericEditField
        sqrttextnormalMCRBphi_0EditFieldLabel  matlab.ui.control.Label
        sqrttextnormalMCRBrho_0EditField  matlab.ui.control.NumericEditField
        sqrttextnormalMCRBrho_0EditFieldLabel  matlab.ui.control.Label
        sqrttextnormalMCRBf_d0EditField  matlab.ui.control.NumericEditField
        sqrttextnormalMCRBf_d0Label     matlab.ui.control.Label
        sqrttextnormalMCRBtau_0EditField  matlab.ui.control.NumericEditField
        sqrttextnormalMCRBtau_0EditFieldLabel  matlab.ui.control.Label
        degLabel_5                      matlab.ui.control.Label
        HzLabel_5                       matlab.ui.control.Label
        sampleLabel_5                   matlab.ui.control.Label
        biasphi_0EditField              matlab.ui.control.NumericEditField
        biasphi_0EditFieldLabel         matlab.ui.control.Label
        biasrho_0EditField              matlab.ui.control.NumericEditField
        biasrho_0EditFieldLabel         matlab.ui.control.Label
        biasf_d0EditField               matlab.ui.control.NumericEditField
        biasf_d0EditFieldLabel          matlab.ui.control.Label
        biastau_0EditField              matlab.ui.control.NumericEditField
        biastau_0EditFieldLabel         matlab.ui.control.Label
        degLabel_4                      matlab.ui.control.Label
        HzLabel_4                       matlab.ui.control.Label
        sampleLabel_4                   matlab.ui.control.Label
        ProcessButton_2                 matlab.ui.control.Button
        SingleSourceCramrRaoBound1SCRBPanel  matlab.ui.container.Panel
        Lamp_1SCRB                      matlab.ui.control.Lamp
        CRB1SphiEditField               matlab.ui.control.NumericEditField
        sqrttextnormalCRBphiLabel       matlab.ui.control.Label
        CRB1SrhoEditField               matlab.ui.control.NumericEditField
        sqrttextnormalCRBrhoLabel       matlab.ui.control.Label
        CRB1Sf_dEditField               matlab.ui.control.NumericEditField
        sqrttextnormalCRBf_dLabel       matlab.ui.control.Label
        CRB1StauEditField               matlab.ui.control.NumericEditField
        sqrttextnormalCRBtauLabel       matlab.ui.control.Label
        degLabel_3                      matlab.ui.control.Label
        HzLabel_3                       matlab.ui.control.Label
        sampleLabel_3                   matlab.ui.control.Label
        ProcessButton                   matlab.ui.control.Button
    end

    
    methods (Access = private)
        
        function updateCommandWindow(app,newString)
            oldString = app.HistoryTextArea.Value; % The string as it is now.
            app.HistoryTextArea.Value = [{newString};oldString];
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            warning('off');
            addpath('src/');
            global considered_signal 
            global fs fc Fc F0 Tcode
            
            considered_signal   = [];
            fs                  = [];
            fc                  = [];
            Fc                  = [];
            F0                  = [];
            Tcode               = [];
            
            global lampColor
            lampColor = [0.39, 0.83, 0.07;... % lampColor(1,:) = green light
                         0.85, 0.33, 0.10;... % lampColor(2,:) = orange light
                         0.85, 0.09, 0.09];   % lampColor(3,:) = red light
            
        end

        % Button pushed function: LoadButton
        function LoadButtonPushed(app, event)
            global considered_signal 
            global fs fc Fc F0 Tcode
            global lampColor
            
            [fileName,pathName] = uigetfile('*.mat','Choose the file ...'); % load the file
            if pathName
                load([pathName fileName]);
                
                F0 = 1023000;
                Fc = 1575.42e6;
                fc = Fc/F0;
    
                original_signal = c_0;
                Tcode = 1e-3;
                
                fs = app.samplingfrenquencyfactorEditField.Value;
                
                considered_signal = repmat(original_signal',fs,1);
                considered_signal = considered_signal(:);
                
                app.basebandsignalsamplesEditField.Value = fileName;
                updateCommandWindow(app, [ fileName,' successfully loaded!']);
                
                app.Lamp_1SCRB.Color = lampColor(3,:);
                app.Lamp_2SCRB.Color = lampColor(3,:);
                app.Lamp_MCRB.Color = lampColor(3,:);
                app.Lamp_CCBR.Color = lampColor(3,:);
                
                app.Lamp_signal.Color = lampColor(1,:);
            else
                updateCommandWindow(app, 'error: no file loaded.');
                
                app.Lamp_signal.Color = lampColor(3,:);
            end
        end

        % Button pushed function: ProcessButton
        function ProcessButtonPushed(app, event)
            global considered_signal 
            global fs fc Fc
            global lampColor
            
            % Process '1S-CRB'
            updateCommandWindow(app, 'Processing 1S-CRB...');
            tic;
            
            % 1- retrieve LOS parameters (assumption: no NLOS)
            tau_0   = app.timedelaytau_0EditField.Value;
            fd_0    = app.Dopplerfrequencyf_d0EditField.Value;
            rho_0   = app.amplituderho_0EditField.Value;
            phi_0   = app.phasephi_0EditField.Value*pi/180; % rad
            
            epsilon = [tau_0; fd_0/Fc; rho_0; phi_0; zeros(4,1)];
            
            SNR_out = 10^(app.SNR_textnormaloutSlider.Value/20); % linear
            attenuationSquared = SNR_out/(sum(abs(considered_signal).^2)*rho_0^2);
            
            % 2- compute Fisher Information Matrix
            FIM = 2*fs*real(generate_FIM(considered_signal, fs, fc, epsilon));

            % 3- FIM inversion leads directly to the Cramer-Rao Bound
            CRB = real(inv(FIM(1:4,1:4)))/attenuationSquared; % the first four terms are required for 1S scenario
            
            CRB_tau_0   = sqrt(CRB(1,1));
            CRB_fd_0    = sqrt(CRB(2,2))*Fc;
            CRB_rho_0   = sqrt(CRB(3,3));
            CRB_phi_0   = sqrt(CRB(4,4))*180/pi;
            
            % 4- display results
            app.CRB1StauEditField.Value = CRB_tau_0;
            app.CRB1Sf_dEditField.Value = CRB_fd_0;
            app.CRB1SrhoEditField.Value = CRB_rho_0;
            app.CRB1SphiEditField.Value = CRB_phi_0;
            
            t=toc;
            updateCommandWindow(app, ['1S-CRB processed (' num2str(t) ' s)']);
            
            app.Lamp_1SCRB.Color = lampColor(1,:);
        end

        % Button pushed function: ProcessButton_3
        function ProcessButton_3Pushed(app, event)
            global considered_signal 
            global fs fc Fc
            global lampColor
            
            % Process '2S-CRB'
            updateCommandWindow(app, 'Processing 2S-CRB...');
            tic;
            
            % 1- retrieve LOS and NLOS parameters
            tau_0   = app.timedelaytau_0EditField.Value;
            fd_0    = app.Dopplerfrequencyf_d0EditField.Value;
            rho_0   = app.amplituderho_0EditField.Value;
            phi_0   = app.phasephi_0EditField.Value*pi/180; % rad
            
            tau_1   = app.timedelaytau_1EditField.Value;
            fd_1    = app.Dopplerfrequencyf_d1EditField.Value; 
            rho_1   = app.amplituderho_1EditField.Value;
            phi_1   = app.phasephi_1EditField.Value*pi/180; % rad
            
            epsilon = [tau_0; fd_0/Fc; rho_0; phi_0; tau_1; fd_1/Fc; rho_1; phi_1];
            
            SNR_out = 10^(app.SNR_textnormaloutSlider.Value/20); % linear
            attenuationSquared = SNR_out/(sum(abs(considered_signal).^2)*rho_0^2);
            
            % 2- compute Fisher Information Matrix
            FIM = 2*fs*real(generate_FIM(considered_signal, fs, fc, epsilon));

            % 3- FIM inversion leads directly to the Cramer-Rao Bound
            CRB = real(inv(FIM))/attenuationSquared;
            
            CRB_tau_0   = sqrt(CRB(1,1));
            CRB_fd_0    = sqrt(CRB(2,2))*Fc;
            CRB_rho_0   = sqrt(CRB(3,3));
            CRB_phi_0   = sqrt(CRB(4,4))*180/pi;

            CRB_tau_1   = sqrt(CRB(5,5));
            CRB_fd_1    = sqrt(CRB(6,6))*Fc;
            CRB_rho_1   = sqrt(CRB(7,7));
            CRB_phi_1   = sqrt(CRB(8,8))*180/pi;

            % 4- display results
            app.CRB2Stau_0EditField.Value = CRB_tau_0;
            app.CRB2Sf_d0EditField.Value = CRB_fd_0;
            app.CRB2Srho_0EditField.Value = CRB_rho_0;
            app.CRB2Sphi_0EditField.Value = CRB_phi_0;

            app.CRB2Stau_1EditField.Value = CRB_tau_1;
            app.CRB2Sf_d1EditField.Value = CRB_fd_1;
            app.CRB2Srho_1EditField.Value = CRB_rho_1;
            app.CRB2Sphi_1EditField.Value = CRB_phi_1;

            t=toc;
            updateCommandWindow(app, ['2S-CRB processed (' num2str(t) ' s)']);
            
            app.Lamp_2SCRB.Color = lampColor(1,:);
        end

        % Button pushed function: ProcessButton_2
        function ProcessButton_2Pushed(app, event)
            global considered_signal 
            global fs fc Fc F0 Tcode
            global lampColor
            OSF = 100;
            
            % Process 'MCRB'
            updateCommandWindow(app, 'Processing MCRB...');
            tic;
            
            % 1- retrieve LOS and NLOS parameters
            tau_0   = app.timedelaytau_0EditField.Value;
            fd_0    = app.Dopplerfrequencyf_d0EditField.Value;
            rho_0   = app.amplituderho_0EditField.Value;
            phi_0   = app.phasephi_0EditField.Value*pi/180; % rad
            
            tau_1   = app.timedelaytau_1EditField.Value;
            fd_1    = app.Dopplerfrequencyf_d1EditField.Value;
            rho_1   = app.amplituderho_1EditField.Value;
            phi_1   = app.phasephi_1EditField.Value*pi/180; % rad
            
            epsilon = [tau_0; fd_0/Fc; rho_0; phi_0; tau_1; fd_1/Fc; rho_1; phi_1];
            
            SNR_out = 10^(app.SNR_textnormaloutSlider.Value/20); % linear
            attenuationSquared = SNR_out/(sum(abs(considered_signal).^2)*rho_0^2);
            
            % 2- compute pseudo-true parameters (biases) using 1S-MLE procedure
            dataSamp         = adapt_signal(considered_signal, OSF);

            time_axis = (0:(length(considered_signal)-1))*Tcode/length(considered_signal);
            Fd_axis   = fd_0 + linspace(-abs(fd_1-fd_0),abs(fd_1-fd_0),201);
            
            doppler_table = exp(-1j*2*pi*Fd_axis.*time_axis'); % first axis is the doppler axis, second is time
            
            c_multipath = rho_0*exp(1j*phi_0)*circshift(considered_signal,tau_0).*exp(-1j*2*pi*fd_0*time_axis') ...
            + rho_1*exp(1j*phi_1)*circshift(considered_signal,tau_1).*exp(-1j*2*pi*fd_1*(time_axis' - (tau_1-tau_0)/(fs*F0)));

            % 2a- Doppler frequency estimation and compensation
            [~, fd_pt, ~] = Single_Source_MLE_tau_b(c_multipath, considered_signal, doppler_table, Fd_axis, 1);
            
            c_multipath_noDoppler = c_multipath.*exp(1j*2*pi*fd_pt.*time_axis');
            
            dataSamp_multipath = adapt_signal(c_multipath_noDoppler,OSF);
            
            % 2b- time delay and complex amplitude estimation
            [tau_pt, ~, rho_pt, phi_pt] = Single_Source_MLE_tau_b(dataSamp_multipath, dataSamp, ones(size(dataSamp)), 0, 1);
            tau_pt = tau_pt/OSF;
            
            % 3- compute information matrices
            epsilon_pt = [tau_pt; fd_pt/Fc; rho_pt; phi_pt];
            
            % 3a- matrix B = FIM_pseudo_true
            B = 2*fs*real(generate_FIM(considered_signal, fs, fc, [epsilon_pt;zeros(4,1)]));
            B = B(1:4,1:4);
            
            % 3b- matrix A: contribution of the bias
            A = 2*fs*real(generate_MCRB_term([epsilon;epsilon_pt], considered_signal, fs, fc));
            A = A-B;
            
            % 4- compute MCRB
            MCRB = real(inv(A))*B*real(inv(A))/attenuationSquared;
            
            MCRB_tau = sqrt(MCRB(1,1));
            MCRB_fd   = sqrt(MCRB(2,2))*Fc;
            MCRB_rho = sqrt(MCRB(3,3));
            MCRB_phi = sqrt(MCRB(4,4))*180/pi;
            
            % 5- display results
            app.biastau_0EditField.Value = tau_pt - tau_0;
            app.biasf_d0EditField.Value  = fd_pt - fd_0;
            app.biasrho_0EditField.Value = rho_pt - rho_0;
            app.biasphi_0EditField.Value = (phi_pt - phi_0)*180/pi;

            app.sqrttextnormalMCRBtau_0EditField.Value = MCRB_tau;
            app.sqrttextnormalMCRBf_d0EditField.Value  = MCRB_fd;
            app.sqrttextnormalMCRBrho_0EditField.Value = MCRB_rho;
            app.sqrttextnormalMCRBphi_0EditField.Value = MCRB_phi;

            t=toc;
            updateCommandWindow(app, ['MCRB processed (' num2str(t) ' s)']);
            
            app.Lamp_MCRB.Color = lampColor(1,:);
        end

        % Button pushed function: ProcessButton_4
        function ProcessButton_4Pushed(app, event)
            global considered_signal 
            global fs fc Fc
            global lampColor
            
            % Process 'CCBR'
            updateCommandWindow(app, 'Processing CCBR...');
            tic;
            
            % 1- retrieve LOS and NLOS parameters
            tau_0   = app.timedelaytau_0EditField.Value;
            fd_0    = app.Dopplerfrequencyf_d0EditField.Value;
            rho_0   = app.amplituderho_0EditField.Value;
            phi_0   = app.phasephi_0EditField.Value*pi/180; % rad
            
            tau_1   = app.timedelaytau_1EditField.Value;
            fd_1    = app.Dopplerfrequencyf_d1EditField.Value; 
            rho_1   = app.amplituderho_1EditField.Value;
            phi_1   = app.phasephi_1EditField.Value*pi/180; % rad
            
            delta_phi = phi_1 - phi_0;
            
            epsilon = [tau_0; fd_0/Fc; rho_0; phi_0; tau_1; fd_1/Fc; rho_1; phi_1];
            
            % 2- compute CCBR
            
            % 2a- composite bound
            PHI = 2*fs*generate_Phi(considered_signal,fs,fc,epsilon);

            PHI1 = PHI(1:2,1:2);
            PHI2 = PHI(3:4,3:4);
            PHI21 = PHI(3:4,1:2);

            A = real(PHI1);
            B = 0.5*(real(PHI21)/real(PHI2)*real(PHI21) + imag(PHI21)/real(PHI2)*imag(PHI21));
            C = 0.5*(real(PHI21)/real(PHI2)*real(PHI21) - imag(PHI21)/real(PHI2)*imag(PHI21));
            D = 0.5*(real(PHI21)/real(PHI2)*imag(PHI21) + imag(PHI21)/real(PHI2)*real(PHI21));

    
            det = (A(1,1) - B(1,1) - cos(2*delta_phi)*C(1,1) + sin(2*delta_phi)*D(1,1)) ...
                 .*(A(2,2) - B(2,2) - cos(2*delta_phi)*C(2,2) + sin(2*delta_phi)*D(2,2)) ...
                 -(A(1,2) - B(1,2) - cos(2*delta_phi)*C(1,2) + sin(2*delta_phi)*D(1,2)) ...
                 .*(A(2,1) - B(2,1) - cos(2*delta_phi)*C(2,1) + sin(2*delta_phi)*D(2,1));

            CRB_tau_0_epsilon = (A(2,2) - B(2,2) - cos(2*delta_phi)*C(2,2) + sin(2*delta_phi)*D(2,2))./det;
            CRB_b_0_epsilon = (A(1,1) - B(1,1) - cos(2*delta_phi)*C(1,1) + sin(2*delta_phi)*D(1,1))./det;
            
            % b- clean bound
            FIM = 2*fs*real(generate_FIM(considered_signal,fs,fc,epsilon));
            CRB00 = real(inv(FIM(1:4,1:4)));
            
            % c- evaluation of the CCBR
            CCBR_tau = CRB00(1,1)/CRB_tau_0_epsilon;
            CCBR_fd  = CRB00(2,2)/CRB_b_0_epsilon;

             %3- display results
            app.CCBRtauEditField.Value  = CCBR_tau;
            app.CCBRf_dEditField.Value  = CCBR_fd;
            
            t=toc;
            updateCommandWindow(app, ['CCBR processed (' num2str(t) ' s)']);
            
            app.Lamp_CCBR.Color = lampColor(1,:);
        end

        % Value changed function: Dopplerfrequencyf_d1EditField, 
        % amplituderho_1EditField, phasephi_1EditField, 
        % timedelaytau_1EditField
        function parameterSignal_1EditFieldValueChanged(app, event)
            global lampColor
            app.Lamp_2SCRB.Color = lampColor(3,:);
            app.Lamp_MCRB.Color = lampColor(3,:);
            app.Lamp_CCBR.Color = lampColor(3,:);
        end

        % Value changed function: timedelaytau_0EditField
        function parameterSignal_0EditFieldValueChanged(app, event)
            global lampColor
            app.Lamp_1SCRB.Color = lampColor(3,:);
            app.Lamp_2SCRB.Color = lampColor(3,:);
            app.Lamp_MCRB.Color = lampColor(3,:);
            app.Lamp_CCBR.Color = lampColor(3,:);
        end

        % Value changed function: samplingfrenquencyfactorEditField
        function samplingfrenquencyfactorEditFieldValueChanged(app, event)
            global lampColor
            % if a signal was loaded, then it should be refreshed to be
            % adapted to the new value of fs.
            if sum(app.Lamp_signal.Color) == sum(lampColor(1,:))
                app.Lamp_signal.Color = lampColor(2,:);
            end
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            clear global considered_signal fs fc Fc F0 Tcode lampColor
            delete(app)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1150 777];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create SingleSourceCramrRaoBound1SCRBPanel
            app.SingleSourceCramrRaoBound1SCRBPanel = uipanel(app.UIFigure);
            app.SingleSourceCramrRaoBound1SCRBPanel.TitlePosition = 'centertop';
            app.SingleSourceCramrRaoBound1SCRBPanel.Title = 'Single Source Cramér-Rao Bound (1S-CRB)';
            app.SingleSourceCramrRaoBound1SCRBPanel.FontName = 'Calibri';
            app.SingleSourceCramrRaoBound1SCRBPanel.FontWeight = 'bold';
            app.SingleSourceCramrRaoBound1SCRBPanel.FontSize = 16;
            app.SingleSourceCramrRaoBound1SCRBPanel.Position = [38 259 530 228];

            % Create ProcessButton
            app.ProcessButton = uibutton(app.SingleSourceCramrRaoBound1SCRBPanel, 'push');
            app.ProcessButton.ButtonPushedFcn = createCallbackFcn(app, @ProcessButtonPushed, true);
            app.ProcessButton.FontName = 'Calibri';
            app.ProcessButton.Position = [450 27 61 22];
            app.ProcessButton.Text = 'Process';

            % Create sampleLabel_3
            app.sampleLabel_3 = uilabel(app.SingleSourceCramrRaoBound1SCRBPanel);
            app.sampleLabel_3.Interpreter = 'latex';
            app.sampleLabel_3.Position = [178 157 55 22];
            app.sampleLabel_3.Text = '[sample]';

            % Create HzLabel_3
            app.HzLabel_3 = uilabel(app.SingleSourceCramrRaoBound1SCRBPanel);
            app.HzLabel_3.Interpreter = 'latex';
            app.HzLabel_3.Position = [178 114 55 22];
            app.HzLabel_3.Text = '[Hz]';

            % Create degLabel_3
            app.degLabel_3 = uilabel(app.SingleSourceCramrRaoBound1SCRBPanel);
            app.degLabel_3.Interpreter = 'latex';
            app.degLabel_3.Position = [178 28 55 22];
            app.degLabel_3.Text = '[deg]';

            % Create sqrttextnormalCRBtauLabel
            app.sqrttextnormalCRBtauLabel = uilabel(app.SingleSourceCramrRaoBound1SCRBPanel);
            app.sqrttextnormalCRBtauLabel.Interpreter = 'latex';
            app.sqrttextnormalCRBtauLabel.Position = [21 156 56 22];
            app.sqrttextnormalCRBtauLabel.Text = '$\sqrt{\textnormal{CRB}}\tau$';

            % Create CRB1StauEditField
            app.CRB1StauEditField = uieditfield(app.SingleSourceCramrRaoBound1SCRBPanel, 'numeric');
            app.CRB1StauEditField.ValueDisplayFormat = '%.4g';
            app.CRB1StauEditField.Editable = 'off';
            app.CRB1StauEditField.Position = [95 152 77 30];

            % Create sqrttextnormalCRBf_dLabel
            app.sqrttextnormalCRBf_dLabel = uilabel(app.SingleSourceCramrRaoBound1SCRBPanel);
            app.sqrttextnormalCRBf_dLabel.Interpreter = 'latex';
            app.sqrttextnormalCRBf_dLabel.Position = [21 113 59 22];
            app.sqrttextnormalCRBf_dLabel.Text = '$\sqrt{\textnormal{CRB}} f_d$';

            % Create CRB1Sf_dEditField
            app.CRB1Sf_dEditField = uieditfield(app.SingleSourceCramrRaoBound1SCRBPanel, 'numeric');
            app.CRB1Sf_dEditField.ValueDisplayFormat = '%.4g';
            app.CRB1Sf_dEditField.Editable = 'off';
            app.CRB1Sf_dEditField.Position = [95 109 77 30];

            % Create sqrttextnormalCRBrhoLabel
            app.sqrttextnormalCRBrhoLabel = uilabel(app.SingleSourceCramrRaoBound1SCRBPanel);
            app.sqrttextnormalCRBrhoLabel.Interpreter = 'latex';
            app.sqrttextnormalCRBrhoLabel.Position = [21 70 58 22];
            app.sqrttextnormalCRBrhoLabel.Text = '$\sqrt{\textnormal{CRB}} \rho$';

            % Create CRB1SrhoEditField
            app.CRB1SrhoEditField = uieditfield(app.SingleSourceCramrRaoBound1SCRBPanel, 'numeric');
            app.CRB1SrhoEditField.ValueDisplayFormat = '%.4g';
            app.CRB1SrhoEditField.Editable = 'off';
            app.CRB1SrhoEditField.Position = [95 66 77 30];

            % Create sqrttextnormalCRBphiLabel
            app.sqrttextnormalCRBphiLabel = uilabel(app.SingleSourceCramrRaoBound1SCRBPanel);
            app.sqrttextnormalCRBphiLabel.Interpreter = 'latex';
            app.sqrttextnormalCRBphiLabel.Position = [21 27 59 22];
            app.sqrttextnormalCRBphiLabel.Text = '$\sqrt{\textnormal{CRB}} \phi$';

            % Create CRB1SphiEditField
            app.CRB1SphiEditField = uieditfield(app.SingleSourceCramrRaoBound1SCRBPanel, 'numeric');
            app.CRB1SphiEditField.ValueDisplayFormat = '%.4g';
            app.CRB1SphiEditField.Editable = 'off';
            app.CRB1SphiEditField.Position = [95 23 77 30];

            % Create Lamp_1SCRB
            app.Lamp_1SCRB = uilamp(app.SingleSourceCramrRaoBound1SCRBPanel);
            app.Lamp_1SCRB.Position = [510 186 10 10];
            app.Lamp_1SCRB.Color = [0.851 0.0902 0.0902];

            % Create MisspecifiedCramrRaoBoundMCRBPanel
            app.MisspecifiedCramrRaoBoundMCRBPanel = uipanel(app.UIFigure);
            app.MisspecifiedCramrRaoBoundMCRBPanel.TitlePosition = 'centertop';
            app.MisspecifiedCramrRaoBoundMCRBPanel.Title = 'Misspecified Cramér-Rao Bound (MCRB)';
            app.MisspecifiedCramrRaoBoundMCRBPanel.FontName = 'Calibri';
            app.MisspecifiedCramrRaoBoundMCRBPanel.FontWeight = 'bold';
            app.MisspecifiedCramrRaoBoundMCRBPanel.FontSize = 16;
            app.MisspecifiedCramrRaoBoundMCRBPanel.Position = [591 259 530 228];

            % Create ProcessButton_2
            app.ProcessButton_2 = uibutton(app.MisspecifiedCramrRaoBoundMCRBPanel, 'push');
            app.ProcessButton_2.ButtonPushedFcn = createCallbackFcn(app, @ProcessButton_2Pushed, true);
            app.ProcessButton_2.FontName = 'Calibri';
            app.ProcessButton_2.Position = [451 27 61 22];
            app.ProcessButton_2.Text = 'Process';

            % Create sampleLabel_4
            app.sampleLabel_4 = uilabel(app.MisspecifiedCramrRaoBoundMCRBPanel);
            app.sampleLabel_4.Interpreter = 'latex';
            app.sampleLabel_4.Position = [178 157 55 22];
            app.sampleLabel_4.Text = '[sample]';

            % Create HzLabel_4
            app.HzLabel_4 = uilabel(app.MisspecifiedCramrRaoBoundMCRBPanel);
            app.HzLabel_4.Interpreter = 'latex';
            app.HzLabel_4.Position = [178 114 55 22];
            app.HzLabel_4.Text = '[Hz]';

            % Create degLabel_4
            app.degLabel_4 = uilabel(app.MisspecifiedCramrRaoBoundMCRBPanel);
            app.degLabel_4.Interpreter = 'latex';
            app.degLabel_4.Position = [178 28 55 22];
            app.degLabel_4.Text = '[deg]';

            % Create biastau_0EditFieldLabel
            app.biastau_0EditFieldLabel = uilabel(app.MisspecifiedCramrRaoBoundMCRBPanel);
            app.biastau_0EditFieldLabel.Interpreter = 'latex';
            app.biastau_0EditFieldLabel.Position = [34 156 45 22];
            app.biastau_0EditFieldLabel.Text = 'bias $\tau_0$';

            % Create biastau_0EditField
            app.biastau_0EditField = uieditfield(app.MisspecifiedCramrRaoBoundMCRBPanel, 'numeric');
            app.biastau_0EditField.ValueDisplayFormat = '%.4g';
            app.biastau_0EditField.Editable = 'off';
            app.biastau_0EditField.Position = [94 152 79 30];

            % Create biasf_d0EditFieldLabel
            app.biasf_d0EditFieldLabel = uilabel(app.MisspecifiedCramrRaoBoundMCRBPanel);
            app.biasf_d0EditFieldLabel.Interpreter = 'latex';
            app.biasf_d0EditFieldLabel.Position = [34 113 51 22];
            app.biasf_d0EditFieldLabel.Text = 'bias $f_{d, 0}$';

            % Create biasf_d0EditField
            app.biasf_d0EditField = uieditfield(app.MisspecifiedCramrRaoBoundMCRBPanel, 'numeric');
            app.biasf_d0EditField.ValueDisplayFormat = '%.4g';
            app.biasf_d0EditField.Editable = 'off';
            app.biasf_d0EditField.Position = [94 109 79 30];

            % Create biasrho_0EditFieldLabel
            app.biasrho_0EditFieldLabel = uilabel(app.MisspecifiedCramrRaoBoundMCRBPanel);
            app.biasrho_0EditFieldLabel.Interpreter = 'latex';
            app.biasrho_0EditFieldLabel.Position = [34 70 46 22];
            app.biasrho_0EditFieldLabel.Text = 'bias $\rho_0$';

            % Create biasrho_0EditField
            app.biasrho_0EditField = uieditfield(app.MisspecifiedCramrRaoBoundMCRBPanel, 'numeric');
            app.biasrho_0EditField.ValueDisplayFormat = '%.4g';
            app.biasrho_0EditField.Editable = 'off';
            app.biasrho_0EditField.Position = [94 66 79 30];

            % Create biasphi_0EditFieldLabel
            app.biasphi_0EditFieldLabel = uilabel(app.MisspecifiedCramrRaoBoundMCRBPanel);
            app.biasphi_0EditFieldLabel.Interpreter = 'latex';
            app.biasphi_0EditFieldLabel.Position = [34 27 48 22];
            app.biasphi_0EditFieldLabel.Text = 'bias $\phi_0$';

            % Create biasphi_0EditField
            app.biasphi_0EditField = uieditfield(app.MisspecifiedCramrRaoBoundMCRBPanel, 'numeric');
            app.biasphi_0EditField.ValueDisplayFormat = '%.4g';
            app.biasphi_0EditField.Editable = 'off';
            app.biasphi_0EditField.Position = [94 23 79 30];

            % Create sampleLabel_5
            app.sampleLabel_5 = uilabel(app.MisspecifiedCramrRaoBoundMCRBPanel);
            app.sampleLabel_5.Interpreter = 'latex';
            app.sampleLabel_5.Position = [408 157 55 22];
            app.sampleLabel_5.Text = '[sample]';

            % Create HzLabel_5
            app.HzLabel_5 = uilabel(app.MisspecifiedCramrRaoBoundMCRBPanel);
            app.HzLabel_5.Interpreter = 'latex';
            app.HzLabel_5.Position = [408 114 55 22];
            app.HzLabel_5.Text = '[Hz]';

            % Create degLabel_5
            app.degLabel_5 = uilabel(app.MisspecifiedCramrRaoBoundMCRBPanel);
            app.degLabel_5.Interpreter = 'latex';
            app.degLabel_5.Position = [408 28 55 22];
            app.degLabel_5.Text = '[deg]';

            % Create sqrttextnormalMCRBtau_0EditFieldLabel
            app.sqrttextnormalMCRBtau_0EditFieldLabel = uilabel(app.MisspecifiedCramrRaoBoundMCRBPanel);
            app.sqrttextnormalMCRBtau_0EditFieldLabel.Interpreter = 'latex';
            app.sqrttextnormalMCRBtau_0EditFieldLabel.Position = [232 156 73 22];
            app.sqrttextnormalMCRBtau_0EditFieldLabel.Text = '$\sqrt{\textnormal{MCRB}} \tau_0$';

            % Create sqrttextnormalMCRBtau_0EditField
            app.sqrttextnormalMCRBtau_0EditField = uieditfield(app.MisspecifiedCramrRaoBoundMCRBPanel, 'numeric');
            app.sqrttextnormalMCRBtau_0EditField.ValueDisplayFormat = '%.4g';
            app.sqrttextnormalMCRBtau_0EditField.Editable = 'off';
            app.sqrttextnormalMCRBtau_0EditField.Position = [323 152 79 30];

            % Create sqrttextnormalMCRBf_d0Label
            app.sqrttextnormalMCRBf_d0Label = uilabel(app.MisspecifiedCramrRaoBoundMCRBPanel);
            app.sqrttextnormalMCRBf_d0Label.Interpreter = 'latex';
            app.sqrttextnormalMCRBf_d0Label.Position = [232 113 79 22];
            app.sqrttextnormalMCRBf_d0Label.Text = '$\sqrt{\textnormal{MCRB}} f_{d, 0}$';

            % Create sqrttextnormalMCRBf_d0EditField
            app.sqrttextnormalMCRBf_d0EditField = uieditfield(app.MisspecifiedCramrRaoBoundMCRBPanel, 'numeric');
            app.sqrttextnormalMCRBf_d0EditField.ValueDisplayFormat = '%.4g';
            app.sqrttextnormalMCRBf_d0EditField.Editable = 'off';
            app.sqrttextnormalMCRBf_d0EditField.Position = [323 109 79 30];

            % Create sqrttextnormalMCRBrho_0EditFieldLabel
            app.sqrttextnormalMCRBrho_0EditFieldLabel = uilabel(app.MisspecifiedCramrRaoBoundMCRBPanel);
            app.sqrttextnormalMCRBrho_0EditFieldLabel.Interpreter = 'latex';
            app.sqrttextnormalMCRBrho_0EditFieldLabel.Position = [232 70 74 22];
            app.sqrttextnormalMCRBrho_0EditFieldLabel.Text = '$\sqrt{\textnormal{MCRB}} \rho_0$';

            % Create sqrttextnormalMCRBrho_0EditField
            app.sqrttextnormalMCRBrho_0EditField = uieditfield(app.MisspecifiedCramrRaoBoundMCRBPanel, 'numeric');
            app.sqrttextnormalMCRBrho_0EditField.ValueDisplayFormat = '%.4g';
            app.sqrttextnormalMCRBrho_0EditField.Editable = 'off';
            app.sqrttextnormalMCRBrho_0EditField.Position = [323 66 79 30];

            % Create sqrttextnormalMCRBphi_0EditFieldLabel
            app.sqrttextnormalMCRBphi_0EditFieldLabel = uilabel(app.MisspecifiedCramrRaoBoundMCRBPanel);
            app.sqrttextnormalMCRBphi_0EditFieldLabel.Interpreter = 'latex';
            app.sqrttextnormalMCRBphi_0EditFieldLabel.Position = [232 27 76 22];
            app.sqrttextnormalMCRBphi_0EditFieldLabel.Text = '$\sqrt{\textnormal{MCRB}} \phi_0$';

            % Create sqrttextnormalMCRBphi_0EditField
            app.sqrttextnormalMCRBphi_0EditField = uieditfield(app.MisspecifiedCramrRaoBoundMCRBPanel, 'numeric');
            app.sqrttextnormalMCRBphi_0EditField.ValueDisplayFormat = '%.4g';
            app.sqrttextnormalMCRBphi_0EditField.Editable = 'off';
            app.sqrttextnormalMCRBphi_0EditField.Position = [323 23 79 30];

            % Create Lamp_MCRB
            app.Lamp_MCRB = uilamp(app.MisspecifiedCramrRaoBoundMCRBPanel);
            app.Lamp_MCRB.Position = [511 186 10 10];
            app.Lamp_MCRB.Color = [0.851 0.0902 0.0902];

            % Create CleantoCompositeBoundRatioCCBRPanel
            app.CleantoCompositeBoundRatioCCBRPanel = uipanel(app.UIFigure);
            app.CleantoCompositeBoundRatioCCBRPanel.TitlePosition = 'centertop';
            app.CleantoCompositeBoundRatioCCBRPanel.Title = 'Clean-to-Composite Bound Ratio (CCBR)';
            app.CleantoCompositeBoundRatioCCBRPanel.FontName = 'Calibri';
            app.CleantoCompositeBoundRatioCCBRPanel.FontWeight = 'bold';
            app.CleantoCompositeBoundRatioCCBRPanel.FontSize = 16;
            app.CleantoCompositeBoundRatioCCBRPanel.Position = [591 15 530 228];

            % Create ProcessButton_4
            app.ProcessButton_4 = uibutton(app.CleantoCompositeBoundRatioCCBRPanel, 'push');
            app.ProcessButton_4.ButtonPushedFcn = createCallbackFcn(app, @ProcessButton_4Pushed, true);
            app.ProcessButton_4.FontName = 'Calibri';
            app.ProcessButton_4.Position = [451 19 61 22];
            app.ProcessButton_4.Text = 'Process';

            % Create CCBRtauEditFieldLabel
            app.CCBRtauEditFieldLabel = uilabel(app.CleantoCompositeBoundRatioCCBRPanel);
            app.CCBRtauEditFieldLabel.Interpreter = 'latex';
            app.CCBRtauEditFieldLabel.Position = [37 154 54 22];
            app.CCBRtauEditFieldLabel.Text = 'CCBR $\tau$';

            % Create CCBRtauEditField
            app.CCBRtauEditField = uieditfield(app.CleantoCompositeBoundRatioCCBRPanel, 'numeric');
            app.CCBRtauEditField.ValueDisplayFormat = '%.4g';
            app.CCBRtauEditField.Editable = 'off';
            app.CCBRtauEditField.Position = [106 150 74 30];

            % Create CCBRf_dEditFieldLabel
            app.CCBRf_dEditFieldLabel = uilabel(app.CleantoCompositeBoundRatioCCBRPanel);
            app.CCBRf_dEditFieldLabel.Interpreter = 'latex';
            app.CCBRf_dEditFieldLabel.Position = [37 111 58 22];
            app.CCBRf_dEditFieldLabel.Text = 'CCBR $f_d$';

            % Create CCBRf_dEditField
            app.CCBRf_dEditField = uieditfield(app.CleantoCompositeBoundRatioCCBRPanel, 'numeric');
            app.CCBRf_dEditField.ValueDisplayFormat = '%.4g';
            app.CCBRf_dEditField.Editable = 'off';
            app.CCBRf_dEditField.Position = [106 107 74 30];

            % Create Lamp_CCBR
            app.Lamp_CCBR = uilamp(app.CleantoCompositeBoundRatioCCBRPanel);
            app.Lamp_CCBR.Position = [511 184 10 10];
            app.Lamp_CCBR.Color = [0.851 0.0902 0.0902];

            % Create DualSourceCramrRaoBound2SCRBPanel
            app.DualSourceCramrRaoBound2SCRBPanel = uipanel(app.UIFigure);
            app.DualSourceCramrRaoBound2SCRBPanel.TitlePosition = 'centertop';
            app.DualSourceCramrRaoBound2SCRBPanel.Title = 'Dual Source Cramér-Rao Bound (2S-CRB)';
            app.DualSourceCramrRaoBound2SCRBPanel.FontName = 'Calibri';
            app.DualSourceCramrRaoBound2SCRBPanel.FontWeight = 'bold';
            app.DualSourceCramrRaoBound2SCRBPanel.FontSize = 16;
            app.DualSourceCramrRaoBound2SCRBPanel.Position = [37 15 530 228];

            % Create ProcessButton_3
            app.ProcessButton_3 = uibutton(app.DualSourceCramrRaoBound2SCRBPanel, 'push');
            app.ProcessButton_3.ButtonPushedFcn = createCallbackFcn(app, @ProcessButton_3Pushed, true);
            app.ProcessButton_3.FontName = 'Calibri';
            app.ProcessButton_3.Position = [451 24 61 22];
            app.ProcessButton_3.Text = 'Process';

            % Create sampleLabel_6
            app.sampleLabel_6 = uilabel(app.DualSourceCramrRaoBound2SCRBPanel);
            app.sampleLabel_6.Interpreter = 'latex';
            app.sampleLabel_6.Position = [178 153 55 22];
            app.sampleLabel_6.Text = '[sample]';

            % Create HzLabel_6
            app.HzLabel_6 = uilabel(app.DualSourceCramrRaoBound2SCRBPanel);
            app.HzLabel_6.Interpreter = 'latex';
            app.HzLabel_6.Position = [178 110 55 22];
            app.HzLabel_6.Text = '[Hz]';

            % Create degLabel_6
            app.degLabel_6 = uilabel(app.DualSourceCramrRaoBound2SCRBPanel);
            app.degLabel_6.Interpreter = 'latex';
            app.degLabel_6.Position = [178 24 55 22];
            app.degLabel_6.Text = '[deg]';

            % Create sqrttextnormalCRBtau_0Label
            app.sqrttextnormalCRBtau_0Label = uilabel(app.DualSourceCramrRaoBound2SCRBPanel);
            app.sqrttextnormalCRBtau_0Label.Interpreter = 'latex';
            app.sqrttextnormalCRBtau_0Label.Position = [19 152 61 22];
            app.sqrttextnormalCRBtau_0Label.Text = '$\sqrt{\textnormal{CRB}} \tau_0$';

            % Create CRB2Stau_0EditField
            app.CRB2Stau_0EditField = uieditfield(app.DualSourceCramrRaoBound2SCRBPanel, 'numeric');
            app.CRB2Stau_0EditField.ValueDisplayFormat = '%.4g';
            app.CRB2Stau_0EditField.Editable = 'off';
            app.CRB2Stau_0EditField.Position = [95 148 78 30];

            % Create sqrttextnormalCRBf_d0Label
            app.sqrttextnormalCRBf_d0Label = uilabel(app.DualSourceCramrRaoBound2SCRBPanel);
            app.sqrttextnormalCRBf_d0Label.Interpreter = 'latex';
            app.sqrttextnormalCRBf_d0Label.Position = [19 109 67 22];
            app.sqrttextnormalCRBf_d0Label.Text = '$\sqrt{\textnormal{CRB}} f_{d, 0}$';

            % Create CRB2Sf_d0EditField
            app.CRB2Sf_d0EditField = uieditfield(app.DualSourceCramrRaoBound2SCRBPanel, 'numeric');
            app.CRB2Sf_d0EditField.ValueDisplayFormat = '%.4g';
            app.CRB2Sf_d0EditField.Editable = 'off';
            app.CRB2Sf_d0EditField.Position = [95 105 78 30];

            % Create sqrttextnormalCRBrho_0Label
            app.sqrttextnormalCRBrho_0Label = uilabel(app.DualSourceCramrRaoBound2SCRBPanel);
            app.sqrttextnormalCRBrho_0Label.Interpreter = 'latex';
            app.sqrttextnormalCRBrho_0Label.Position = [19 66 62 22];
            app.sqrttextnormalCRBrho_0Label.Text = '$\sqrt{\textnormal{CRB}} \rho_0$';

            % Create CRB2Srho_0EditField
            app.CRB2Srho_0EditField = uieditfield(app.DualSourceCramrRaoBound2SCRBPanel, 'numeric');
            app.CRB2Srho_0EditField.ValueDisplayFormat = '%.4g';
            app.CRB2Srho_0EditField.Editable = 'off';
            app.CRB2Srho_0EditField.Position = [95 62 78 30];

            % Create sqrttextnormalCRBphi_0Label
            app.sqrttextnormalCRBphi_0Label = uilabel(app.DualSourceCramrRaoBound2SCRBPanel);
            app.sqrttextnormalCRBphi_0Label.Interpreter = 'latex';
            app.sqrttextnormalCRBphi_0Label.Position = [19 23 64 22];
            app.sqrttextnormalCRBphi_0Label.Text = '$\sqrt{\textnormal{CRB}} \phi_0$';

            % Create CRB2Sphi_0EditField
            app.CRB2Sphi_0EditField = uieditfield(app.DualSourceCramrRaoBound2SCRBPanel, 'numeric');
            app.CRB2Sphi_0EditField.ValueDisplayFormat = '%.4g';
            app.CRB2Sphi_0EditField.Editable = 'off';
            app.CRB2Sphi_0EditField.Position = [95 19 78 30];

            % Create sampleLabel_7
            app.sampleLabel_7 = uilabel(app.DualSourceCramrRaoBound2SCRBPanel);
            app.sampleLabel_7.Interpreter = 'latex';
            app.sampleLabel_7.Position = [406 153 55 22];
            app.sampleLabel_7.Text = '[sample]';

            % Create HzLabel_7
            app.HzLabel_7 = uilabel(app.DualSourceCramrRaoBound2SCRBPanel);
            app.HzLabel_7.Interpreter = 'latex';
            app.HzLabel_7.Position = [406 110 55 22];
            app.HzLabel_7.Text = '[Hz]';

            % Create degLabel_7
            app.degLabel_7 = uilabel(app.DualSourceCramrRaoBound2SCRBPanel);
            app.degLabel_7.Interpreter = 'latex';
            app.degLabel_7.Position = [406 24 55 22];
            app.degLabel_7.Text = '[deg]';

            % Create sqrttextnormalCRBtau_1Label
            app.sqrttextnormalCRBtau_1Label = uilabel(app.DualSourceCramrRaoBound2SCRBPanel);
            app.sqrttextnormalCRBtau_1Label.Interpreter = 'latex';
            app.sqrttextnormalCRBtau_1Label.Position = [244 152 61 22];
            app.sqrttextnormalCRBtau_1Label.Text = '$\sqrt{\textnormal{CRB}} \tau_1$';

            % Create CRB2Stau_1EditField
            app.CRB2Stau_1EditField = uieditfield(app.DualSourceCramrRaoBound2SCRBPanel, 'numeric');
            app.CRB2Stau_1EditField.ValueDisplayFormat = '%.4g';
            app.CRB2Stau_1EditField.Editable = 'off';
            app.CRB2Stau_1EditField.Position = [323 148 78 30];

            % Create sqrttextnormalCRBf_d1Label
            app.sqrttextnormalCRBf_d1Label = uilabel(app.DualSourceCramrRaoBound2SCRBPanel);
            app.sqrttextnormalCRBf_d1Label.Interpreter = 'latex';
            app.sqrttextnormalCRBf_d1Label.Position = [244 109 67 22];
            app.sqrttextnormalCRBf_d1Label.Text = '$\sqrt{\textnormal{CRB}} f_{d,1}$';

            % Create CRB2Sf_d1EditField
            app.CRB2Sf_d1EditField = uieditfield(app.DualSourceCramrRaoBound2SCRBPanel, 'numeric');
            app.CRB2Sf_d1EditField.ValueDisplayFormat = '%.4g';
            app.CRB2Sf_d1EditField.Editable = 'off';
            app.CRB2Sf_d1EditField.Position = [323 105 78 30];

            % Create sqrttextnormalCRBrho_1Label
            app.sqrttextnormalCRBrho_1Label = uilabel(app.DualSourceCramrRaoBound2SCRBPanel);
            app.sqrttextnormalCRBrho_1Label.Interpreter = 'latex';
            app.sqrttextnormalCRBrho_1Label.Position = [244 66 62 22];
            app.sqrttextnormalCRBrho_1Label.Text = '$\sqrt{\textnormal{CRB}} \rho_1$';

            % Create CRB2Srho_1EditField
            app.CRB2Srho_1EditField = uieditfield(app.DualSourceCramrRaoBound2SCRBPanel, 'numeric');
            app.CRB2Srho_1EditField.ValueDisplayFormat = '%.4g';
            app.CRB2Srho_1EditField.Editable = 'off';
            app.CRB2Srho_1EditField.Position = [323 62 78 30];

            % Create sqrttextnormalCRBphi_1Label
            app.sqrttextnormalCRBphi_1Label = uilabel(app.DualSourceCramrRaoBound2SCRBPanel);
            app.sqrttextnormalCRBphi_1Label.Interpreter = 'latex';
            app.sqrttextnormalCRBphi_1Label.Position = [244 23 64 22];
            app.sqrttextnormalCRBphi_1Label.Text = '$\sqrt{\textnormal{CRB}} \phi_1$';

            % Create CRB2Sphi_1EditField
            app.CRB2Sphi_1EditField = uieditfield(app.DualSourceCramrRaoBound2SCRBPanel, 'numeric');
            app.CRB2Sphi_1EditField.ValueDisplayFormat = '%.4g';
            app.CRB2Sphi_1EditField.Editable = 'off';
            app.CRB2Sphi_1EditField.Position = [323 19 78 30];

            % Create Lamp_2SCRB
            app.Lamp_2SCRB = uilamp(app.DualSourceCramrRaoBound2SCRBPanel);
            app.Lamp_2SCRB.Position = [508 184 10 10];
            app.Lamp_2SCRB.Color = [0.851 0.0902 0.0902];

            % Create SignalparametersPanel
            app.SignalparametersPanel = uipanel(app.UIFigure);
            app.SignalparametersPanel.TitlePosition = 'centertop';
            app.SignalparametersPanel.Title = 'Signal parameters';
            app.SignalparametersPanel.FontName = 'Calibri';
            app.SignalparametersPanel.FontWeight = 'bold';
            app.SignalparametersPanel.FontSize = 16;
            app.SignalparametersPanel.Position = [37 501 1083 268];

            % Create basebandsignalsamplesEditFieldLabel
            app.basebandsignalsamplesEditFieldLabel = uilabel(app.SignalparametersPanel);
            app.basebandsignalsamplesEditFieldLabel.Interpreter = 'latex';
            app.basebandsignalsamplesEditFieldLabel.Position = [723 145 159 22];
            app.basebandsignalsamplesEditFieldLabel.Text = 'baseband signal samples';

            % Create basebandsignalsamplesEditField
            app.basebandsignalsamplesEditField = uieditfield(app.SignalparametersPanel, 'text');
            app.basebandsignalsamplesEditField.Editable = 'off';
            app.basebandsignalsamplesEditField.HorizontalAlignment = 'center';
            app.basebandsignalsamplesEditField.Position = [894 145 100 22];
            app.basebandsignalsamplesEditField.Value = 'none';

            % Create LoadButton
            app.LoadButton = uibutton(app.SignalparametersPanel, 'push');
            app.LoadButton.ButtonPushedFcn = createCallbackFcn(app, @LoadButtonPushed, true);
            app.LoadButton.FontName = 'Calibri';
            app.LoadButton.Position = [1005 145 61 22];
            app.LoadButton.Text = 'Load';

            % Create Signal1reflectedNLOSPanel
            app.Signal1reflectedNLOSPanel = uipanel(app.SignalparametersPanel);
            app.Signal1reflectedNLOSPanel.BorderType = 'none';
            app.Signal1reflectedNLOSPanel.TitlePosition = 'centertop';
            app.Signal1reflectedNLOSPanel.Title = 'Signal 1 (reflected, NLOS)';
            app.Signal1reflectedNLOSPanel.FontName = 'Calibri';
            app.Signal1reflectedNLOSPanel.FontWeight = 'bold';
            app.Signal1reflectedNLOSPanel.FontSize = 14;
            app.Signal1reflectedNLOSPanel.Position = [383 18 319 215];

            % Create timedelaytau_1EditField_2Label
            app.timedelaytau_1EditField_2Label = uilabel(app.Signal1reflectedNLOSPanel);
            app.timedelaytau_1EditField_2Label.Interpreter = 'latex';
            app.timedelaytau_1EditField_2Label.HorizontalAlignment = 'right';
            app.timedelaytau_1EditField_2Label.Position = [78 151 82 22];
            app.timedelaytau_1EditField_2Label.Text = 'time-delay $\tau_1$';

            % Create timedelaytau_1EditField
            app.timedelaytau_1EditField = uieditfield(app.Signal1reflectedNLOSPanel, 'numeric');
            app.timedelaytau_1EditField.ValueChangedFcn = createCallbackFcn(app, @parameterSignal_1EditFieldValueChanged, true);
            app.timedelaytau_1EditField.Position = [175 147 58 30];
            app.timedelaytau_1EditField.Value = 2;

            % Create Dopplerfrequencyf_d1EditField_2Label
            app.Dopplerfrequencyf_d1EditField_2Label = uilabel(app.Signal1reflectedNLOSPanel);
            app.Dopplerfrequencyf_d1EditField_2Label.Interpreter = 'latex';
            app.Dopplerfrequencyf_d1EditField_2Label.HorizontalAlignment = 'right';
            app.Dopplerfrequencyf_d1EditField_2Label.Position = [27 108 133 22];
            app.Dopplerfrequencyf_d1EditField_2Label.Text = 'Doppler frequency $f_{d,1}$';

            % Create Dopplerfrequencyf_d1EditField
            app.Dopplerfrequencyf_d1EditField = uieditfield(app.Signal1reflectedNLOSPanel, 'numeric');
            app.Dopplerfrequencyf_d1EditField.ValueChangedFcn = createCallbackFcn(app, @parameterSignal_1EditFieldValueChanged, true);
            app.Dopplerfrequencyf_d1EditField.Position = [175 104 58 30];
            app.Dopplerfrequencyf_d1EditField.Value = 100;

            % Create amplituderho_1EditField_2Label
            app.amplituderho_1EditField_2Label = uilabel(app.Signal1reflectedNLOSPanel);
            app.amplituderho_1EditField_2Label.Interpreter = 'latex';
            app.amplituderho_1EditField_2Label.HorizontalAlignment = 'right';
            app.amplituderho_1EditField_2Label.Position = [81 65 79 22];
            app.amplituderho_1EditField_2Label.Text = 'amplitude $\rho_1$';

            % Create amplituderho_1EditField
            app.amplituderho_1EditField = uieditfield(app.Signal1reflectedNLOSPanel, 'numeric');
            app.amplituderho_1EditField.ValueChangedFcn = createCallbackFcn(app, @parameterSignal_1EditFieldValueChanged, true);
            app.amplituderho_1EditField.Position = [175 61 58 30];
            app.amplituderho_1EditField.Value = 0.5;

            % Create phasephi_1EditField_2Label
            app.phasephi_1EditField_2Label = uilabel(app.Signal1reflectedNLOSPanel);
            app.phasephi_1EditField_2Label.Interpreter = 'latex';
            app.phasephi_1EditField_2Label.HorizontalAlignment = 'right';
            app.phasephi_1EditField_2Label.Position = [103 22 57 22];
            app.phasephi_1EditField_2Label.Text = 'phase $\phi_1$';

            % Create phasephi_1EditField
            app.phasephi_1EditField = uieditfield(app.Signal1reflectedNLOSPanel, 'numeric');
            app.phasephi_1EditField.ValueChangedFcn = createCallbackFcn(app, @parameterSignal_1EditFieldValueChanged, true);
            app.phasephi_1EditField.Position = [175 18 58 30];
            app.phasephi_1EditField.Value = 15;

            % Create sampleLabel
            app.sampleLabel = uilabel(app.Signal1reflectedNLOSPanel);
            app.sampleLabel.Interpreter = 'latex';
            app.sampleLabel.Position = [238 152 55 22];
            app.sampleLabel.Text = '[sample]';

            % Create HzLabel
            app.HzLabel = uilabel(app.Signal1reflectedNLOSPanel);
            app.HzLabel.Interpreter = 'latex';
            app.HzLabel.Position = [238 109 55 22];
            app.HzLabel.Text = '[Hz]';

            % Create degLabel
            app.degLabel = uilabel(app.Signal1reflectedNLOSPanel);
            app.degLabel.Interpreter = 'latex';
            app.degLabel.Position = [238 23 55 22];
            app.degLabel.Text = '[deg]';

            % Create Signal0directLOSPanel
            app.Signal0directLOSPanel = uipanel(app.SignalparametersPanel);
            app.Signal0directLOSPanel.BorderType = 'none';
            app.Signal0directLOSPanel.TitlePosition = 'centertop';
            app.Signal0directLOSPanel.Title = 'Signal 0 (direct, LOS)';
            app.Signal0directLOSPanel.FontName = 'Calibri';
            app.Signal0directLOSPanel.FontWeight = 'bold';
            app.Signal0directLOSPanel.FontSize = 14;
            app.Signal0directLOSPanel.Position = [25 18 319 215];

            % Create timedelaytau_0EditField_2Label
            app.timedelaytau_0EditField_2Label = uilabel(app.Signal0directLOSPanel);
            app.timedelaytau_0EditField_2Label.Interpreter = 'latex';
            app.timedelaytau_0EditField_2Label.HorizontalAlignment = 'right';
            app.timedelaytau_0EditField_2Label.Position = [79 148 82 22];
            app.timedelaytau_0EditField_2Label.Text = 'time-delay $\tau_0$';

            % Create timedelaytau_0EditField
            app.timedelaytau_0EditField = uieditfield(app.Signal0directLOSPanel, 'numeric');
            app.timedelaytau_0EditField.ValueChangedFcn = createCallbackFcn(app, @parameterSignal_0EditFieldValueChanged, true);
            app.timedelaytau_0EditField.Position = [176 144 58 30];

            % Create Dopplerfrequencyf_d0EditField_3Label
            app.Dopplerfrequencyf_d0EditField_3Label = uilabel(app.Signal0directLOSPanel);
            app.Dopplerfrequencyf_d0EditField_3Label.Interpreter = 'latex';
            app.Dopplerfrequencyf_d0EditField_3Label.HorizontalAlignment = 'right';
            app.Dopplerfrequencyf_d0EditField_3Label.Position = [28 105 133 22];
            app.Dopplerfrequencyf_d0EditField_3Label.Text = 'Doppler frequency $f_{d,0}$';

            % Create Dopplerfrequencyf_d0EditField
            app.Dopplerfrequencyf_d0EditField = uieditfield(app.Signal0directLOSPanel, 'numeric');
            app.Dopplerfrequencyf_d0EditField.Position = [176 101 58 30];

            % Create amplituderho_0EditField_2Label
            app.amplituderho_0EditField_2Label = uilabel(app.Signal0directLOSPanel);
            app.amplituderho_0EditField_2Label.Interpreter = 'latex';
            app.amplituderho_0EditField_2Label.HorizontalAlignment = 'right';
            app.amplituderho_0EditField_2Label.Position = [82 62 79 22];
            app.amplituderho_0EditField_2Label.Text = 'amplitude $\rho_0$';

            % Create amplituderho_0EditField
            app.amplituderho_0EditField = uieditfield(app.Signal0directLOSPanel, 'numeric');
            app.amplituderho_0EditField.Position = [176 58 58 30];
            app.amplituderho_0EditField.Value = 1;

            % Create phasephi_0EditField_2Label
            app.phasephi_0EditField_2Label = uilabel(app.Signal0directLOSPanel);
            app.phasephi_0EditField_2Label.Interpreter = 'latex';
            app.phasephi_0EditField_2Label.HorizontalAlignment = 'right';
            app.phasephi_0EditField_2Label.Position = [104 19 57 22];
            app.phasephi_0EditField_2Label.Text = 'phase $\phi_0$';

            % Create phasephi_0EditField
            app.phasephi_0EditField = uieditfield(app.Signal0directLOSPanel, 'numeric');
            app.phasephi_0EditField.Position = [176 15 58 30];

            % Create sampleLabel_2
            app.sampleLabel_2 = uilabel(app.Signal0directLOSPanel);
            app.sampleLabel_2.Interpreter = 'latex';
            app.sampleLabel_2.Position = [239 149 55 22];
            app.sampleLabel_2.Text = '[sample]';

            % Create HzLabel_2
            app.HzLabel_2 = uilabel(app.Signal0directLOSPanel);
            app.HzLabel_2.Interpreter = 'latex';
            app.HzLabel_2.Position = [239 106 55 22];
            app.HzLabel_2.Text = '[Hz]';

            % Create degLabel_2
            app.degLabel_2 = uilabel(app.Signal0directLOSPanel);
            app.degLabel_2.Interpreter = 'latex';
            app.degLabel_2.Position = [239 20 55 22];
            app.degLabel_2.Text = '[deg]';

            % Create HistoryTextArea
            app.HistoryTextArea = uitextarea(app.SignalparametersPanel);
            app.HistoryTextArea.Editable = 'off';
            app.HistoryTextArea.Position = [723 14 343 63];
            app.HistoryTextArea.Value = {'Welcome to the CRB overviewer!'};

            % Create SNR_textnormaloutSliderLabel
            app.SNR_textnormaloutSliderLabel = uilabel(app.SignalparametersPanel);
            app.SNR_textnormaloutSliderLabel.Interpreter = 'latex';
            app.SNR_textnormaloutSliderLabel.HorizontalAlignment = 'center';
            app.SNR_textnormaloutSliderLabel.Position = [723 95 48 22];
            app.SNR_textnormaloutSliderLabel.Text = 'SNR$_{\textnormal{out}}$';

            % Create SNR_textnormaloutSlider
            app.SNR_textnormaloutSlider = uislider(app.SignalparametersPanel);
            app.SNR_textnormaloutSlider.Limits = [0 60];
            app.SNR_textnormaloutSlider.MajorTicks = [0 6 12 18 24 30 36 42 48 54 60];
            app.SNR_textnormaloutSlider.Position = [783 113 272 3];
            app.SNR_textnormaloutSlider.Value = 24;

            % Create samplingfrenquencyfactorEditFieldLabel
            app.samplingfrenquencyfactorEditFieldLabel = uilabel(app.SignalparametersPanel);
            app.samplingfrenquencyfactorEditFieldLabel.Interpreter = 'latex';
            app.samplingfrenquencyfactorEditFieldLabel.Position = [723 179 160 22];
            app.samplingfrenquencyfactorEditFieldLabel.Text = 'sampling frenquency factor';

            % Create samplingfrenquencyfactorEditField
            app.samplingfrenquencyfactorEditField = uieditfield(app.SignalparametersPanel, 'numeric');
            app.samplingfrenquencyfactorEditField.ValueChangedFcn = createCallbackFcn(app, @samplingfrenquencyfactorEditFieldValueChanged, true);
            app.samplingfrenquencyfactorEditField.HorizontalAlignment = 'left';
            app.samplingfrenquencyfactorEditField.Position = [894 179 100 22];
            app.samplingfrenquencyfactorEditField.Value = 1;

            % Create Lamp_signal
            app.Lamp_signal = uilamp(app.SignalparametersPanel);
            app.Lamp_signal.Position = [1031 186 10 10];
            app.Lamp_signal.Color = [0.851 0.0902 0.0902];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = main

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