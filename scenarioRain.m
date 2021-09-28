%% General sim preparations for Rain Conditions
% clear up
clear all variables; clc;

dofig  = false;
dosavedata = true;
dosavefig  = true;

% timeline
Sim.dt       = 0.1;             % sec
Sim.Tsim     = 15*60;           % sec
Sim.timeline = 0:Sim.dt:Sim.Tsim-Sim.dt;    % sec
Sim.numsteps = numel(Sim.timeline);
% update progress every ... steps
Sim.dtprogress = 450;
% road stretch
Sim.xstart =  0;
Sim.xmax   = 6e3;  % exit at 3 km
% max nr of vehicles generated:
Sim.numvehicles = 350;
% Safe zone: we cannot apply reaction time logic close to origin.
Sim.safezone = 5*40; % 5 (sec) * 40 (m/s) (ultra long reaction time * very high speed)

%% Scenario preparations

% base OD profile (sine curve)
q = zeros(size(Sim.timeline));
sinebase  = 1+sin(1.5*pi/Sim.Tsim * Sim.timeline);
% qbase = 900 * sinebase;
qbase = [900  * ones(1,0.1*Sim.numsteps), ...
    2200 * ones(1,0.3*Sim.numsteps), ...
    900  * ones(1,0.6*Sim.numsteps)];

%qbase = 1400  * ones(1,Sim.numsteps);
%q = poissrnd(qbase);
Sim.q = qbase;
%% here's the list of scenario's
scenarios = { ...
%      'scen0_basecase_dry',...
%      'scen0_basecase_wetTarmac',...
% %      'scen0_basecase_wetTarmac_mu=0.45',...
%       'scen1_I=0.5',...
        'scen1_I=1.0',...
%       'scen1_I=1.5',...
%       'scen1_I=2.0',...
%       'scen1_I=2.5',...
%       'scen1_I=3.0',...
%       'scen1_I=3.5',...
%       'scen1_I=4.0',...
%       'scen1_I=4.5',...
%       'scen1_I=5.0'...
};

% run over all of them
numscen = numel(scenarios);

% 
for itau = 1
    Sim.base_tau    = itau/10;
    tau_str = sprintf('tau=%.2d',itau);
    
    fprintf('\n-------------------------\n');
    fprintf('BASE REACTION TIME = %.1f  \n',Sim.base_tau);
    fprintf('-------------------------\n');
    
    % over all the scenarios
    for iscenario = 1 : numscen
        
        % BASE Settings for the scenarios
        
       % A bottleneck at 4 and 4.5
        Sim.xDisturbance  = [4e3,4.5e3];    % disturbance location
        Sim.tDisturbance  = [1, 15] * 60;   % disturbance
        Sim.disturbanceStrength = 0.36;      % Severity (qcap = 2100 veh/h, u = 33.33 m/s)
    
         % Base scenario:
        Sim.doHF        = true;      % do we compute HF?
        Sim.doDistraction = false;   % do we simulate a distraction?
        Sim.doDisturbance = true;    % do we force an additional disturbances, lets assume a flow conserving bottleneck?
        Sim.doHeterogeneity = false; % do we vary with critical task saturation?
        %Sim.base_tau    = itau;     % do we give drivers a base reaction time?
        Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
        Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
        Sim.doHF_tau    = false;      % do we have HF affect reaction time (attention lag)
        Sim.doHF_v0     = false;     % do we have HF affect v0
        Sim.doHF_T      = false;     % do we have HF affect T
        Sim.doHF_TD     = false;     % do we have HF affect Task Difficulty
        
        Sim.doHF_perc_s_over = 0.5;  % fraction of drivers that overestimates
        % distance gaps/speeds if set to 0, 100%
        % underestimates, if set to 1 100%
        % overestimates, and any number in between
        % yields a mix (default: 50-50)
        % selected vehicles for HF plots
        Sim.iHF = [60,130];
        
        % select scenario
        Sim.scenario = scenarios{iscenario};
        switch Sim.scenario
		%% The following set of scenarios, regard a face calibration of the model when...
		...only the appearance of wet surface is on!!!
			case 'scen0_basecase_dry' % like martini ho ho ho!
                Sim.mu = 0.6;
                Sim.ri = 0.0;
                Sim.R  = 1.0;
                Sim.doHeterogeneity = true; % do we vary with critical task saturation?
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;      % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;     % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
%             case 'scen0_basecase_wetTarmac'
%                 Sim.mu = 0.4;
%                 Sim.ri = 0.0;
%                 Sim.R  = 1.0;
%                 Sim.doHeterogeneity = true; % do we vary with critical task saturation?
%                 Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
%                 Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
%                 Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
%                 Sim.doHF_v0     = true;     % do we have HF affect v0
%                 Sim.doHF_T      = true;     % do we have HF affect T
            case 'scen0_basecase_wetTarmac_mu=0.45'
                Sim.mu = 0.40;
                Sim.ri = 0.0;
                Sim.R  = 1.0;
                Sim.doHeterogeneity = true; % do we vary with critical task saturation?
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;      % do we have HF affect reaction time (attention lag)
%                 Sim.doHF_v0     = true;     % do we have HF affect v0
%                 Sim.doHF_T      = true;     % do we have HF affect T
           case 'scen1_I=0.5'
                Sim.mu = 0.40;
                Sim.ri = 0.5;
%                 Sim.doHeterogeneity = true; % do we vary with critical task saturation?
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;      % do we have HF affect reaction time (attention lag)
%               Sim.doHF_v0     = true;     % do we have HF affect v0
%               Sim.doHF_T      = true;     % do we have HF affect T
                
           case 'scen1_I=1.0'
                Sim.mu = 0.40;
                Sim.ri = 1.0;
                Sim.doHeterogeneity = true; % do we vary with critical task saturation?
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;      % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;     % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
                 inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T,Sim.doHF_TD);
            case 'scen1_I=1.5'
                Sim.mu = 0.40;
                Sim.ri = 1.5;
%                 Sim.doHeterogeneity = true; % do we vary with critical task saturation?
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;      % do we have HF affect reaction time (attention lag)
%                 Sim.doHF_v0     = true;     % do we have HF affect v0
%                 Sim.doHF_T      = true;     % do we have HF affect T   
           case 'scen1_I=2.0'
                Sim.mu = 0.40;
                Sim.ri = 2.0;
                Sim.doHeterogeneity = true; % do we vary with critical task saturation?
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;      % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;     % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T 
           case 'scen1_I=2.5'
                Sim.mu = 0.40;
                Sim.ri = 2.5;
                Sim.doHeterogeneity = true; % do we vary with critical task saturation?
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;      % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;     % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T 
                case 'scen1_I=3.0'
                Sim.mu = 0.40;
                Sim.ri = 3.0;
                Sim.doHeterogeneity = true; % do we vary with critical task saturation?
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;      % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;     % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T 
                case 'scen1_I=3.5'
                Sim.mu = 0.40;
                Sim.ri = 3.5;
                Sim.doHeterogeneity = true; % do we vary with critical task saturation?
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;      % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;     % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T 
                case 'scen1_I=4.0'
                Sim.mu = 0.40;
                Sim.ri = 4.0;
                Sim.doHeterogeneity = true; % do we vary with critical task saturation?
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;      % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;     % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
                case 'scen1_I=4.5'
                Sim.mu = 0.40;
                Sim.ri = 4.5;
                Sim.doHeterogeneity = true; % do we vary with critical task saturation?
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;      % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;     % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T 
                case 'scen1_I=5.0'
                Sim.mu = 0.40;
                Sim.ri = 5.0;
%                 Sim.doHeterogeneity = true; % do we vary with critical task saturation?
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;      % do we have HF affect reaction time (attention lag)
%                 Sim.doHF_v0     = true;     % do we have HF affect v0
%                 Sim.doHF_T      = true;     % do we have HF affect T 
                
        end
   
     % close current figs
        close all;
        
        % run simulations
        Sim = simulation_fcn(Sim);
        
        % recieve information from the detectors
        Sim = virtDet(Sim);
        
        % figures
        if dofig
            % visualisations
            figs = visualisation_fcn(Sim);
        end
        
        % save data?
        if dosavedata
            disp('saving everything ...');
            
            if ~isdir(fullfile(pwd,'RainResults'))
                mkdir(pwd,'RainResults')
            end;
            fname = sprintf('RainResults/%s_%s.mat',Sim.scenario, tau_str);
            save(fname,'Sim');
        end;
        
        % save figures?
        if dosavefig
            % we need the figs ...
            if ~dofig
                figs = visualisation_fcn(Sim);
            end;
            % export graphs to png
            if ~isdir(fullfile(pwd,'RainImages'))
                mkdir(pwd,'RainImages');
            end;
            saveas(figs.Traj  ,sprintf('RainImages/%s_%s_traj.png',Sim.scenario,tau_str,'-r300'));
            saveas(figs.HF    ,sprintf('RainImages/%s_%s_hf.png',Sim.scenario,tau_str,'-r300'));
            saveas(figs.FDTD  ,sprintf('RainImages/%s_%s_fdtd.png',Sim.scenario,tau_str,'-r300'));
            saveas(figs.FD    ,sprintf('RainImages/%s_%s_fd.png',Sim.scenario,tau_str,'-r300'));
            saveas(figs.FDsel ,sprintf('RainImages/%s_%s_fdsel.png',Sim.scenario,tau_str,'-r300'));
            saveas(figs.Demand,sprintf('RainImages/%s_%s_demand.png',Sim.scenario,tau_str,'-r300'));
			saveas(figs.HEADACC,sprintf('RainImages/%s_%s_accheadspeed.png',Sim.scenario,tau_str,'-r300'));
% 			saveas(figs.PPD,sprintf('RainImages/%s_%s_virtualDetectors.png',Sim.scenario,tau_str));
        end
        
        disp('DONE');
    end
end

Sim.numscen = numscen;
return;
%% Below follow cells with code to generate the two-column result graphs 
%  and the results table used in the manuscript

%% Compute TTS and num inc table
tbl_TTS = zeros(5, numscen);
tbl_INC = zeros(5, numscen);
tbl_TDC = zeros(5, numscen); % density
taus = [1,2,3,4];
for itau = 1:numel(taus)
    aTaustr = sprintf('tau=%.2d',taus(itau));
    for iscenario = 1:numscen
        aScenario = scenarios{iscenario};
        fname = sprintf('RainResults/%s_%s.mat',aScenario, aTaustr);
        load(fname);
        tbl_TTS(itau, iscenario) = Sim.TTS;
        tbl_INC(itau, iscenario) = Sim.numCollisions;
        tbl_TDC(itau, iscenario) = Sim.TDC;
    end
end;

fname = 'tbls2.xls'
xlswrite(fname,scenarios,'A1');
xlswrite(fname,[tbl_TTS; tbl_INC; tbl_TDC],'A2');

    
