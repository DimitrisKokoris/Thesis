function [Sim] = simulation_fcn(Sim, scenario)
%[Sim] = simulation(Sim):
% Simulate with IDMplus model for Tr-B paper%% prepare all variables
% For some reason , some matrices are still visible for the next siulation, hence:
tic;
if nargin<2
    scenario = '';
end

% first argument is a structure or a filename
if ~isstruct(Sim)
    Sim = load(Sim);
end

% Copy all vars that are in Sim structure to root shortcuts (e.g. Sim.dt
% becomes dt)
flds = fieldnames(Sim);
for i = 1:numel(flds)
    evalstr = sprintf('%s = Sim.%s;',flds{i}, flds{i});
    eval(evalstr);
end

% add fields that are generated in this function
flds{end+1} = 'VEH';
flds{end+1} = 'TTS';
flds{end+1} = 'H';
flds{end+1} = 'S';
flds{end+1} = 'V';
flds{end+1} = 't';
flds{end+1} = 'numCollisions';
flds{end+1} = 'Collisions';
flds{end+1} = 'TTS_bar';
flds{end+1} = 'TTC_crit';
flds{end+1} = 'RBR_crit';
flds{end+1} = 'mean_Speed_Edie';
flds{end+1} = 'TDC';

%% prepare all variables

% compute vehicle generation profile in units of dt
h_demand  = dt*round((3600./q)./dt);

% create veh array
VEH(numvehicles) = createDefVehicle(Sim);
t = 0; k=1;
for n = 1:numvehicles
    VEH(n) = createDefVehicle(Sim);
    % make sure we tage the vehicle so that we find it also when we select
    % subsections of this list
    VEH(n).ID = n;
    % we now compute the desired generation time (in units of dt). If n==1
    % we generate a vehicle at t=0. In all other cases, we generate
    % according to the demand pattern
    if n > 1
        % find time instant in time line
        k = find(timeline >= t);
        if isempty(k)
            k = numsteps;
        end;
    end;
    % set values
    VEH(n).kgendes = k(1);
    VEH(n).tgendes = t;
    % advance to that time instant
    t = t + h_demand(k(1));
    
    % we now preallocate the arrays for trajectory data (we truncate
    % afterwards)
%     VEH(n).traj = zeros(numsteps,6);
%     if doHF
%         VEH(n).trajHF = zeros(numsteps,7);
%     end;
    
    % check if we're done with creating vehicles
    if k == numsteps
        break;
    end
end


%% Actual Simulation
fprintf('SCENARIO: %s...\n',scenario);
numCollisions = 0;
Collisions = [];

iVeh   = 1;     % index of next vehicle
doGenerate = false; % if set true: we can generate it

% if set, all vehicles after this ID experience a shorter track;
IsIncident = false;

% Here generate the estimation errors!
[matrix, actProb] = estimation(numvehicles, doHF_prob_overestimation);
fprintf('Probability of underestimation is: %0.2f...\n',1 - actProb);
fprintf('Probability of overestimation is: %0.2f...\n',actProb);

% run the simulation
for k= 1:numsteps
    
    % current time instant
    t = timeline(k);
    
    % indices for vehicles on track
    active   = [VEH.active];
    finished = [VEH.finished];
    ontrack  = active & ~finished;
    
    %% ---- Vehicle Generation ---------------------------------------------
    
    % if we're out of vehicles somebody (me) should do something ;-)
    if iVeh>numvehicles
        doGenerate = false;
        % if this is the first vehicle OR if the previous vehicle is no longer
        % on the track we simply generate it
    elseif iVeh == 1 || ~ontrack(iVeh-1)
        doGenerate = true;
    else
        % (1) is it time for this vehicle to enter?
        shouldenter = k>=VEH(iVeh).kgendes;
        % if it should enter we check a few things
        if shouldenter
            % the NET gap to the leading vehicle equals
            VEH(iVeh).s = (VEH(iVeh-1).x - VEH(iVeh).par.s0);
            % the NET headway to the leading vehicle equals
            VEH(iVeh).h = VEH(iVeh).s/VEH(iVeh).v;
            % desired headway of the present vehicle
            h_des = VEH(iVeh).par.T;
            % if there's no space: no vehicle
            collision = 1.1*VEH(iVeh).par.s0 >= VEH(iVeh-1).x;
            % we will generate this vehicle ONLY if this fits
            % with the vehicles own preferences AND if no collision occurs
            doGenerate = ~collision && VEH(iVeh).h >= h_des;
            % if the vehicle cannot be generated we add queueing time
            VEH(iVeh).tqueue = VEH(iVeh).tqueue + ~doGenerate * dt;
        else
            doGenerate = false;
        end
    end
    
    % if doGenerate is set true, we activate the next vehicle and update
    % the other variables
    if doGenerate
        % activate vehicle
        VEH(iVeh).active    = true;
        % make sure all fields that require the right starting time are
        % filled ...
        VEH(iVeh).kgen      = k;
        VEH(iVeh).tgen      = t;
        VEH(iVeh).t         = t;
        VEH(iVeh).x         = xstart;
        % we adjust the speed to avoid crashes / smooth entrance
        if iVeh>1
            VEH(iVeh).v     = VEH(iVeh-1).v;
            % the GROSS gap to the leading vehicle equals
            VEH(iVeh).s     = VEH(iVeh-1).x;
        else
            VEH(iVeh).v     = VEH(iVeh).par.v0;
            % no vehicle on the track yet ...
            VEH(iVeh).s     = xmax;
        end
        % the GROSS headway to the leading vehicle equals
        VEH(iVeh).h = VEH(iVeh).s/VEH(iVeh).v;
        % acceleration is zero. 
        
        % what is the mix of over- and under perception errors?
        % doHF_perc_s_over is a number between 0 and 1, so if it is 0,
        % every one underestimates (percbias=-1), if it is 1, every driver
        % overestimates (percbias=1), and otherwise a mix.
        % VEH(iVeh).percbias = sign(doHF_perc_s_over - rand);
          VEH(iVeh).percbias =  matrix(iVeh);
        
        % Heterogeneity?
        if doHeterogeneity
            % desired speed
            %VEH(iVeh).par.v0      = min(40,  max(30,  35 + randn));
            % desired gap
            %VEH(iVeh).par.T       = min(1.6, max(0.8, 1.2 + 0.1*randn));
            % critical saturation (proxy for skill level)
            % VEH(iVeh).SAFn.tscrit = min(1.0, max(0.5, 0.75 + 0.1*randn));
            % we vary with task capacity (proxy for skill level)
            %VEH(iVeh).TC = min(1.2, max(0.8, 1 + 0.1*randn));
            VEH(iVeh).TC = min(1.1, max(0.9, 1 + 0.05*randn));
        end     
        % reaction time!!!
        VEH(iVeh).par.tau = base_tau;
        
        % Finally, we now construct trajectory:
        % movement
        VEH(iVeh).traj = [VEH(iVeh).t,VEH(iVeh).x, ...
            VEH(iVeh).v, VEH(iVeh).a, VEH(iVeh).s,VEH(iVeh).h];
        % HF
        if doHF
            VEH(iVeh).trajHF = [ ...
                    VEH(i).t,VEH(i).x,VEH(i).TDcf,VEH(i).TDdi, VEH(i).SA, ...
                    VEH(i).ctl.tauextra + VEH(i).par.tau, inf];
        end
        
%         VEH(iVeh).traj(k,:) = [VEH(iVeh).t,VEH(iVeh).x, ...
%             VEH(iVeh).v, VEH(iVeh).a, VEH(iVeh).s,VEH(iVeh).h];
%         % HF
%         if doHF
%             VEH(iVeh).trajHF(k,:) = [ ...
%                     VEH(i).t,VEH(i).x,VEH(i).TDcf,VEH(i).TDdi, VEH(i).SA, ...
%                     VEH(i).ctl.tauextra + VEH(i).par.tau, inf];
%         end;
        
        % update counter
        iVeh = iVeh + 1;
        % update flag
        doGenerate = false;
    end
    
    % redo indices for vehicles on track
    active   = [VEH.active];
    ontrack  = active & ~finished;
    
    % ------ Behaviour & MOVEMENT OF VEHICLES -------
    % we bother with the selection on track only
    iVEH = find(ontrack);
    % furthest vehicle on the track ...
    iMin = min(iVEH);
    
    % we run in increasing vehicle number order to determine
    % (a) HF dynamics
    % (b) accelerations
    for i = iVEH
        
        %% ----- PERCEPTION ------------------------------------------------
        
        % first vehicle goes free
        if i == iMin
            ds = inf; dv = inf;
            dv_gt = inf; ds_gt = inf;
            
            % all others have a leader and need to perceive the stimuli (gap
            % and speed difference)
        else
            % ground-truth NET gaps and speed differences
            ds_gt = VEH(i-1).x - (VEH(i).x + VEH(i).par.s0);
            dv_gt = VEH(i).v   - VEH(i-1).v;
            
            % We cannot apply the reaction time logic close to vehicle
            % generation, so we only do this after passing the safezone
            if VEH(i).x <= safezone
                ds = ds_gt;
                dv = dv_gt;
            else
                % get reaction time (which might be already longer than normal)
                tau = VEH(i).par.tau + VEH(i).ctl.tauextra;
                % translate reaction time to simsteps for VEH(i)
                ktau = max(0, round(tau / dt));
                % figure out what VEH(i) and VEH(i-1) were doing at the
                % time: get speed and gap differences of that time slice.
%                 ds = VEH(i-1).traj(k-ktau, 2) - (VEH(i).traj(k-ktau, 2) + VEH(i).par.s0);
%                 dv = VEH(i).traj(k-ktau, 3)   - VEH(i-1).traj(k-ktau, 3);
                ds = VEH(i-1).traj(end-ktau, 2) - (VEH(i).traj(end-ktau, 2) + VEH(i).par.s0);
                dv = VEH(i).traj(end-ktau, 3)   - VEH(i-1).traj(end-ktau, 3);
            end;
        end;
        
        % current perceived headway
        h = ds / VEH(i).v;
        
        % actual gap and headway
        VEH(i).s = ds_gt;
        VEH(i).h = ds_gt / VEH(i).v;
        
        %% ----------- HF Dynamics ----------------------------------------
        
        % distraction?
        if doDistraction && t>=tDistraction(1) && t<tDistraction(2)
            xdist = VEH(i).x - xDistraction;
        else
            xdist = inf;
        end;
        
        % Do we apply human factors? these are all based on perceived
        % stimuli!
        if doHF
            % we first compute the HF state (task demands, task saturation,
            % SA levels, etc)
            VEH(i) = computeHF(VEH(i), h, xdist);
            
            % SA deterioration factor (between 0-1)
            sa_fac = min(VEH(i).SAFn.samax - VEH(i).SA, 1);
            
            % response adaptation factor (between 0-1)
            rs_adapt = min( max(VEH(i).TS - VEH(i).SAFn.tscrit, 0), 1);
            
            % we then apply these to ...
            
            % (1) perception quality as a function of deterioration in SA.
            % Drivers either OVER- or UNDERESTIMATE gaps and
            % speeds, we assume they get worse at it under pressure
            perc_fac = 1;
            if doHF_perc_s
                ds = ds + VEH(i).percbias * sa_fac * perc_fac * ds;
            end;
            if doHF_perc_v
                dv = dv + VEH(i).percbias * sa_fac * perc_fac * dv;
            end;
            
            % (2) reaction time increase...? This extra attention
            % lag is assumed proportional to deterioration in SA
            if doHF_tau
                VEH(i).ctl.tauextra = round((sa_fac + 1) * VEH(i).SAFn.taumax,1);
            end;
            
            % (3) response adaptation: v0? Assumption: proportional to task
            % saturation with a maximum of 90% (which is full in the
            % brakes)
            v0_fac = .9;
            if doHF_v0
                VEH(i).ctl.v0fac = 1 - v0_fac*rs_adapt;
            end;
            
            % (4) response adaptation: T?
            % we assume that drivers compensate proportionally to task
            % saturation with a maximum of 100%.
            T_fac = 1;
            if doHF_T
                VEH(i).ctl.Tfac = 1 + T_fac*rs_adapt;
            end;
            
        end;
        
        %% ----- Incidents, bottlenecks, etc ------------------------------
        % time to collision
        if dv_gt > 0
            ttc = ds_gt/dv_gt;
        else
            ttc = NaN;
        end;

        % check for collisions
        if ds_gt<=0 && ~VEH(i).incident
            VEH(i).incident = true;
            VEH(i).inc.x = VEH(i).x;
            VEH(i).inc.t = t;
            % report
            numCollisions = numCollisions + 1;
            Collisions(numCollisions).rear = i-1;
            Collisions(numCollisions).front = i;
            Collisions(numCollisions).x = VEH(i).x;
            Collisions(numCollisions).k = k;
            Collisions(numCollisions).t = t;
        end;
        
        % if bottleneck ...
        if doDisturbance
            if VEH(i).x >= xDisturbance(1) && VEH(i).x < xDisturbance(2) ...
                    && t >= tDisturbance(1) && t < tDisturbance(2)
                VEH(i).par.T = 1 + disturbanceStrength; % inmpose a fictional bottleneck
            else
                VEH(i).par.T = 1.2;  % initial value
            end;
        end;
        
        if ~VEH(i).incident
            %% determine acceleration in IDM model (can be replaced with any
            % other model that uses ds and dv as stimuli; but extendeding these
            % stimuli with other ones is easy)
            VEH(i).a = IDMplus(VEH(i), ds, dv);
        
        else
            % that's the end of this vehicle ...
            VEH(i).v = 0;
            VEH(i).a = 0;
        end;
        
    end;
    
    %% ---- NOW THE ACTUAL MOVEMENT ---------------------------------------
    % we have all the information to determine the actual movement
    for i = iVEH
        
%         % only propagate if there is no incident ...
%         if ~VEH(i).incident
            
            % update speed and location: vehicle kinematics
            VEH(i).v = max(0, VEH(i).v + dt*VEH(i).a);
            VEH(i).x = VEH(i).x + max(0, VEH(i).v*dt + 0.5*VEH(i).a*dt^2);
            VEH(i).t = t;
            
            % new ground-truth GROSS gaps and headways (to compute densities
            % and flows)
            if i>iMin
                VEH(i).s = VEH(i-1).x - VEH(i).x;
                VEH(i).h = VEH(i).s/VEH(i).v;
            end;
            
            % save vehicle record:
%             VEH(i).traj(k,:) = [ ...
%                 VEH(i).t,VEH(i).x,VEH(i).v,VEH(i).a,VEH(i).s,VEH(i).h];
            VEH(i).traj = [VEH(i).traj; ...
                VEH(i).t,VEH(i).x,VEH(i).v,VEH(i).a,VEH(i).s,VEH(i).h];

            % save vehicle HF records:
            if doHF
%                 VEH(i).trajHF(k,:) = [ ...
%                     VEH(i).t,VEH(i).x,VEH(i).TDcf,VEH(i).TDdi, VEH(i).SA, ...
%                     VEH(i).ctl.tauextra + VEH(i).par.tau, ttc];
                VEH(i).trajHF = [VEH(i).trajHF; ...
                    VEH(i).t,VEH(i).x,VEH(i).TDcf,VEH(i).TDdi, VEH(i).SA, ...
                    VEH(i).ctl.tauextra + VEH(i).par.tau, ttc];
            end;
            
            % is vehicle finished or time is up ???
            VEH(i).finished = VEH(i).x >= xmax;
            VEH(i).kfin = k;

%         else
%             
%             % we remove the accident immediately
%             VEH(i).finished = true;
%             VEH(i).kfin = k;
%         end;
        
    end;
    
    % update progress
    if rem(k,dtprogress)==0 || k==1
        if k==1
            fprintf('Simulation %3.0f%%\n',100*k/numsteps);
        else
            fprintf('\b\b\b\b\b%3.0f%%\n',100*k/numsteps);
        end;
    end;
    % In the case there is a collision break the scenario!
    if IsIncident
         break;
     end;
    
end

%%  Finalise trajectories and then compute TTS and FD
 active   = [VEH.active];
 iVEH = find(active);
 TTS = 0; 
 V = []; S = []; H = [];
 TTC_crit = 0;
 RBR_crit = 0;
 mean_Speed_Edie = 0; 
 TDC = 0;
 crit_threshold = 4; % [s]
 crit_br_threshold = -2.0; % [m/s2]
% 
 for i=iVEH
% %     % first truncate trajectories
% %     VEH(i).traj = VEH(i).traj( VEH(i).kgen:VEH(i).kfin, :);
% %     if doHF
% %         VEH(i).trajHF = VEH(i).trajHF( VEH(i).kgen:VEH(i).kfin, :);
% %     end;
% %   then compute statistics
     TTS = TTS + (VEH(i).tqueue + VEH(i).t-VEH(i).tgen)/60;
     TDC = TDC + (VEH(i).traj(end,2) - VEH(i).traj(1,2))/1000;
     V = [V; VEH(i).traj(:,3)];
     S = [S; VEH(i).traj(:,5)];
     H = [H; VEH(i).traj(:,6)];
     % find id of vehicles
     idtcrit = VEH(i).trajHF(:,7) < crit_threshold;
     % calculate TTC crit
     TTC_crit = TTC_crit + sum(VEH(i).trajHF(idtcrit,7));
     % First find the id of the deceleration in matrix
      idRBR = VEH(i).traj(:,4) < crit_br_threshold ;
     % Secondly, compute RBR
     RBR_crit = RBR_crit + sum(VEH(i).traj(idRBR,4));
 end
% 
% % calculate mean speed for the whole segment
 mean_Speed_Edie = TDC / (TTS/60);
% 
TTC_crit = TTC_crit/3600;
% % we had one problem with the previous aparatus of TTS. The NUMBER of the vehicles is not the same under all the...
% % the simulations, therefore something Must to be done. Well, we normalize the TTS by dividin' it with the #vehicles.
 nActiveVehicles = length(find([VEH.active] == 1));
 TTS_bar = TTS / nActiveVehicles;
% fprintf('DONE\n'); 

%% Copy all updated vars back in Sim structure:
for i = 1:numel(flds)
    evalstr = sprintf('Sim.%s = %s;',flds{i}, flds{i});
    eval(evalstr);
end
toc;
end
