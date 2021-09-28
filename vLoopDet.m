function [Sim] = vLoopDet(Sim)
% Copy all vars that are in Sim structure to root shortcuts (e.g. Sim.dt
% becomes dt)
flds = fieldnames(Sim);
for i = 1:numel(flds)
    evalstr = sprintf('%s = Sim.%s;',flds{i}, flds{i});
    eval(evalstr);
end

flds{end+1} = 'Qvlp';
flds{end+1} = 'Vvlp';
flds{end+1} = 'Rhovlp';
flds{end+1} = 'detPos'; 
detPos      = 100 : 200 : xmax;  % Position we set the pseudo detectors
aggPeriod   = 30;                % Aggregation Period is 30s
period      = Tsim / aggPeriod;
nDet        = length(detPos);

%% Find all the active vehicles 
active = [VEH.active];     % which vehicles
iVEH   = find(active);     % vector of active vehicles
%% 
Qvlp(nDet, period)   = 0;
Vvlp(nDet, period)   = 0;
Rhovlp(nDet, period) = 0;

%% Do the Main Loop
% Expand over time
for i = 1 : period
    % Expand over space
    for j = 1 : nDet
        % Iterate over the vehicles' trajectories
        for k = 1 : length(iVEH)  % Number of vehivcle
            id  = find(VEH(k).traj(:,2)>=detPos(j));
            if any(id > 0)
                if (VEH(k).traj(id(1),1) >= (i-1) * aggPeriod && VEH(k).traj(id(1),1) <= i * aggPeriod)
                    % Flow for the period i is:
                    Qvlp(j,i) = Qvlp(j,i) + 1;
                    % Interpolate linearly for obtaining the speed!
                    v = interp1([VEH(k).traj(id(1) - 1,2), VEH(k).traj(id(1),2)],[ VEH(k).traj(id(1) - 1,3), VEH(k).traj(id(1),3)],detPos(j));
                    % Estimate Harmonic Speed!!!
                    Vvlp(j,i) = Vvlp(j,i) + 1 / ( v * 3.6 ); % Harmonic speed
                end
            end
        
        end 
    end   
end
%% Final Calculations
% https://en.wikipedia.org/wiki/Harmonic_mean
         Vvlp   = Qvlp ./ Vvlp;
         Qvlp   = Qvlp .* 3600  ./ aggPeriod ;
         Rhovlp = Qvlp ./ Vvlp;


%% Copy all updated vars back in Sim structure:
for i = 1:numel(flds)
    evalstr = sprintf('Sim.%s = %s;',flds{i}, flds{i});
    eval(evalstr);
end

end
