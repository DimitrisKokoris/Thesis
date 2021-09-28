function [Sim] = virtDet(Sim)
% Copy all vars that are in Sim structure to root shortcuts (e.g. Sim.dt
% becomes dt)
flds = fieldnames(Sim);
for i = 1:numel(flds)
    evalstr = sprintf('%s = Sim.%s;',flds{i}, flds{i});
    eval(evalstr);
end

flds{end+1} = 'virtDn';
flds{end+1} = 'virtDv';
flds{end+1} = 'virtDrho';
flds{end+1} = 'detPos';
detPos = 100 : 500 : xmax;
aggPeriod   = 10; % Aggregation Period is 1 min(10s)
period      = Tsim / aggPeriod;

%% finalise trajectories 
active = [VEH.active];     % which vehicles
iVEH   = find(active);     % vector of acitve vehicles
virtDn = [];
virtDv = [];

for dp = 1:length(detPos)                                    % Detector Position
    for i = 1 : period                                      % Aggregation Period
         virtDn(dp,i) = 0;                                   % # of vehicles passing
         virtDv(dp,i) = 0;                                  % Speeds
        for j = 1 : length(iVEH)                            % Number of vehivcle
                id = find(VEH(j).traj(:,2) >= detPos(dp));  % Find all the id'
                if any(id > 0)
                    if VEH(j).traj(id(1),1) <= aggPeriod * i          %
                        virtDn(dp,i) = virtDn(dp,i) + 1;
                        virtDv(dp,i) = virtDv(dp,i) + 1/(VEH(j).traj(id(1),3)*3.6); % Harmonic speed
                    end
                end
        end
        virtDv(dp,i) = virtDn(dp,i)  / virtDv(dp,i);
        virtDn(dp,i) = virtDn(dp,i) /(i * aggPeriod) * 3600 ;
        virtDrho(dp,i) = virtDn(dp,i) ./ virtDv(dp,i);  
    end
end



%% Copy all updated vars back in Sim structure:
for i = 1:numel(flds)
    evalstr = sprintf('Sim.%s = %s;',flds{i}, flds{i});
    eval(evalstr);
end

end
