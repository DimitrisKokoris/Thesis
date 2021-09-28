function Sim = Edie(Sim)
%% Copy all vars that are in Sim structure to root shortcuts (e.g. Sim.dt
% becomes dt)
flds = fieldnames(Sim);
for i = 1:numel(flds)
    evalstr = sprintf('%s = Sim.%s;',flds{i}, flds{i});
    eval(evalstr);
end
% Copy all vars that are in Sim structure to root shortcuts (e.g. Sim.dt
% becomes dt)
flds = fieldnames(Sim);
for i = 1:numel(flds)
    evalstr = sprintf('%s = Sim.%s;',flds{i}, flds{i});
    eval(evalstr);
end
% add fields that are generated in this function
flds{end+1} = 'Edie';


%% The Actual code goes here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % We want to find the capacity of a predefined @area@ by using Edie's
% % Generalized definitions for the fundamental parameters (q,k,u). 
% % Thus, we defined an area on a t-x plane with boundaries:
% % tmin = 150s
% % tmax = 800s
% % xmin = 1000m
% % xmax = 5000m
% % The equation for calculating capacity is:
% % q = d(wmega) / |wmega|, wmega is the area of calculating the fundamental parameter
% % k = t(wmega) / |wmega|, -||-
% % v = q / k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finalise trajectories and then compute q,k,u
active   = [VEH.active];
iVEH = find(active);
TTS_EDIE = 0; TDC_EDIE = 0;
% Dedine the Area ?
tmine = 100;
tmaxe = 750;
xmine = 500;
xmaxe = 4500; % this is bad programming, change the variable name, interfere with xmax, the main parameter of simulation
wmega = (tmaxe - tmine) * (xmaxe - xmine);
wmega = wmega / 3600000;
for i=iVEH
    % Check the trajectorie if its is inside the area ?
    %    if (any(VEH(i).traj(:,1) >= tmin) && any(VEH(i).traj(:,1) <= tmax)...
    %      && any(VEH(i).traj(:,2) >= xmin) && any(VEH(i).traj(:,2) <= xmax))
    % filter
       st = [ VEH(i).traj(:,1), VEH(i).traj(:,2)];
       these_xmin = (st(:,2) >= xmine);
       st = st(these_xmin,:);
       these_tmin = (st(:,1) >= tmine);
       st = st(these_tmin,:);
       these_xmax = (st(:,2) <= xmaxe);
       st = st(these_xmax, :);
       these_tmax = (st(:,1) <= tmaxe);
       st = st(these_tmax,:);
       if ~isnan(st)
       TTS_EDIE = TTS_EDIE + st(end,1) - st(1,1);
       TDC_EDIE = TDC_EDIE + st(end,2) - st(1,2);
       end
       st = [];
  %  end
end
Edie.st = st;
%% Consistency of units
TDC_EDIE = TDC_EDIE / 1000;
TTS_EDIE = TTS_EDIE/3600;
Edie.TTS = TTS_EDIE;
Edie.TDC = TDC_EDIE;
Edie.q = TDC_EDIE / wmega;
Edie.k = TTS_EDIE / wmega;
Edie.u = TDC_EDIE / TTS_EDIE;

%% Copy all updated vars back in Sim structure:
for i = 1:numel(flds)
    evalstr = sprintf('Sim.%s = %s;',flds{i}, flds{i});
    eval(evalstr);
end
return;

%% Additional code TO DRAW AREA OMEGA
x = [100 750 750 100];
y = [500 500 4500 4500];
patch(x,y,'yellow', 'FaceAlpha',0.2 );
txt = '|Ω|';
text(700,1000,txt,'Fontsize', 22);

end
