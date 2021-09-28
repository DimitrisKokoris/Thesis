function ttc = TTC(Sim)
% Copy all vars that are in Sim structure to root shortcuts (e.g. Sim.dt
% becomes dt)
flds = fieldnames(Sim);
for i = 1:numel(flds)
    evalstr = sprintf('%s = Sim.%s;',flds{i}, flds{i});
    eval(evalstr);
end
ttc = 0;
%%  Finalise trajectories and then compute TTS and FD
 active   = [VEH.active];
 iVEH = find(active);
 
 for i=iVEH
     
     
 ttc = ttc + sum(~isnan(VEH(i).trajHF(:,7)));
 
 end

ttc = ttc / (60);
end