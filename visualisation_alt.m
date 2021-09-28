function [fig, figFDsel] = visualisation_alt(SimStruct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [fig,fig2] = visualisation_alt(Sim):
% visualise results for Tr-B paper and get handle to figures in return
%
% A generic multi-level framework for microscopic traffic simulation?theory
% and an example case in modelling driver distraction  
% J.W.C. van Lint* & S.C. Calvert
% Submitted to Transportation Rearch Part B, 2018
%
% Copyright Hans van Lint, 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a ='abcdefghijklmnopqrstuvwxyz';
numSim = numel(SimStruct);
fsz = 14;
k=1;

fig = figure;
set(fig,'pos',[285 114 1920 1040]);
for iSim=1:numSim
    
    % select Sim structure
    Sim=SimStruct(iSim);
    
    % Copy all vars are in Sim structure:
    flds = fieldnames(Sim);
    for i = 1:numel(flds)
        evalstr = sprintf('%s = Sim.%s;',flds{i}, flds{i});
        eval(evalstr);
    end;
    
    % trajectories
    pl = [iSim:numSim:numSim*4];
    
    subplot(4,numSim,pl(1));
    
    % trajectories
    set(gca,'xlim', [0,Tsim(end)], 'ylim', [xstart, xmax], ...
        'box','on', 'clim',[0,140/3.6],'fontsize',fsz);
    xlabel('time  (s) \rightarrow');
    ylabel('space (m) \rightarrow');
    %colormap jet; colorbar;
    title(sprintf('(%s) %s',a(k), Title));
    k=k+1;
    
    % run through all vehicles
    veh = createDefVehicle;
    active   = [VEH.active];
    iVEH = find(active);
    for i = iVEH
        if any(iHF == i)
            VEH(i).color = 'b'; VEH(i).width =1;
        else
            VEH(i).color = veh.color; VEH(i).width =veh.width;
        end;
        VEH(i).hobj = line(gca, VEH(i).traj(:,1), VEH(i).traj(:,2), ...
            VEH(i).traj(:,3),'color',VEH(i).color,'linewidth',VEH(i).width);
        if VEH(i).incident
            hobj= line(VEH(i).inc.t, VEH(i).inc.x, 'markersize', 12, ...
                'Marker','d','markerfacecolor','r','markeredgecolor','r');
        end;
    end;
    
    % draw incident / BN location
    line(gca,[0,Tsim],[xDistraction,xDistraction],'color','k', 'linestyle','--','linewidth',2);
    
    % HF plots
    cols = 'kbrgmc';
    icol = 1;
    
    i = iHF(1);
    
    % trajectories
    trajHF = VEH(i).trajHF;
    traj   = VEH(i).traj;
    
    % first plot: task demand
    subplot(4,numSim,pl(2));
    
    yyaxis right
    plot(trajHF(:,2), trajHF(:,5), ... %'color', cols(icol), ...
        'linestyle','-','linewidth',2);
    set(gca, 'ylim', [0.4,1.1], 'xlim', [xstart, xmax]);
    ylabel('SA (-) \rightarrow');
    
    yyaxis left
    line(trajHF(:,2), trajHF(:,3), ...%'color', cols(icol), ...
        'linestyle','-','linewidth',2);
    line(trajHF(:,2), trajHF(:,3)+trajHF(:,4), 'color', cols(icol), ...
        'linestyle','-','linewidth',2);
    line([xstart,xmax], [VEH(i).SAFn.tscrit,VEH(i).SAFn.tscrit], ...%'color', 'k', ...
        'linestyle',':','linewidth',2);
    set(gca, 'ylim', [0,1.7], 'xlim', [xstart, xmax], ...
        'box','on','ygrid','on','fontsize',fsz);
    ylabel('TD (-) \rightarrow');
    
    title(sprintf('(%s) Task demand & Awareness (veh %d) ',a(k), i));
    k=k+1;
    
    legend({'TD_{cf}','TD_{tot}','TS_{cr}','SA'},...
        'location','eastoutside','orientation','vertical');
    
    % second plot: speed and acceleration
    subplot(4,numSim,pl(3));
    
    set(gca, 'box','on','ygrid','on','fontsize',fsz);
    
    yyaxis left
    plot(traj(:,2), traj(:,3), ... %'color', cols(icol), ...
        'linestyle','-','linewidth',2);
    set(gca, 'xlim', [xstart, xmax],'ylim',[0,40]);
    ylabel('Spd (m/s) \rightarrow');
    
    yyaxis right
    plot(traj(:,2), traj(:,4), ... %'color', cols(icol), ...
        'linestyle','-','linewidth',2);
    set(gca, 'xlim', [xstart, xmax],'ylim',[-8,8]);
    ylabel('Acc (m/s^2) \rightarrow');
    
    title(sprintf('(%s) Speed & acceleration (veh %d) ',a(k), i));
    k=k+1;
    legend({'Spd','Acc'},...
        'location','eastoutside','orientation','vertical');
    
    % third plot: reaction time and TTC
    subplot(4,numSim,pl(4));
    set(gca,'box','on','ygrid','on','fontsize',fsz);
    
    yyaxis left
    plot(trajHF(:,2), trajHF(:,6), ... %'color', cols(c), ...
        'linestyle','-','linewidth',2);
    set(gca, 'xlim', [xstart, xmax],'ylim',[0,2]);
    ylabel('\tau (s) \rightarrow');
    
    yyaxis right
    plot(trajHF(:,2), trajHF(:,7), ... %'color', cols(c), ...
        'linestyle','-','linewidth',2);
    set(gca, 'xlim', [xstart, xmax],'ylim',[0,15]);
    ylabel('TTC (s) \rightarrow');
    
    title(sprintf('(%s) Reaction time & TTC (veh %d)',a(k), i));
    k=k+1;
    xlabel('Distance (m) \rightarrow');
    legend({'Tau','TTC'},...
        'location','eastoutside','orientation','vertical');
    
end;

figFDsel = figure;
figFDsel.Name = sprintf('FD Selection %s',scenario);
figFDsel.Position = [50,550,900,300];

for iSim=1:numSim
    
    cols = 'kbrm';
    mrkr = '*.sd+';
    
    % select Sim structure
    Sim=SimStruct(iSim);
    
    % Copy all vars are in Sim structure:
    flds = fieldnames(Sim);
    for i = 1:numel(flds)
        evalstr = sprintf('%s = Sim.%s;',flds{i}, flds{i});
        eval(evalstr);
    end;
    
    % plot(S,V,'linestyle','none','marker','*');
    % set(gca, 'xlim',[0,500],'ylim',[0,40]);
    % xlabel('Distance gap (m) \rightarrow');
    % ylabel('Headway (s) \rightarrow');
    
    i=1; % we use the params of VEH 1
    
    % compute eq FD
    rmax = 1e3/VEH(i).par.s0;
    qmax = 3.6e3/(VEH(i).par.T*VEH(i).ctl.Tfac);
    v0=VEH(i).ctl.v0fac*VEH(i).par.v0*3.6;
    rr = 0:rmax;
    vv = min(v0, qmax*(1./rr-1./rmax));
    qq = vv.*rr;
    
    subplot(1,2,1);
    hold on;
    plot(rr, qq, 'linestyle',':','marker','none',...
        'linewidth',1,'color',cols(iSim));
    subplot(1,2,2);
    hold on;
    plot(rr, vv,'linestyle',':','marker','none',...
        'linewidth',1,'color',cols(iSim));
    
    i = iHF(1);
    
    % get speed, gap and headway
    v = VEH(i).traj(:,3);
    s = VEH(i).traj(:,5);
    h = VEH(i).traj(:,6);
    
    subplot(1,2,1);
    hold on;
    h1(iSim) = plot(1e3./s,3600./h,'linestyle','-','marker',mrkr(iSim),...
        'linewidth',.25,'color',cols(iSim),...
        'displayname',sprintf('Driver %d',i));
    
    subplot(1,2,2);
    hold on;
    h2(iSim)=plot(1e3./s,3.6*v,'linestyle','-','marker',mrkr(iSim),...
        'linewidth',.25,'color',cols(iSim),...
        'displayname',sprintf('Driver %d',i));
    %'markerfacecolor',cols(ii),
    
    if iSim==numSim
        subplot(1,2,1);
        hold off;
        set(gca, 'xlim',[0,rmax+10],'ylim',[0,max(qq)+100],...
            'box','on', 'xgrid','on','ygrid','on','fontsize',12);
        xlabel('Density (veh/km) \rightarrow');
        ylabel('Flow (veh/h) \rightarrow');
        title('Q^e(\rho)')
        legend(h1);
        
        subplot(1,2,2);
        hold off;
        set(gca, 'xlim',[0,rmax+10],'ylim',[0,v0+10],...
            'box','on', 'xgrid','on','ygrid','on','fontsize',12);
        xlabel('Density (veh/km) \rightarrow');
        ylabel('Speed (km/h) \rightarrow');
        title('V^e(\rho)')
        legend(h2);
    end;
end

end

