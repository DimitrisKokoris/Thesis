function [figs] = visualisation_fcn(Sim, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % [figs] = visualisation_fcn(Sim):
% % visualise results for Tr-B paper and get handles to figures in return
%
% % A generic multi-level framework for microscopic traffic simulation?theory
% % and an example case in modelling driver distraction  
% % J.W.C. van Lint* & S.C. Calvert
% % Submitted to Transportation Rearch Part B, 2018
%
% % Copyright Hans van Lint, 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first argument is a structure or a filename
if ~isstruct(Sim)
    Sim = load(Sim);
end;

%% Copy all vars are in Sim structure:
flds = fieldnames(Sim);
for i = 1:numel(flds)
    evalstr = sprintf('%s = Sim.%s;',flds{i}, flds{i});
    eval(evalstr);
end

q_cap(1:length(q)) = 2000; 

%% missing vars?
if isempty(whos('xmax'))
    xmax = 6e3;
end;

% 'overwrite' vars? (requires a minimum of 2 extra (par+val) arguments)
if nargin>3
    for j = 1:2:(nargin-2)
        if iscell(varargin{j+1})
            varargin{j+1} = sprintf('''%s''',varargin{j+1}{1});
        end
        evalstr = sprintf('%s = %s;',varargin{j},varargin{j+1});
        eval(evalstr);
    end;
end;

%% plot HF functions
disp('FD for HF and demand ...');
fsz = 14; %'fontsize'

veh = VEH(1); 
hh = 0:.1:5; xxdist = -500:50:450; ts = 0:.1:2.5;
veh = computeHF(veh,hh,xxdist,ts);

figFDTD = figure; %figFDTD.WindowStyle = 'docked'; 
figFDTD.Name = sprintf('TDFD %s',scenario);
figFDTD.Position = [50,900,900,300];
cla; 

if doHeterogeneity
    % we 
    active   = [VEH.active];
    iVEH = find(active);
    numveh = numel(iVEH);
    for i = iVEH
        if any (iHF==i)
            c = 'b'; w = 1;
        else
            c = 'k'; w = .1;
        end;
        veh = VEH(i);
        veh = computeHF(veh,hh,xxdist,ts);
        
        % car following
        subplot(1,2,1);
        if i ==1
            hold on;
        end;
        plot(hh,veh.TDcf,'color',c,'linewidth',w)
        if i ==numveh
            hold off;
            set(gca,'xlim', [0,5], 'ylim', [0, 1.2], ...
                'box','on', 'xgrid','on','ygrid','on');
            xlabel('headway  (s) \rightarrow');
            ylabel('task demand (-) \rightarrow');
            title('(a) FD task demand car following');
        end;
        
        % distraction
%         subplot(1,2,2);
%         if i ==1
%             hold on;
%         end;
%         plot(xxdist,veh.TDdi,'color',c,'linewidth',w)
%         if i ==numveh
%             hold off;
%             set(gca,'xlim', [-550,450], 'ylim', [0, 1.2], ...
%                 'box','on', 'xgrid','on', 'ygrid','on');
%             xlabel('distance to distraction (m) \rightarrow');
%             ylabel('task demand (-) \rightarrow');
%             title('(b) FD task demand distraction');
%         end;
        % awareness
        subplot(1,3,3);
        if i ==1
            hold on;
        end;
        plot(veh.TS,veh.SA,'color','k','linewidth',.1)
        set(gca,'xlim', [0,2.5], 'ylim', [0, 1.2], ...
            'box','on','xgrid','on', 'ygrid','on');
        if i ==numveh
            hold off;
            xlabel('task saturation (-) \rightarrow');
            ylabel('awareness (-) \rightarrow');
            title('(c) Awareness (level 2: understanding)');
        end;
    end;
    
else
    % car following
    subplot(1,3,1);
    plot(hh,veh.TDcf,'color','k','linewidth',2)
    set(gca,'xlim', [0,5], 'ylim', [0, 1.2], ...
        'box','on', 'xgrid','on','ygrid','on');
    xlabel('headway  (s) \rightarrow');
    ylabel('task demand (-) \rightarrow');
    title('FD task demand car following');
%     % distraction
%     subplot(1,3,2);
%     plot(xxdist,veh.TDdi,'color','k','linewidth',2)
%     set(gca,'xlim', [-550,450], 'ylim', [0, 1.2], ...
%         'box','on', 'xgrid','on', 'ygrid','on');
%     xlabel('distance to distraction (m) \rightarrow');
%     ylabel('task demand (-) \rightarrow');
%     title('FD task demand distraction');
    % awareness
    subplot(1,3,3);
    plot(veh.TS,veh.SA,'color','k','linewidth',2)
    set(gca,'xlim', [0,2.5], 'ylim', [0, 1.2], ...
        'box','on','xgrid','on', 'ygrid','on');
    xlabel('task saturation (-) \rightarrow');
    ylabel('awareness (-) \rightarrow');
    title('Awareness (level 2: understanding)');
end;

%% demand profile
figdemand = figure; %figdemand.WindowStyle = 'docked';
figdemand.Name = sprintf('Demand %s',scenario);
figdemand.Position = [50,900,500,500];
cla; 
hold on;
plot(timeline,q,'color','k','linewidth',2)
plot(timeline, q_cap,'color','r','linewidth',2)
set(gca,'xlim', [0,Tsim(end)], 'ylim', [0, 2500], ...
    'box','on', 'ygrid','on', 'xgrid', 'on', 'clim',[0,140/3.6]);
legend({'Demand Pattern','Capacity of Bottleneck'},'location','northeast');
xlabel('time  (s) \rightarrow');
ylabel('demand (veh/h) \rightarrow');
%colormap jet; colorbar;
title(sprintf('Demand pattern %s',scenario));

%% FD plots (original version)
figFD = figure;
figFD.Name = sprintf('FD %s',scenario);
figFD.Position = [50,50,900,300];
% plot(S,V,'linestyle','none','marker','*');
% set(gca, 'xlim',[0,500],'ylim',[0,40]);
% xlabel('Distance gap (m) \rightarrow');
% ylabel('Headway (s) \rightarrow');
subplot(1,2,1);
plot(1e3./S,3600./H,'linestyle','none','marker','.','color','k');
set(gca, 'xlim',[0,150],'ylim',[0,3000],'box','on', 'xgrid','on','ygrid','on');
xlabel('Density (veh/km) \rightarrow');
ylabel('Flow (veh/h) \rightarrow');
title('Q-K relation')

subplot(1,2,2);
plot(1e3./S,3.6*V,'linestyle','none','marker','.','color','k');
set(gca, 'xlim',[0,150],'ylim',[0,135],'box','on', 'xgrid','on','ygrid','on');
xlabel('Density (veh/km) \rightarrow');
ylabel('Speed (km/h) \rightarrow');
title('V-K relation')

%% FD for just selected vehicles

figFDsel = figure;
figFDsel.Name = sprintf('FD Selection %s',scenario);
figFDsel.Position = [50,550,900,300];
cols = 'kbrm';
% plot(S,V,'linestyle','none','marker','*');
% set(gca, 'xlim',[0,500],'ylim',[0,40]);
% xlabel('Distance gap (m) \rightarrow');
% ylabel('Headway (s) \rightarrow');

i=1; % we use the params of VEH 1
ii = 1; % color for plotting

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
    'linewidth',1,'color',cols(ii));
subplot(1,2,2);
hold on;
plot(rr, vv,'linestyle',':','marker','none',...
    'linewidth',1,'color',cols(ii));

for i=iHF
    
    % get speed, gap and headway
    v = VEH(i).traj(:,3);
    s = VEH(i).traj(:,5);
    h = VEH(i).traj(:,6);

    subplot(1,2,1);
    hold on;
    h1(ii) = plot(1e3./s,3600./h,'linestyle','-','marker','.',...
        'linewidth',1,'color',cols(ii),...
        'displayname',sprintf('Driver %d',i));

    subplot(1,2,2);
    hold on;
    h2(ii)=plot(1e3./s,3.6*v,'linestyle','-','marker','.',...
        'linewidth',1,'color',cols(ii),...
        'displayname',sprintf('Driver %d',i));
    %'markerfacecolor',cols(ii),
    ii = ii+1;
end;

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


%% trajectory plotting
disp('Trajectory figure ...');
figtraj = figure; %figtraj.WindowStyle = 'docked'; 
figtraj.Name = sprintf('Trajectories %s',scenario);
figtraj.Position = [50,50,1280,520];
cla; axSim = gca;
set(axSim,'xlim', [0,Tsim(end)], 'ylim', [xstart, xmax], ...
    'box','on', 'clim',[0,140/3.6],'fontsize',fsz);
xlabel('time  (s) \rightarrow');
ylabel('space (m) \rightarrow');
%colormap jet; colorbar;
title(sprintf('Trajectories %s',scenario));

% run through all vehicles
veh = createDefVehicle(Sim);
active   = [VEH.active];
iVEH = find(active);
for i = iVEH
    if any(iHF == i)
       VEH(i).color = 'red'; VEH(i).width =1;
       color_line(Sim.VEH(i).traj(:,1),Sim.VEH(i).traj(:,2),Sim.VEH(i).traj(:,3),VEH(i).color);
    else
       color_line(Sim.VEH(i).traj(:,1),Sim.VEH(i).traj(:,2),Sim.VEH(i).traj(:,3));
       % VEH(i).color = veh.color; VEH(i).width =veh.width;
    end;
    
    % VEH(i).hobj = line(axSim, VEH(i).traj(:,1), VEH(i).traj(:,2), ...
    %     VEH(i).traj(:,3),'color',VEH(i).color,'linewidth',VEH(i).width);
   
    % Here put a marker for accident!!!
    if VEH(i).incident
         hobj= line(VEH(i).inc.t, VEH(i).inc.x, 'markersize', 12, ...
            'Marker','d','markerfacecolor','r','markeredgecolor','r');
    end;
end
  c = jet;
  c = flipud(c);
  colormap(c);
  colorbar;
  hcb = colorbar;
  hcb.Title.String = "m/s";
% Color specific vehicle
for i = iVEH
  if any(iHF == i)
      
  end 
end

% add TTS annotation
figure(figtraj);
hTxt = annotation('textbox');
set(hTxt,'position',[.15,.7,.3,.2],'FitBoxToText','on', ...
    'String',sprintf('TTS: %.3f (h)',Edie.TTS), 'fontsize',18, ...
    'backgroundcolor',[1,1,1], 'facealpha',0.9);
% draw incident / BN location
% line(axSim,[0,Tsim],[xDistraction,xDistraction],'color','k', 'linestyle','--','linewidth',2);
% and update title
title(axSim, sprintf('(a) Simulation (%.0f/%.0f)',t,Tsim));
drawnow;

%% selected HF trajectories
disp('HF figures ...');

% figure for HF trajectories
figHF = figure; %figHF.WindowStyle = 'docked'; 
figHF.Name = sprintf('HF dynamics %s',scenario);
figHF.Position = [50,900,1280,520];
cla; 

numveh = numel(iHF);
numsubplots = 3*numveh;
captitle ='bcdefghijklmnopqrstuvwxyz';
cols = 'kbrgmc';
icol = 1;

kk =1;
for k = 1:numveh
    % vehicle id
    i = iHF(k);
    
    % create plot indices
    plotnr = k:numveh:numsubplots;
    trajHF = VEH(i).trajHF;
    traj   = VEH(i).traj;
    
    % first plot: task demand
    subplot(3,numveh,plotnr(1)); cla;
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
    
    
    title(sprintf('(%s) Task demand & Awareness (veh %d)',...
      captitle(kk), i));
    legend({'TD_{cf}','TD_{tot}','TS_{cr}','SA'},...
        'location','eastoutside','orientation','vertical');    
    kk=kk+1;
    
    % second plot: speed and acceleration
    subplot(3,numveh,plotnr(2)); cla;
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

    title(sprintf('(%s) Speed & acceleration (veh %d)',...
        captitle(kk),i));
    legend({'Spd','Acc'},...
        'location','eastoutside','orientation','vertical');    
    kk=kk+1;
    
    % third plot: reaction time and TTC
    subplot(3,numveh,plotnr(3)); cla;
    set(gca,'box','on','ygrid','on','fontsize',fsz);

    yyaxis left
    plot(trajHF(:,2), trajHF(:,6), ... %'color', cols(c), ...
        'linestyle','-','linewidth',2);
    set(gca, 'xlim', [xstart, xmax],'ylim',[0,1.5]);
    ylabel('\tau (s) \rightarrow');
    
    yyaxis right
    plot(trajHF(:,2), trajHF(:,7), ... %'color', cols(c), ...
        'linestyle','-','linewidth',2);
    set(gca, 'xlim', [xstart, xmax],'ylim',[0,15]);
    ylabel('TTC (s) \rightarrow');

    title(sprintf('(%s) Reaction time & TTC (veh %d)',...
        captitle(kk), i));
    xlabel('Distance (m) \rightarrow');
    legend({'Tau','TTC'},...
        'location','eastoutside','orientation','vertical');    
    kk=kk+1;
    
end;

%% Plot headway, acceleration histogram for all vehicles in simulation
disp('Headway and Acceleration Distribution ...');
figacchead = figure; %figtraj.WindowStyle = 'docked'; 
figacchead.Name = sprintf('Acceleartion and Headway Distribution %s',scenario);
figacchead.Position = [50,50,1280,520];
title(sprintf('Acceleartion and Headway Distribution %s',scenario));

headways = [];
% find the active vehicles
active   = [VEH.active]; % find active vehicles
iVEH = find(active);     % number of active vechicles
for i = iVEH
    headways = vertcat(headways,VEH(i).traj(:,6)); % time heaways
end
% rmv Inf headways
% thesehead   = find(headways == Inf);
thesehead10 = find(headways >= 25);
% consider only headways <= 10 [s]
headwayfiltered = headways;
headwayfiltered(thesehead10,:)= [];
% Fist plot: Histogramm of time headways.
subplot(1,3,1)
edges = 0:0.2:25;
[n, xout] = hist(headwayfiltered,edges);
bar(xout, n, 'barwidth', 1, 'basevalue', 1);
set(gca, 'xlim', [0, 25]);
set(gca,'YScale','log','box','on', 'xgrid','on','ygrid','on');
ylabel('Frequency');
xlabel('Headway (s)');

acc = [];
% find the active vehicles
active   = [VEH.active]; % find active vehicles
iVEH = find(active);     % number of active vechicles
for i = iVEH
    acc = vertcat(acc,VEH(i).traj(:,4)); % acc
end
% Second subplot: Histogram of aceeleraitons.
subplot(1,3,2)
edges = -8:0.2:8;
[n, xout] = hist(acc,edges);
bar(xout, n, 'barwidth', 1, 'basevalue', 1);
set(gca,'YScale','log','box','on', 'xgrid','on','ygrid','on');
set(gca, 'xlim', [-8.0, 8]);
ylabel('Frequency');
xlabel('Acceleration (m/s^2)');

sp = [];
% find the active vehicles
active   = [VEH.active]; % find active vehicles
iVEH = find(active);     % number of active vechicles
for i = iVEH
    sp = vertcat(sp,VEH(i).traj(:,3)); % speed
end
% Third subplot : Histogram of Speeds.
subplot(1,3,3)
edges = 0:1:40;
[n, xout] = hist(sp,edges);
bar(xout, n, 'barwidth', 1, 'basevalue', 1);
set(gca,'YScale','log','box','on', 'xgrid','on','ygrid','on');
set(gca, 'xlim', [0,40]);
ylabel('Frequency');
xlabel('Speed (m/s)');

%% Plot Pseudo -Perfect- Detectors Colleted data
disp('Preudo -perfect- Detectors Data ...');
figdet= figure; %figtraj.WindowStyle = 'docked'; 
figdet.Name = sprintf('Preudo-Perfect-Detectors %s',scenario);
figdet.Position = [50,50,1366,512];
title(sprintf('Preudo-Perfect-Detectors %s',scenario));

i=1; % we use the params of VEH 1
ii = 1; % color for plotting

% compute equilibrium coditions FD
rmax = 1e3/VEH(i).par.s0; % jam density
qmax = 3.6e3/(VEH(i).par.T*VEH(i).ctl.Tfac); % max flow
v0=VEH(i).ctl.v0fac*VEH(i).par.v0*3.6; % max speed
rr = 0:rmax;
vv = min(v0, qmax*(1./rr-1./rmax));
qq = vv.*rr;

vv0(1:length(rr)) = v0; 

subplot(1,3,1);
hold on;
plot(rr, qq, 'linestyle',':','marker','none',...
    'linewidth',1,'color',cols(ii));

subplot(1,3,2);
hold on;
plot(rr, vv,'linestyle',':','marker','none',...
    'linewidth',1,'color',cols(ii));

subplot(1,3,3);
hold on;
plot(qq, vv0,'linestyle',':','marker','none',...
    'linewidth',1,'color',cols(ii));
plot(qq, vv,'linestyle',':','marker','none',...
    'linewidth',1,'color',cols(ii));

% Flow - Density
subplot(1,3,1);
hold on;
plot(Rhovlp,Qvlp,'linestyle','none','marker','.','color','r')
set(gca, 'xlim',[0,130],'ylim',[0,3000],'box','on', 'xgrid','on','ygrid','on');
xlabel('Density (veh/km) \rightarrow');
ylabel('Flow (veh/h) \rightarrow');
title('Q-K relation');

% Density - Speed
subplot(1,3,2);
hold on;
plot(Rhovlp, Vvlp,'linestyle','none','marker','.','color','r')
set(gca,'xlim',[0,130],'ylim',[0,135],'box','on', 'xgrid','on','ygrid','on');
xlabel('Density (veh/km) \rightarrow');
ylabel('Speed (km/h) \rightarrow');
title('V-K relation');

%Flow - speed
subplot(1,3,3);
hold on;
plot(Qvlp, Vvlp,'linestyle','none','marker','.','color','r')
set(gca,'xlim',[0,3000],'ylim',[0,135],'box','on', 'xgrid','on','ygrid','on');
xlabel('Flow (veh/h) \rightarrow');
ylabel('Speed (km/h) \rightarrow');
title('V-Q relation');

% output
figs.FDTD   = figFDTD;
figs.Demand = figdemand;
figs.FD     = figFD;
figs.FDsel  = figFDsel;
figs.Traj   = figtraj;
figs.HF     = figHF;
figs.HEADACC= figacchead;
figs.PPD    = figdet;

end

