function  xlsPlot(scenarios)
 % FileList = dir(fullfile(uigetdir(),'**','*.mat')); % Wrap the directory
 numscen = length(scenarios); % get size
 %numscen = numscen(1,1); % size
 
 %% Compute TTS and num inc table
taus = 5;
tbl_TTS = zeros(length(taus), numscen);
tbl_INC = zeros(length(taus), numscen);
tbl_TDC = zeros(length(taus), numscen); % density
tbl_SPD = zeros(length(taus), numscen); % speed
tbl_RBR = zeros(length(taus), numscen); % speed
for itau = 1:numel(taus)
    aTaustr = sprintf('tau=%.2d',taus(itau));
    for iscenario = 1:numscen
        aScenario = scenarios{iscenario};
        fname = sprintf('RainIntResults/%s_%s.mat',aScenario, aTaustr);
        load(fname,'Sim');
        tbl_TTS(itau, iscenario) = Sim.TTS;
        tbl_INC(itau, iscenario) = Sim.numCollisions;
        tbl_TDC(itau, iscenario) = Sim.TDC;
        tbl_SPD(itau, iscenario) = Sim.mean_Speed_Edie;
        tbl_RBR(itau, iscenario) = Sim.RBR_crit;     
    end
end

%% Pops up gui to manually choose the outputdir and name of file.xls
[filename, pathname] = uiputfile('*.xls', 'Choose a file name');
outname = fullfile(pathname, filename);
% % Write data to specific sheets in excel file!
xlswrite(outname,scenarios,'A1');
xlswrite(outname,[tbl_TTS; tbl_INC; tbl_TDC; tbl_SPD; tbl_RBR],'A2');
% xlswrite(outname, M);
 
end