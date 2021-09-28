function inParPlot(dohf, dodist, doheterogeneity, dotau, doperc_ds, doperc_v, dov0, doT, doTD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Parameters that govern the conceptual framework,
% (1) inputs: 
%           doHF      ; Human factors
%           dodist    ; Impose flow-conserving bottleneck
%           doTD      ; Impose the concept of induced task difficulty
%           dotau     ; Reaction time
%           doperc_ds ; Relative distance perception errors
%           doperc_v  ; Speed perception error
%           dov0      ; Behavioural adaptation (desired speed)
%           doT       ; Behavioural adaptation (desired (time) headway)
%           doH       ; Heterogeneity in terms of (mental) capacity
% (2) output:
%          Plot the aforementiond parameters at this same order (T/F)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO HUMAN FACTORS AND IMPOSE A FLOW-CONSERVING DISTURBANCE
%     if dohf == true
%         fprintf('Human Factors: %s  \n', 'true');
%     else
%         fprintf('Human Factors: %s  \n', 'false');
%     end
%     if dodist == true
%         fprintf('Disturbance: %s  \n', 'true');
%     else
%         fprintf('Disturbance: %s  \n', 'false');
%     end
%     % INDUCED TASK DIFFICULTY
%         if doTD == true
%         fprintf('Task Difficulty: %s  \n', 'true');
%     else
%         fprintf('Task Difficulty: %s  \n', 'false');
%         end
%     % DO WE HAVE IMPACT ON REACTION TIME?
%     if dotau == true
%         fprintf('Reaction Time: %s  \n', 'true');
%     else
%         fprintf('Reaction Time: %s  \n', 'false');
%     end
%     % DO WE HAVE IMPACTS ON PERCEPTION??
%     if doperc_ds == true
%         fprintf('Distance perception: %s  \n', 'true');
%     else
%         fprintf('Distance perception: %s  \n', 'false');
%     end
%     
%     if doperc_v == true
%         fprintf('Speed perception: %s  \n', 'true');
%     else
%         fprintf('Speed perception: %s  \n', 'false');
%     end
%     % DO WE HAVE ANY BEHAVIOURAL ADAPTATION??
%     if dov0 == true
%         fprintf('Adaptation of Desired Speed: %s  \n', 'true');
%     else
%         fprintf('Adaptation of Desired Speed: %s  \n', 'false');
%     end
%     
%     if doT == true
%         fprintf('Adaptation of Desired Time Headway: %s  \n', 'true');
%     else
%         fprintf('Adaptation of Desired Time Headway: %s  \n', 'false');
%     end
%     
%     % DO WE HAVE ANY HETEROGENEITY OF THE DRIVING POPULATION???
%     if doheterogeneity == true
%         fprintf('Heterogeneity: %s  \n', 'true');
%     else
%         fprintf('Heterogeneity: %s  \n', 'false');
%     end
    %% Plot the parameters in the cmd in a table form
    T = table([dohf; dodist; doTD; dotau; doperc_ds; doperc_v; dov0; doT; doheterogeneity],...
        'VariableNames',{'State'},'RowName',{'Human Factors','Disturbance',...
        'Induced Task Difficulty','Reaction Time','Distance perception',...
        'Speed perception','Adaptation of Desired Speed','Adaptation of Desired Time Headway',...
        'Heterogeneity'});
    disp(T);
end