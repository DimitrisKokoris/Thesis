function [uoMat,actualProbability] = estimation(nveh,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Detailed Description goes here. 
% % This function returns the overestimation/underestimation capabilities
% % of the driving population.
% % parameters
% % p = 0.75 % Is the probability that the dirvers overestimation
% % 1 - p    % Is the probability that the drivers underestimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw a value from a uniform distribution.
% We use a seed here to guarantee that at each simulation scenaria, the
% same vehicles will have the same percpetion biases as the previous ones.
% In case we want more stochasticity, one can comment the following row.
rand('seed',1);
uoMat = rand([nveh,1]) <= p;
actualProbability = mean(uoMat);
% Now, the matrix is composed by 0 and 1 values, we need to substitute 0 to
% -1 for our simulation to work. Hence:
uoMat = sign(-0.1 + uoMat);
% More details of how this function works can be found to:
%  https://nl.mathworks.com/help/matlab/ref/sign.html
end