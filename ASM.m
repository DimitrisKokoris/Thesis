function [z, u, x, t, w, zcong, zfree, ucong, ufree, Ncong, Nfree] = ...
    ASM(Z,U,X,T,varargin)

% spatiotemporal traffic filter on the basis of Treiber and Helbing (2002)
% 
% function [z, u, x, t, w, zcong, zfree, ucong, ufree, Ncong, Nfree] ...
%           = ASM(Z,U,X,T)
%       OR
% function [z, u, x, t, w, zcong, zfree, ucong, ufree, Ncong, Nfree] ...
%           = rl_asm(Z,U,X,T,'parameter','value', ...)
%
% inputs   %%[numDets: number of detectors]
%     Z:   numDets x numTimeSteps matrix of traffic observations (speeds[km/h],
%          flows[veh/h], etc). Z may contain NaN values (missing or corrupted data). 
%
%     U:   numDets x numTimeSteps matrix of speeds (can be off course be the
%          same as Z in case a speedmap is requested). Also U may contain
%          NaN values.
%
%     X:   1 x numDets vector of increasing detector locations [m]. 
%
%     T:   1 x numTimeSteps vector of increasing time instants [s] 
%
% allowed param/value pairs (defaults in brackets)
%     'dx'     :  spatial resolution of resulting map [m] [default 100]. dx
%                 may also be a monotonically increasing vector of specific
%                 locations for which the filter must make an estimate
%     'dt'     :  temporal resolution of resulting map [s] (default 30). dt
%                 may also be a monotonically increasing vector of specific
%                 time instants for which the filter must make an estimate
%     'c_cong' :  characteristic speed for congestion [m/s] (default -18/3.6)
%     'c_free' :  characteristic speed for free-flow [m/s] (default 80/3.6)  
%                 (80% of desired speed V0 on empty road[25])
%     'vc'     :  parameter in weighting model [m/s] (default 60/3.6). Vc
%                 determines the boundary between congested and free-flow
%                 conditions
%     'dv'     :  parameter in weighting model [m/s] (default 10/3.6), dV
%                 determines the slope of the weighting function. A large
%                 value depicts more gradual change from free-flow to
%                 congested
%     'sig'    :  width of the smoothing kernel over space [m] (default 300);
%     'tau'    :  width of the smoothing kernel over time [s] (default 30);
%     'showprogress' : if set to a positive number then the procedure echos
%                 progress in percentages (with steps of size 'showprogress')
%                 %%[only for external display???!!!]
%     'xstart'  : start of resulting map [m], used in combination with dx,
%                 (default X(1))
%     'xend'    : end of resulting map [m], used in combination with dx,
%                 (default X(end))
%     'tstart'  : start of resulting map [s], used in combination with dt,
%                 (default T(1))
%     'tend'    : end of resulting map [s], used in combination with dt,
%                 (default T(end))
%     'Xmax'    : maximum width [m] of the filter - useful in
%                 case of very large datasets
%     'Tmax'    : maximum timeperiod [s] of the filter - useful in
%                 case of very large datasets
%     'phi'     : kernel, (default exponential), exponential or gauss
%     'direction': 'R' (default) ascending order, 'L': descending
%     'method'  : 'conventional' (default) method as implemented by Serge
%                 and Hans 
%                 OR 'crosscorrel' for cross-correlation method by Thomas
%                 OR 'fft' for fft method by Thomas
%
% output:
%     z:        n x m matrix of filtered traffic observations [in same units as Z]
%               where n = (X(end)-X(1)) / dx and  m = T(end)-T(1) / dt
%     u:        n x m matrix of filtered speeds [in km/h]
%               where n = (X(end)-X(1)) / dx and  m = T(end)-T(1) / dt
%     x:        1 x n vector of cell locations [m]
%     t:        1 x m vector of time instants [s]
% these are optional:
%     u:        n x m matrix of filtered traffic speeds 
%     w:        n x m matrix of filter weights              
%     zcong:    n x m matrix of z according to congested filter
%     zfree:    n x m matrix of z according to free-flow filter
%     ucong:    n x m matrix of u according to congested filter
%     ufree:    n x m matrix of u according to free-flow filter
%     Ncong:    n x m matrix of normalisation factor congested filter
%     Nfree:    n x m matrix of normalisation factor free-flow filter

% helper flags for updating progress
PROGRESSFLAG_NEW    = 0;
PROGRESSFLAG_UPDATE = 1;
PROGRESSFLAG_FINISH = 2;

% timing
try
    dummy = toc;
catch       % construction looks weird but it works
    tic;
end
starttime = toc;

% default parameters
dt       = 20;          % s
dx       = 50;          % m
c_cong   = -18/3.6;     % m/s ===>change to -16 as mean[Many literature: -18]
c_free   = 80/3.6;      % m/s
Vc       = 80/3.6;      % m/s
dV       = 10/3.6;      % m/s
sig      = 300;         % m
tau      = 30;          % s
showprogress = 1;      % percent
%alpha    = ones(size(Z));
phi      = @phi_gauss;
method   = 'conventional';
xtsource = 0;
Xmax = 1000;
Tmax = 360;

% other not so relevant defaults
%alpha       = ones(size(Z));
direction   = 'r';
% xstart, xend, tstart, tend are set later
doHeterogeneities = false;
debug = false;
FCDres = [50,10];

args = varargin;   %% varargin:set of orginal parameters[defined one by one][can be default(no input)!!!]
if length(args)>1
    for i = 1:2:length(args)
        switch lower(args{i})
            case 'dt'
                dt       = args{i+1};
            case 'dx'
                dx       = args{i+1};
            case 'c_cong'
                c_cong   = args{i+1};
            case 'c_free'
                c_free   = args{i+1};
            case 'vc'
                Vc       = args{i+1};
            case 'dv'
                dV       = args{i+1};
            case 'sig'
                sig      = args{i+1};
            case 'tau'
                tau      = args{i+1};
            case 'fcdres'
                FCDres   = args{i+1};
                if length(FCDres)~=2
                    FCDres = [50,10];
                end;
            case 'showprogress'
                showprogress = args{i+1}(1);
%             case 'xtsource'
%                 xtsource = args{i+1}(1);
            case 'xmax'
                Xmax = args{i+1};
            case 'tmax'
                Tmax = args{i+1};
            case 'alpha'
                alpha = args{i+1};
            case 'phi'
                if strcmpi(args{i+1}, 'gauss')
                    phi = @phi_gauss;
                end
            case 'direction'
                direction = lower(args{i+1});
            case 'xstart'
                xstart = args{i+1};
            case 'xend'
                xend = args{i+1};
            case 'tstart'
                tstart = args{i+1};
            case 'tend'
                tend = args{i+1};
            case 'heterogeneities'
                Heterogeneities = args{i+1};
                doHeterogeneities = false;
            case 'debug'
                debug = args{i+1};
                debugSaveFile = 'RL_ASM.MAT';
            case 'method'
                method = args{i+1};
        end;
    end;
end;

% progress monitoring
if ishandle(showprogress)
    h = showprogress;
    showprogress = 10;
else
    h =-1;
end;

doshow = showprogress>0;

% check important helper variables that determine start and end of the x-t
% grid for which we make estimates
if ~exist('xstart', 'var')
    xstart = X(1);
end
if ~exist('xend', 'var')
    xend = X(end);
end
if ~exist('tstart', 'var')
    tstart = T(1);
end
if ~exist('tend', 'var')
    tend = T(end);
end

% make sure Xmax fits with data
if numel(X)>=2
    dX = X(2:end) - X(1:end-1);
    Xmax = max([Xmax; dX(:)]);
end;

% set output variables for all methods
x = xstart:dx:xend;
t = tstart:dt:tend;

% construct sparse matrices [1:max(X),1:max(T)] in case the data is not
% provided in matrices but in tuples Z,X,T (and by implication U,X,T)
% TODO !!! klopt nog niet veel van
% if isvector(Z) 
%     Src.Z=Z;
%     Src.U=U;
%     Src.X=X;
%     Src.T=T;
%     % place FCD measurements on dx, dt grid
%     iX = ceil(X/FCDres(1));
%     iT = ceil(T/FCDres(2));
%     % construct matrices
%     Z = sparse(iX,iT,Z);
%     U = sparse(iX,iT,U);
%     % construct new rounded off space and time vectors
%     X = 0:FCDres(1):max(X);
%     T = 0:FCDres(2):max(T);
% end;
% ------------

if doshow, dowb(0,'filtering data: ',PROGRESSFLAG_NEW);end; 

% let's go!

    % conventional algorithm by Hans and Serge
    %disp('Xmax ='),disp(Xmax)
    %disp('Tmax ='),disp(Tmax)
    %disp('Vc ='),disp(Vc)
    
    % if Z and U are equal (if they both reflect the same matrix of
    % observations - "speeds") this speeds up calculations
    Ztst = Z(~isnan(Z));    %%Omit NaN; isnan = return not a number
    Utst = U(~isnan(U));
    dospd = (numel(Ztst) == numel(Utst)) && all( (Ztst(:)-Utst(:)) == 0);   
    
    % Now prepare the x and t variables
    % either we output the same resolution as source data ...
    
    % %%xtsource ~= 0, then the same resolution as X & T. Usually x&t are
    % %%determined by dx & dt
    if xtsource
        x = X;    %related to n   X: detectors   N   !!!row vector==>column
        t = T;    %related to m   T: timesteps   M   !!!row vector
        % ... or the resolution is determined by dx and dt
    else
        if length(dx) > 1
            x = dx;   
        end;
        if length(dt) > 1
            t = dt;
        end;
    end;
    
    % handy shortcuts for sizes of the grids
    n = length(x); m = length(t);
    N = length(X); M = length(T);
    
    % Then look for detectors which are structurally failing. If this is the
    % case in more than 95% of the time, the entire detector location / time
    % period is removed for efficiency reasons
    % ADD HVL 2015-06-02 no real need to do this ...
    gooddets = ones(size(X));  
    for i = 1:N
        if sum(isnan(Z(i,:))) >= 1.1 * M
            gooddets(i) = 0;
        end;
    end;
    idx = find(gooddets);
    Z = Z(idx,:); U = U(idx,:); X = X(idx);
    
    goodtimes = ones(size(T));
    for i = 1:M
        if sum(isnan(Z(:,i))) >= 1.1 * N
            goodtimes(i) = 0;
        end;
    end;
    idx = find(goodtimes);
    Z = Z(:,idx); U = U(:,idx); T = T(idx);
    
    % store positions of NaN values (these get zero weight in the filter)
    % Important: set the corresponding data to a non-NaN value - otherwise the
    % reconstruction doesn't work (NaN * zero is still NaN!)
    nans = isnan(Z) | isnan(U);     %%Second filter
    idxnan = find(nans);            % look for positions of original NaNs.
    Z(idxnan) = -1;
    U(idxnan) = -1;
    %alpha(idxnan) = 0;              %Similar to the function of previous Nidx{}
    
    % prepare for faster filtering with the fields Tmax and Xmax. If these are
    % set (>0) then "each" measurement Z(X,T) will affect only areas X-Xmax<x<X+Xmax
    % and T-Tmax<t<T+Tmax. This implies that we also have to index the NaN values
    % in these smaller space-time regions
    
    Xidx = cell(n,1);  % empty matrices  'NaN' elements
    for j = 1:n
        if Xmax ~= 0       
            % [Calculating n of z(x) values,several Z(X) values are needed for each z==]
            % In each cell of Xidx related to x(1*n), there are several Z values for a specific x cal.
            dist = abs(x(j)-X);  %% Generate a matrix with all subtract?!
            Xidx{j} = find(dist<=Xmax);  %%when 'dist <= Xmax' is true, the element is 1; otherwise it is NaN.
            %%==> return the "sequent order" of ture values
        else
            Xidx{j} = 1:N;
        end;
    end;
    Tidx = cell(m,1);  %In each cell,contains several Z(T) values for calculating z(t)[??????????z?,????Z?]
    for j = 1:m
        if Tmax ~= 0
            dist = abs(t(j)-T);
            Tidx{j} = find(dist<=Tmax);  %%!!Each small cell needs which measured Z values for smoothing
        else
            Tidx{j} = 1:M;
        end;
    end;
    
    
    % prepare grids of sizes "space x time" for faster indexing (otherwise we
    % would need n x m calls to "repmat", now just 2):
    XX = repmat(X(:), 1, M);   % X(:)  N*1 column vector   ==> N*M  corresponding to N*M U+Z
    TT = repmat(T(:)', N, 1);  % T(:)' 1*M row vector  ==> Improving calculation
    % xx = repmat(x(:), 1, m);
    % tt = repmat(t(:)', n, 1);
    
    
    % prepare output vars  [used for save output]
    zfree = zeros(n,m);
    zcong = zeros(n,m);
    ucong = zeros(n,m);
    ufree = zeros(n,m);
    Ncong = zeros(n,m);
    Nfree = zeros(n,m);
    z = zeros(n,m);
    u = zeros(n,m);    %% New variable
    w = zeros(n,m);
    
    % init counters [element by element!!!]
    K = n*m;
    k = 1;
    prevperc = 0; % for displaying progress
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % algorithm:
    % we run through each z(x,t) in the output space-time grid and calculated a
    % weighted sum of the filtered measurements Z in a region +/- {Xmax, Tmax}
    % around it. In case this region is not set ALL measurements Z are used for
    % each calculated value z  %%(Within the range of Xmax+Tmax, the related Z values are used for z)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:n % space
        % ADD 2009-10-23 ThS: Heterogeneities (Tristan Paper)
        % find heterogeneitie
        
        for j= 1:m % time
            
            % free flow filter (note: NaNs get zero weight)
            % %%Using Nidx{} to realize send NaN with 0 weight
            pfree = phi(XX(Xidx{i},Tidx{j})-x(i), TT(Xidx{i},Tidx{j})-t(j), ...
                sig, tau, c_free); %return a matric!!!
            
            % weighting by alpha
            wpfree = pfree;
            
            % normalizing factor
            Nfree(i,j) = sum(sum(wpfree,2));  %%First sum 2nd dimension,then sum 1st dimension!!!
            % the value according to the free flow then equals
            zfree(i,j) = sum(sum(wpfree .* Z(Xidx{i}, Tidx{j}), 2)) / Nfree(i,j);
            if ~dospd     %%If U~=Z
                ufree(i,j) = sum(sum(wpfree .* U(Xidx{i}, Tidx{j}), 2)) / Nfree(i,j);
            else
                ufree(i,j) = zfree(i,j);
            end;
            
            % congested filter (NaNs get zero weight)
            pcong = phi(XX(Xidx{i},Tidx{j})-x(i), TT(Xidx{i},Tidx{j})-t(j), ...
                sig, tau, c_cong);
            
            % weighting by alpha
            wpcong = pcong;
            
            % normalizing factor
            Ncong(i,j) = sum(sum(wpcong,2));
            % the value according to the congested flow filter then equals
            zcong(i,j) = sum(sum(wpcong .* Z(Xidx{i}, Tidx{j}), 2)) / Ncong(i,j);
            if ~dospd
                ucong(i,j) = sum(sum(wpcong .* U(Xidx{i}, Tidx{j}), 2)) / Ncong(i,j);
            else
                ucong(i,j) = zcong(i,j);
            end;
            
            % calculate weights depending on prevailing speeds (if the speed is high
            % the weight will be predominantly free-flow, and congested
            % otherwise)
            w(i,j) = 0.5 * (1 + tanh( (Vc - min(ufree(i,j),ucong(i,j))/3.6) / dV ));
            %  Be careful, if all the input speed values [km/h], then transform!!!
            % In textbook,u(x,t) is only used for calculated weighting factor
            % It denotes the minimum of the estimates of speed derived from two filters
            % u(i,j)=min(ufree(i,j),ucong(i,j))
            % Extension, output 'u' and/or (ucong+ufree)
            
            % display progress
            perc = 100*k/K;
            if (showprogress>0) && (perc-prevperc)>=showprogress
                dowb(perc,'filtering data: ',PROGRESSFLAG_UPDATE);
                prevperc = floor(perc);
            end;
            k=k+1;
            
        end;
    end;
    
    % This line calculates the final result:
    z = w .* zcong + (1-w) .* zfree;
    u = w .* ucong + (1-w) .* ufree;
 
if doshow, dowb(100,'filtering data: ',PROGRESSFLAG_UPDATE);end;

% runtime measurement
endtime = toc;
runtime = endtime - starttime;

% debug => save intermediate results
if debug
    save(debugSaveFile);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% auxilary function: the anisotropic filters (sub-function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linear exponential kernel (as used in THF paper equations)
function [p] = phi_exp(xx, tt, sig, tau, c)
    p = exp(- abs(xx)/sig - abs(tt - xx/c)/tau);
end

% quadratic exponential (=Gaussian) kernel (as used in the THF paper
% figures)
function [p] = phi_gauss(xx, tt, sig, tau, c)
    p = exp(- ((xx).^2)/(2*sig^2) - ((tt - xx/c).^2)/(2*tau^2));
end

% find the right zeta of kernel for detectector location x_det and filter
% output location x_filt
function modifier = getKernelModifier(Heterogeneities, x_det, x_filt, sigma)

    xr = Heterogeneities.location - x_filt;
    x = x_det - x_filt;
    zeta = Heterogeneities.zeta;
    alpha = ones(size(x));
    % adjust alpha, if there is a Heterogeneity xr
    if exist('xr', 'var')
        %wr = exp( -abs(xr)/sigma );
        for r = 1:numel(xr)
            alpha = alpha .* ( ...
                ((xr(r)<0) == (x<xr(r)))  .* (exp( -abs(xr(r)) / (zeta(r)*sigma) )) ...
                + ((xr(r)<0) == (x>=xr(r))) * 1 ...
                );
        end
    end
    modifier = alpha;

    % not so old, but still...
    %[w, modifier] = kernel1D(x_det-x_filt, sigma, Heterogeneities.location-x_filt);

    % % old %%%%%%%%%%%
    % find road section
    % sec_det  = find(Heterogeneities.location > x_det,  1, 'first');
    % sec_filt = find(Heterogeneities.location > x_filt, 1, 'first');
    %
    % % find involved heterogeneities (if sec_filt and sec_det in same section,
    % % then h_end < h_start and for loop is not executed)
    % h_start = min(sec_det, sec_filt);
    % h_end = max(sec_det, sec_filt) - 1;
    %
    % % calc. modifier by multiplying zetas over heterogeneities between x_det
    % % and x_filt
    % modifier = 1;
    % for h = h_start:h_end
    %     modifier = modifier * Heterogeneities.zeta(h);
    % end
    % %%%%%%%%%%%
end

function snapindex = snapToComb(comb, pointstosnap)
% snap points to an equidistant comb
% snapindex = snapToComb(comb, pointstosnap): comb is a vector of
% equidistant points, pointstosnap is a vector of points that are snapped
% to the comb, snapindex is a vector of same size of pointstosnap
% containing the index of closest point of the comb to the pointtosnap
%% begin
% assert that comb contains equidistant points
assert(all(abs(diff(diff(comb))) < 1e-10));     
% distance of comb points
dc = abs(comb(2)-comb(1));
[C,S] = meshgrid(comb,pointstosnap);
% ADD HvL 2017-06-17: Thomas did not account for points outside the comb,
% this is easily remedied. Note that the min operator implicitely chooses
% the first index it encounters (similar to the problem identified below)
diffSnap = abs(C-S);
[~,snapindex] = min(diffSnap,[],2);

% A = ( abs(C-S) <= dc/2 );            % find closest comb point
%     [row, col] = find(A);
%     % rare exception: pointosnap is precisely in the middle of two neighboring
%     % comb points, => snap to left neighbor
%     [~,m,~] = unique(row','first');
%     snapindex = col(m)';
% %% end
end

function dowb(pos,msg,updateflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% auxilary function: progress
% CALL
%     dowb(pos,msg,updateflag)
%
% INPUT
%     pos: 0..100 percentual progress
%     msg: string (e.g.: 'computing averages ...')
%     updateflag (defined in dittlab_loadConstants):
%       PROGRESSFLAG_NEW:    New message
%       PROGRESSFLAG_UPDATE: Only update the percentage (ignore message) 
%       PROGRESSFLAG_FINISH: Like update and add newline character
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper flags for updating progress
PROGRESSFLAG_NEW    = 0;
PROGRESSFLAG_UPDATE = 1;
PROGRESSFLAG_FINISH = 2;

% uniform handling: you can call this function with fractions. We
% don't check wether pos is between 0 and 100, that costs too much
% resources
if pos<=1
    pos=pos*100;
end;

switch updateflag
    case PROGRESSFLAG_NEW
        % new message
        fprintf('%s %3.0f%%\n',msg,pos);
    case PROGRESSFLAG_UPDATE
        % just update the percentage (message is assumed identical). Means 5 backspaces
        fprintf('\b\b\b\b\b%3.0f%%\n',pos);
    case PROGRESSFLAG_FINISH
        % Like update and now add newline
        fprintf('\b\b\b\b\b%3.0f%%\n',pos);
end

end
