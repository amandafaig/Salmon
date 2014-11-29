% STEP 4: Simple fishery harvest problem, post-decision ADP
%
% Run this code after running either "salmon_preVFI" (Step 1),
% "salmon_postVFI" (Step 2), or "salmon_postVFI2" (Step 3).  
%
%
% ------------------------------------------------------------------------
% Created: Nov 6, 2014
% Last updated: Nov 20, 2014
% Amanda Faig
% ------------------------------------------------------------------------

clear all

% 1A. Specify the economic parameters
% -----------------------------------
p       = 10;           % price per kg
c       = 0.075;        % cost of harvest
delta   = 1/1.03;       % discount factor

% Specify the profit function
% ---------------------------
pi      = @(S,Fz,A) p.*(Fz.*S.*(1-A)) - c.*(Fz.*S.*(1-A)).^2;

% 1B. Specify the biological parameters
% -------------------------------------
R       = 201;     % avg offspring per female 
K       = 100;     % carrying capacity
alpha   = (R-1)/K;       

% Specify the stock transition function
% -------------------------------------
SN      = @(S,Fz,A) round((R.*A.*Fz.*S)./(1+alpha*A.*Fz.*S));
 
% 1C. Specify stochastic space
% ----------------------------
% z = 0 implies that the fish arrived to the ocean on the perfect day;
% production gets a boom from that most fortunate event.  The effect on
% production decreases as the arrival date deviates further from 1,
% becoming negative as arrival gets far enough from z=0.

numz    = 5;            % maximum distance from 0 you can go
Zv      = linspace(-numz,numz,2*numz+1);   

% Discretize shock distribution and calculate absolute probability of any 
% given bin being chosen.
% ----------------------------------------------------------------------
numshock    = 10;           % number of possible shock value
meanshock   = 0;            % mean for distribution of shocks
sdshock     = 2;            % SD for distribution of shocks 
shockvec    = linspace(meanshock - 3*sdshock,meanshock + 3*sdshock, numshock); 
shockbin    = shockvec - abs(shockvec(1)-shockvec(2))/2; 
                            % discretize normal pdf into bins 
shockbin(numshock+1)...
            = shockvec(numshock) + abs(shockvec(1)-shockvec(2))/2;
cdfshock    = normcdf(shockbin,meanshock,sdshock);
cdfshock2   = [0, cdfshock];
cdfshock    = [cdfshock, 1];
pshock      = (cdfshock - cdfshock2);
pshock      = pshock(2:numshock+1);
pshock      = pshock./sum(pshock);

% Construct Markov transition matrix
% ----------------------------------
% The date of arrival is naturally tending towards 0 at the rate (1-sigma);
% random shocks may push the arrival date further or closer to 0.
% z_t+1 = z_t*(1-sigma) + shock(t)*sigma

sigma               = 0.4;   
                        % sigma must be between 0 & 1. 
                        % 1 implies shocks are all just random draws.
                        % 0 implies no shocks or regression.
[z0_zs,shock_zs]    = ndgrid(Zv,shockvec);
                        % a grid of all possible z's one might be at today
                        % (z0_sz) and all possible shocks that might hit 
                        % during the transition (shock_sz)
                        % _zs is to remind that each row represents a start
                        % date and each column represents a possible shock
[~,pshock_zs]       = ndgrid(Zv,pshock);                    
z1_zs               = round(z0_zs*(1-sigma) + shock_zs*sigma);
                        % what tomorrow's z will be (z1) given today's (z0)
                        % and the shock
zi1_zs              = z1_zs + numz + 1;
                        % index value of tomorrow's z given today's z.
pz1_z0z1            = zeros(length(Zv));

for t=1:length(Zv)
    for j=1:length(Zv)
        pz1_z0z1(t,j) = pz1_z0z1(t,j) + sum(pshock(zi1_zs(t,:)==j));
    end
end

CDFz1_z0z1 = cumsum(pz1_z0z1,2);
                

% Construct bonus/penalty matrix
% ------------------------------
% to approximate likelihood of arriving at any given day, run simulated
% path and count 
[~,eigvalues,lefteig]   = eig(pz1_z0z1);
[~,eig1]                = min(abs(sum(eigvalues-1)));
freqvisited             = lefteig(:,eig1);
freqvisited             = freqvisited/sum(freqvisited);

% Let the recruitment boom f(z) when z=0 be a 10% increase. 
% Assume linearity for now.
% f(z) = f(z=0) + beta*abs(z)
% --------------------------------------------------------
fzavg   = @(beta) (1.1 + beta*abs(Zv))*freqvisited - 1;  
beta    = fsolve(fzavg,-0.01);
fz      = 1.1 + beta*abs(Zv); 	% Describes the recruitment shock that comes 
                            % from being so many days away from z=0.

% 1D. Specify state space and action space
% ----------------------------------------
Sv    = linspace(0,K,K+1);          % possible endogenous states 
Si    = Sv + 1;                     % note that the index is just 1+value        
Av    = linspace(0,1,100);          % possible actions 
                                    % (% escapement is between 0 and 1)
[S_sza,Fz_sza,A_sza]   = ndgrid(Sv,fz,Av);    % possible combinations  
[S_sz, Z_sz]           = ndgrid(Sv,Zv);

% 1E. Specify solution method parameters
% --------------------------------------
T       = 10000;        % number of periods per iteration (inner loop, t)
n       = 1;            % iteration count (outer loop, n)
conv    = .01;          % convergence check: how close in percent the 
                        % estimates from two successive iterations must be 
                        % in order for the loop to stop
check   = 2*conv;       % ensure that the loop will start
stepsize= 0.05;          % stepsize parameter
% stepsize = @(
                      
% 2. Initial conditions
% ---------------------
Vbar        = 1.1*10^4*ones(length(Sv),length(Zv));    
                                    % sets initial guess of value function for 
                                    % each stock and each shock
Vbar(1,:)   = 0;                    % the value of zero fish is zero
V           = Vbar;                 % save a copy of initial V to use as
                                    % comparison in convergence check
Action      = 0.5*ones(length(Sv),length(Zv));
                                    % mostly a placeholder for the action
                                    % (aka policy) function
% Set up regression information
% -----------------------------
valuefun    = @(b,x) b(1) + b(2)*x(:,1).^b(3) + b(4)*abs(x(:,2)).^b(5);
                                    % the functional form I believe the
                                    % value function will take, so we can
                                    % run regressions to use the data we
                                    % collect each iteration.
                                    % x(:,1) represents the states
                                    % x(:,2) represents the arrival dates
beta0       = [10501,157.45,.275,-3.7,1.55];   % initial guess for beta
weightfun   = @(numvisit) (numvisit+1)/sum(sum(numvisit+1));
                                    % weights for the regression
numvisit    = 0*S_sz;               % placeholder for numvisit counter
% -------
% Iterate
% -------

% Today's Profit and Next Period's Stock for every (S,Z,A) combo
% --------------------------------------------------------------
pi_sza      = pi(S_sza,Fz_sza,A_sza);   % profit for every s_t, z_t-1, a_t               
SN_sza      = SN(S_sza,Fz_sza,A_sza);   % s_t+1 for every s_t, z_t-1, a_t
SNi_sza      = 0*SN_sza;
%%%% Is there a non-loop way to do this?
for t=1:size(SN_sza,1)
    for j=1:size(SN_sza,2)
        for k=1:size(SN_sza,3)
            SNi_sza(t,j,k) = Si(Sv==SN_sza(t,j,k));
        end
    end
end

while check > conv          % keep going until convergence criterion is met
                            % (max dev less than .01%)
    
    si              = Si(max(round(rand()*length(Si)),1));  % random starting stock
    
    % 3. Set up shock path for iteration
    % ----------------------------------
    zi_path         = zeros(T,1);   % create a simulation path of length T
    zi_path(1)      = max(1,round(rand()*length(Zv)));       
                                    % randomly choose starting point of 
                                    % simulation path
    randpath        = rand(T);      % pick T random draws from uniform dist 
                                    % between 0 and 1
    for t = 2:T
        r = randpath(t);         
        if r <= CDFz1_z0z1(zi_path(t-1),1)   
            zi_path(t) = 1;
                                % if r is smaller than the first value of
                                % the CDF, set tomorrow's index to 1.
        else
            diff          = CDFz1_z0z1(zi_path(t-1),:) - r;
            zi_path(t)    = find(diff(1:end-1).*diff(2:end)<0) + 1;
                                % otherwise, find the first value where the
                                % CDF is greater than r by finding where
                                % the negative switches to the positive
                                % (add one because otherwise it gives the
                                % last value where the CDF is less than r)
        end
    end

    % Iterate over T periods
    % ----------------------
    for t = 2:T
        zlasti          = zi_path(t-1); 
        zi              = zi_path(t);
        Si_a            = squeeze(SNi_sza(si,zlasti,:));  
                            % given s_t and z_t-1, 
                            % s_t+1 is a function of a_t (following index)     
        Vn_s            = V(:,zi);
                            % given z_t, 
                            % the possible values for V(s_t+1,z_t)
        Vn_a            = Vn_s(Si_a);  
                            % given s_t+1(a_t), the values for V(a_t,z_t) 
        V_a             = reshape(pi_sza(si,zi,:),100,1) + delta*Vn_a;            
        [Vtilde, Ai]    = max(V_a);     
                            % 4A 
        Vbar(si,zlasti) = (1-stepsize)*Vbar(si,zlasti) + stepsize*Vtilde;  
                            % 4B 
        Action(si,zi)   = Ai;       
                            % update action function
        numvisit(si,zlasti)...
                        = numvisit(si,zlasti)+1;
        Vbar(si,length(Zv)+1-zlasti)...    
                        = Vbar(si,zlasti); 
                            % exploit the symmetry
        Action(si,length(Zv)+1-zi)      = Ai;
        si              = max(2,Si_a(Ai));   
                            % When si=1 the system crashes.
    end
  
    dev         = abs((Vbar - V)./V)*100; % calculate the maximum 
                                          % deviation (in percent) 
                                          % between iterations
    check       = max(max(dev));    % check convergence criterion, 5A 
    n           = n+1;              % add 1 to iteration counter, 5B 
    
    % Run regressions to use information collected this iteration more
    % efficiently (5C)
    % ----------------------------------------------------------------
    % First we must eliminate all entries were S=0.  If we don't, this will
    % result in a regression that just won't run (since the function is
    % discontinuous at S=0.  
    
    S_sz0       = S_sz(2:end,:);
    Z_sz0       = Z_sz(2:end,:);
    Vbar_0      = Vbar(2:end,:);
    numvisit0   = numvisit(2:end,:);
    x(:,1)      = reshape(S_sz0,numel(S_sz0),1);  % states
    x(:,2)      = reshape(Z_sz0,numel(Z_sz0),1);  % shocks
    y           = reshape(Vbar_0,numel(Vbar_0),1);  % value
    weights     = weightfun(numvisit0);
    w           = reshape(weights,numel(weights),1);
    mdl         = fitnlm(x,y,valuefun,beta0,'Weights',w); 
                                                % regress
    Vpred       = predict(mdl,x);               % predict Vbar(s,z)
    Vbar_n0     = reshape(Vpred,length(Sv),length(Zv));
    Vbar        = [Vbar(1,:); Vbar_n0];         % add the rows
    beta0       = mdl{:,{'Estimate'}};          % use the betas from this 
                                                % iteration as next 
                                                % iterations' initial guess
    V           = Vbar;     % set V(t=1, n+1) = V(t=T, n), 5D 
    
    
    % Plot the value function and policy function of each iteration
    % -------------------------------------------------------------
    
    colorvec    = [0.2 ,1 ,1 ; 1, 0, 0.5];  % creates colormap for plot
    
    if check > conv                         % so long as the loop will do 
                                            % another iteration, plot using
                                            % the first color of the
                                            % colormap
        color = colorvec(1,:);
        
    else
        color = colorvec(2,:);              % otherwise, plot with second 
                                            % color
    end
        Vplot  = V*freqvisited;
        Aplot  = Action*freqvisited;
        
        for i = 1:6
            entry = num2str(Zv(i));
            subplot(2,7,i)
            plot(Sv,V(:,i),'Color',color)
            xlabel('stock')
            ylabel(['Value Function, abs(z)=',entry])
            hold on

            subplot(2,7,i+7)
            plot(Sv,Action(:,i),'Color',color)
            xlabel('stock')
            ylabel(['optimal escapement (%), abs(z)=',entry])
            hold on
        end       
        subplot(2,7,7)
        plot(Sv,Vplot,'Color',color)
        xlabel('stock')
        ylabel('Value Function, AVERAGE')
        hold on

        subplot(2,7,14)
        plot(Sv,Aplot,'Color',color)
        xlabel('stock')
        ylabel('optimal escapement (%), AVERAGE')
        hold on  
        pause(.1)
end

      
        