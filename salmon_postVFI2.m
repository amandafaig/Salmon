% STEP 3: Simple fishery harvest problem, post-decision state VFI.
%
% Run this code after running "salmon_preVFI" (Step 1).  The graphs will
% superimpose, allowing you to see how the two solve differently.  If you
% prefer they don't superimpose, save and close the results from Step 1
% before running Step 2.
%
% ------------------------------------------------------------------------
% There software calls on the following files:
% ------------------------------------------------------------------------
% 1. "salmon_postVFI.m"    Computes the VFI solution in the form of 
%                           investment policy function.
% 2. "negpayoff.m"          Calculates the value function at each given state.
%                           This file also includes the hatchery transformation
%                           (from investment to fish) specification.
% ------------------------------------------------------------------------
% Created: Nov 4, 2014
% Last updated: Nov 8, 2014
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
zv      = linspace(-numz,numz,2*numz+1);   

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
[z0_zs,shock_zs]    = ndgrid(zv,shockvec);
                        % a grid of all possible z's one might be at today
                        % (z0_sz) and all possible shocks that might hit 
                        % during the transition (shock_sz)
                        % _zs is to remind that each row represents a start
                        % date and each column represents a possible shock
[~,pshock_zs]       = ndgrid(zv,pshock);                    
z1_zs               = round(z0_zs*(1-sigma) + shock_zs*sigma);
                        % what tomorrow's z will be (z1) given today's (z0)
                        % and the shock
zi1_zs              = z1_zs + numz + 1;
                        % index value of tomorrow's z given today's z.
pz1_z0z1            = zeros(length(zv));

for t=1:length(zv)
    for j=1:length(zv)
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
fzavg   = @(beta) (1.1 + beta*abs(zv))*freqvisited - 1;  
beta    = fsolve(fzavg,-0.01);
fz      = 1.1 + beta*abs(zv); 	% Describes the recruitment shock that comes 
                            % from being so many days away from z=0.

% 1D. Specify state space and action space
% ----------------------------------------
Sv    = linspace(0,K,K+1);          % possible endogenous states 
Si    = Sv + 1;                     % note that the index is just 1+value        
Av    = linspace(0,1,100);          % possible actions 
                                    % (% escapement is between 0 and 1)
[S_sza,Fz_sza,A_sza]   = ndgrid(Sv,fz,Av);    % possible combinations     
[S_sz,Z_sz]            = ndgrid(Sv,zv);
% 1E. Specify solution method parameters
% --------------------------------------
conv    = 0.01; %(1-delta)*100;      % convergence check: how close in percent the 
                        % estimates from two successive iterations must be 
                        % in order for the loop to stop
check   = 2*conv;       % ensure that the loop will start

                      
% 2. Initial conditions
% ---------------------
V        = 1.1*10^4*ones(length(Sv),length(zv));    

                                % sets initial guess of value function for 
                                % each stock and each shock
V(1,:)   = 0;                   % the value of zero fish is zero
% initialize matrices for loop
% ----------------------------
Vnew     = 0*V;                 
Action   = zeros(length(Sv),length(zv));
Vn_sza   = 0*S_sza;
n=0;
% 3. Value Function Iteration
% ---------------------------

% Today's Profit and Next Period's Stock for every (S,Z,A) combo
% --------------------------------------------------------------
pi_sza      = pi(S_sza,Fz_sza,A_sza);   % profit for every s_t, z_t-1, a_t               
SN_sza      = SN(S_sza,Fz_sza,A_sza);   % s_t+1 for every s_t, z_t-1, a_t

while check > conv          % keep going until max dev less than .01%
        for k=1:length(Av)
            Vn_sza(:,:,k) = interp2(Z_sz,S_sz,V,Z_sz,SN_sza(:,:,k),'spline');  
                            % Hypothetical value function next period for 
                            % each (S,Z,A)
        end
        V_sza           = pi_sza + delta*Vn_sza;            
                            % Value function today for each (S,Z,A)
        [V_sz, Ai_sz]   = max(V_sza,[],3);
                            % Value function 
        A_sz            = Av(Ai_sz);

        
    for z = 1:length(zv)
        Vnew(:,z)   = V_sz*pz1_z0z1(z,:)';
        Action(:,z) = A_sz*pz1_z0z1(z,:)';
    end


        dev             = abs((Vnew - V)./V)*100;   % calculate the maximum
                                                    % deviation (in percent) 
                                                    % between iterations
        check           = max(max(dev));      
        V               = Vnew;
        n   = n+1;
        dbstop if n==1000;
    
    % Plot the value function and policy function of each iteration
    % -------------------------------------------------------------
    
        colorvec    = [0.1 ,0.5 ,1 ; 1, 0, 0];  % creates colormap for plot
    
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
            entry = num2str(zv(i));
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

      
        