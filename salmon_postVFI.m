% STEP 2: Simple fishery harvest problem, post-decision state VFI.
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
R       = 201;          
K       = 100;      
alpha   = (R-1)/K;   

% Specify the stock transition function
% -------------------------------------
SN      = @(S,Fz,A) (R.*A.*Fz.*S)./(1+alpha*A.*Fz.*S);
 
% 1C. Specify stochastic space
% ----------------------------
numz        = 10;       % number of possible shock value
mean        = 1;        % mean for distribution of shocks
sd          = .1;       % SD for distribution of shocks 
Zvec        = linspace(mean - sd,mean + sd, numz); 
                        % 100 points from 99.7% of the distribution
zbin        = Zvec - abs(Zvec(1)-Zvec(2))/2; 
                        % discretize normal pdf into bins 
zbin(numz+1)= Zvec(numz) + abs(Zvec(1)-Zvec(2))/2;
cdfz        = normcdf(zbin,mean,sd);
cdfz2       = [0, cdfz];
cdfz        = [cdfz, 1];
pz          = (cdfz - cdfz2);
pz          = pz(2:numz+1);
pz          = pz./sum(pz);
pz          = reshape(pz,length(pz),1);

% 1D. Specify state space and action space
% ----------------------------------------
Svec                    = linspace(0,K,K+1);        % possible states     
Avec                    = linspace(0,1,100);        % possible actions
[S_sza, Z_sza, A_sza]   = ndgrid(Svec,Zvec,Avec);   % possible combinations


% 1E. Specify solution method parameters
% --------------------------------------
conv    = 0.01;         % convergence check: how close in percent the 
                        % estimates from two successive iterations must be 
                        % in order for the loop to stop
check   = 2*conv;       % ensure that the loop will start
                        
% 2. Initial condition
% --------------------
V    = linspace(10^4,10^4,K+1);    % sets initial guess of value 
                                        % function for each stock leven in 
                                        % each period
V(1) = 0;                            % the value of no stock is 0
V    = reshape(V,length(V),1);

% 3. Value Function Iteration
% ---------------------------

% Today's Profit and Next Period's Stock for every (S,Z,A) combo
% --------------------------------------------------------------
pi_sza      = pi(S_sza,Z_sza,A_sza);   % profit for every s_t, z_t-1, a_t               
SN_sza      = SN(S_sza,Z_sza,A_sza);   % s_t+1 for every s_t, z_t-1, a_t

while check > conv          % keep going until max dev less than .01%
    
        Vn_sza          = interp1(Svec,V,SN_sza,'spline');  
                            % Hypothetical value function next period for 
                            % each (S,Z,A)
        V_sza           = pi_sza + delta*Vn_sza;            
                            % Value function today for each (S,Z,A)
        [V_sz, Ai_sz]   = max(V_sza,[],3);
                            % Value function 
        Vnew            = V_sz*pz;                      
 
        dev             = abs((Vnew - V)./V)*100;   % calculate the maximum 
                                                    % deviation (in percent) 
                                                    % between iterations
        check           = max(dev);                 
        V               = Vnew;
        A_sz            = Avec(Ai_sz);     
        A               = A_sz*pz;
    
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
         
        subplot(1,2,1)
        plot(Svec,V,'Color',color)
        xlabel('stock')
        ylabel('Value Function')
        hold on

        subplot(1,2,2)
        plot(Svec,A,'Color',color)
        xlabel('stock')
        ylabel('optimal escapement (%)')
        hold on
        pause(.1)
        
end

      
        