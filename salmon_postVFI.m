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
%dbstop in negpayoff

% 1A. Specify the economic parameters
% -----------------------------------
p       = 10;           % price per kg
c       = 0.075;        % cost of harvest
delta   = 1/1.03;       % discount factor
                        

% 1B. Specify the biological parameters
% -------------------------------------
R       = 201;          
alpha   = 2;          
K       = round((R-1)/alpha);   % carrying capacity

% 1C. Specify stochastic space
% ----------------------------
numz        = 10;       % number of possible shock value
mean        = 0;        % mean for distribution of shocks
sd          = .1;       % SD for distribution of shocks 
Zvec        = linspace(mean - 3*sd,mean + 3*sd, numz); 
                        % 100 points from 99.7% of the distribution
zbin        = Zvec - abs(Zvec(1)-Zvec(2))/2; 
                        % discretize normal pdf into bins 
zbin(numz+1)= Zvec(numz) + abs(Zvec(1)-Zvec(2))/2;
cdfz        = normcdf(zbin,mean,sd);
cdfz2       = [0, cdfz];
cdfz        = [cdfz, 1];
pz          = (cdfz - cdfz2);
pz          = reshape(pz,length(pz),1);

% 1D. Specify state space and action space
% ----------------------------------------
Svec                    = linspace(0,K,K+1);        % possible states     
Avec                    = linspace(0,1,50);         % possible actions
[S_sza, Z_sza, A_sza]   = meshgrid(Svec,Zvec,Avec); % possible combinations


% 1E. Specify solution method parameters
% --------------------------------------
check   = 2;
conv    = 0.01;         % convergence check: how close in percent the 
                        % estimates from two successive iterations must be 
                        % in order for the loop to stop
                        
% 2. Initial condition
% --------------------
Vold    = linspace(10^4,10^4+K,K+1);    % sets initial guess of value 
                                        % function for each stock leven in 
                                        % each period
Vold(1) = 0;                            % the value of no stock is 0

% 3. Value Function Iteration
% ---------------------------

% Today's Profit and Next Period's Stock for every (S,Z,A) combo
% --------------------------------------------------------------
pi_sza      = p.*(S_sza-A_sza.*S_sza) - c.*(S_sza-A_sza.*S_sza).^2;
                    % profit function
SN_sza      = (1.+Z_sza).*(R.*A_sza.*S_sza./(1.+alpha.*A_sza.*S_sza));
                    % stock 


% Placeholders
% ------------


while check > conv          % keep going until max dev less than .01%
    
      Vn_sza            = interp1(Svec,V,SN_sza,'spline');  
                            % Hypothetical value function next period for each (S,Z,A)
      V_sza             = pi_sza + delta*Vn_sza;            
                            % Value function today for each (S,Z,A)
      [V_sz, Ai_sz]     = max(V_sza,[],3);
                            % Value function 
                            
 
    
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
        plot(Svec,Vold,'Color',color)
        xlabel('stock')
        ylabel('Value Function')
        hold on

        subplot(1,2,2)
        plot(Svec,ESCstar,'Color',color)
        xlabel('stock')
        ylabel('optimal escapement (%)')
        hold on
        pause(.1)
        
end

      
        