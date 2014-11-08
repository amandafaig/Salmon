% STEP 1: Simple fishery harvest problem, pre-decision state VFI
%
% Run this code before running "salmon_postVFI" (Step 2).  The graphs will
% superimpose, allowing you to see how the two solve differently.  If you
% prefer they don't superimpose, save and close the results from Step 1
% before running Step 2.
%
% ------------------------------------------------------------------------
% There software calls on the following files:
% ------------------------------------------------------------------------
% 1. "salmon_preVFI.m"         Computes the VFI solution in the form of investment
%                       policy function.
% 2. "negpayoff.m"      Calculates the value function at each given state.
%                       This file also includes the hatchery transformation
%                       (from investment to fish) specification.
% ------------------------------------------------------------------------
% Created: September 30, 2014
% Last Updated: Nov 8, 2014
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
K       = round((R-1)/alpha);       % carrying capacity

% 1C. Specify stochastic space
% ----------------------------
numz        = 10;                   % number of possible shock value
mean        = 0;                    % mean for distribution of shocks
sd          = .1;                   % SD for distribution of shocks 
z           = linspace(mean - 3*sd,mean + 3*sd, numz); 
                                    % 100 points from 99.7% of the 
                                    % distribution
zbin        = z - abs(z(1)-z(2))/2; % discretize normal pdf into bins 
                                    % around shock values
zbin(numz+1)= z(numz) + abs(z(1)-z(2))/2;
cdfz        = normcdf(zbin,mean,sd);
cdfz2       = [0, cdfz];
cdfz        = [cdfz, 1];
pz          = (cdfz - cdfz2);


% 1D. Specify state space
% -----------------------
Svec    = linspace(0,K,K+1);      

% 1E. Specify solution method parameters
% --------------------------------------
check   = 2;
conv    = 0.01;         % convergence check: how close in percent the 
                        % estimates from two successive iterations must be 
                        % in order for the loop to stop
% 2. Initial condition
% --------------------
Vold    = linspace(10^4,10^4+K,K+1);  % sets initial guess of value 
                                        % function for each stock leven in 
                                        % each period
Vold(1) = 0;                            % the value of no stock is 0

% 2. Initial condition
% --------------------
Vold    = linspace(10^4,10^4+K,K+1);  % sets initial guess of value 
                                        % function for each stock leven in 
                                        % each period
Vold(1) = 0;                            % the value of no stock is 0
                                               
% 3. Value Function Iteration
% ---------------------------

% Placeholders
% ------------
V       = zeros(1,length(Svec));
ESCstar = V;

while check > conv          % keep going until max dev less than .01%
    for i = 1:length(Svec)  % loop over stocks
        S               = Svec(i);
        [ESCtmp, Vtmp]  = ...
            fminbnd(@(esc)negpayoff(esc,p,c,delta,R,alpha,Svec,S,Vold,z,pz),0,1);
                            % temporary investment decision and value 
                            % function is equal to the levels that maximize 
                            % the payoff function (minimize -1*payoff)
        ESCstar(i)      = ESCtmp;   
                            % investment decisions for each level of stock
        V(i)            = -Vtmp;   
                            % make current value function the starting 
                            % point for the next iteration
    end
    
    dev         = abs((V - Vold)./V)*100;   % calculate the maximum 
                                            % deviation (in percent) 
                                            % between iterations
    check       = max(dev);                 
    Vold        = V;
    
    % Plot the value function and policy function of each iteration
    % -------------------------------------------------------------
    
    colorvec    = [0.5 ,1 ,0.8 ; 1, 0, 1];  % creates colormap for plot
    
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
        ylabel('optimal escapement (fraction)')
        hold on
        pause(.1)
        
end

      
        