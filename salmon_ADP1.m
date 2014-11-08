% STEP 3: Simple fishery harvest problem, post-decision ADP
%
% Run this code after running either "salmon_preVFI" (Step 1) or
% "salmon_postVFI" (Step 2) or both.  I recommend, unless you happen to
% have run Step 2 anyway, to run Step 1 as it converges quickly.  If you
% prefer the graphs do not superimpose, save and close the results from
% Step 1 and Step 2 before running Step 3.
%
% ------------------------------------------------------------------------
% There software includes the following files:
% ------------------------------------------------------------------------
% 1. "salmon_ADP.m"     Computes the ADP solution in the form of investment
%                       policy function.
% 2. "negpayoff.m"      Calculates the value function at each given state.
%                       This file also includes the hatchery transformation
%                       (from investment to fish) specification.
%
% ------------------------------------------------------------------------
% Created: Nov 6, 2014
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
numshock    = 10;                   % number of possible shock value
mean        = 1;                    % mean for distribution of shocks
rho         = 0.95;                 % persistence
sd          = .3;                   % SD for distribution of shocks 
z           = linspace(mean - 3*sd,mean + 3*sd, numshock); 
                                     
% Construct Markov transition matrix
% ----------------------------------
        a   = (1+rho)/2;
        b   = (1+rho)/2;
        pi  = [a; 1-a; 1-b; b];
        A   = [pi; zeros(length(pi),1); 1; zeros(length(pi)+1,1)];
        B   = [zeros(length(pi),1); pi; 1; zeros(length(pi)+1,1)];
        C   = [zeros(length(pi)+1,1); pi; zeros(length(pi),1); 1];
        D   = [zeros(length(pi)+1,1); zeros(length(pi),1); pi; 1];
        pi  = a*A + (1-a)*B + (1-q)*C + q*E;
        pi  = [pi(1); pi(2:    

% 1D. Specify state space
% -----------------------
Svec    = linspace(0,K,K+1);      

% 1E. Specify solution method parameters
% --------------------------------------
T       = 1000;         % number of periods per iteration (inner loop, t)
n       = 1;            % iteration count (outer loop, n)
check   = 2;
conv    = 0.01;         % convergence check: how close in percent the 
                        % estimates from two successive iterations must be 
                        % in order for the loop to stop
                                                
% 2. Initial conditions
% ---------------------
Vold    = 10^4*ones(length(Svec),length(z));    
                                % sets initial guess of value function for 
                                % each stock and each shock
Vold(1) = 0;                    % the value of zero fish is zero
S       = Svec(round(K/2));     % initial S value    

% Placeholders
% ------------
Vtmp    = Vold;
ESCstar = zeros(1,length(Svec));

while check > conv          % keep going until convergence criterion is met
                            % (max dev less than .01%)
    
    % 3. Set up sample shock
    % ----------------------
    
    zindex      = randsample(numshock,T);
    zsample(1)  = z(zindex);
    
   for i = 1:T
       z            = z(i); % This period's shock is from the randomly 
                            % generated shocks (3)
       [ESCtmp, V]  = ...
            fminbnd(@(esc)negpayoff(esc,p,c,delta,R,alpha,Svec,S,Vtmp,z,pz),0,1);
                            % temporary investment decision and value 
                            % function is equal to the levels that maximize 
                            % the payoff function (minimize -1*payoff)
       [~,SN]       = negpayoff(ESCtmp,p,c,delta,R,alpha,Svec,S,Vtmp,z,pz);                     
       S            = round(SN); % S(t+1)
       ESCstar(i)   = ESCtmp;   
                            % investment decisions for each level of stock
                            
       Vtmp(i)      = -V;   
                            % make current value function the starting 
                            % point for the next iteration
                            
  end
    dev         = abs((Vtmp - Vold)./Vold)*100; % calculate the maximum 
                                                % deviation (in percent) 
                                                % between iterations
    check       = max(dev); % check convergence criterion, 5A in algorithm 
    n           = n+1;      % add 1 to iteration counter, 5B in algorithm
    
    % -------------------------------
    % REGRESS HERE -- 5C in ALGORITHM
    % -------------------------------
    
    Vold        = Vtmp;     % set V(t=1, n+1) = V(t=T, n), 5D in algorithm
    
    
    % Plot the value function and policy function of each iteration
    % -------------------------------------------------------------
    
    colorvec    = [0.2 ,1 ,1 ; 1, 0, 0.5];  % creates colormap for plot
    
    if check > 0.01                         % so long as the loop will do 
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
        ylabel('Value Function value')
        hold on

        subplot(1,2,2)
        plot(Svec,invstar,'Color',color)
        xlabel('stock')
        ylabel('optimal investment')
        hold on
        pause(.1)
       
end

      
        