% Second step towards programming economics of salmon hatchery/trucking
% problem: programming ADP problem with a single state variable (population
% size and a single choice variable (hatchery output), given a harvest
% rule, and stocasticity.
%
% Eventually we need to have population size, and some phenotypic (or 
% genotypic) variables as the sate variables, and hatchery output, trucking
% choices, and the harvest rule as choice variables.
%
% ------------------------------------------------------------------------
% There software includes the following files:
% ------------------------------------------------------------------------
% 1. "salmon_ADP.m"         Computes the ADP solution in the form of investment
%                       policy function.
% 2. "negpayoff.m"      Calculates the value function at each given state.
%                       This file also includes the hatchery transformation
%                       (from investment to fish) specification.
%
% ------------------------------------------------------------------------
% The formulas I use are:
% ------------------------------------------------------------------------
% 1. current profit     pi = p*harv*S - c*(harv*S)^2 - inv
%                       today's profits are equal to the price per pound of
%                       fish, minus the cost of fishing each pound, minus
%                       the amount of money invested into the hatchery.
% 2. stock growth       s_(t+1) = 
% 3. hatchery transformation
%                       f(inv) = 5*log(inv)
%                       This was an arbitrary decision, and should be
%                       changed once we have a good estimate.
% 4. harvest            harvest = S_t - escapement*S_t
%                       This should also get updated
% ------------------------------------------------------------------------
% Created: October 15, 2014
% Last updated: October 20, 2014
% Amanda Faig
% ------------------------------------------------------------------------

clear all
%dbstop in negpayoff
rng(2061);  % Seed the random number generator

% 1A. Determine simulation length
% -------------------------------
T       = 1000;

% 1B. Set up the economic parameters
% ----------------------------------
p       = 10;           % price per kg
c       = 0.075;        % cost of harvest
delta   = 1/1.03;       % discount factor
esc     = 0.4;          % excapement rule: for now just what portion 
                        % of the stock will be harvested
                        
% 1C. Set up the biological parameters
% ------------------------------------
R       = 201;          
alpha   = 2;          
K       = round((R-1)/alpha);          % carrying capacity

% 1D. Set up stochasticity
% ------------------------
mean        = 0;                    % mean for distribution of shocks
sd          = .1;                    % SD for distribution of shocks 
z           = linspace(mean - 3*sd,mean + 3*sd, 100); 
                                % 100 points from 99.7% of the distribution
cz          = normcdf(z,mean,sd);
cz2         = [0, cz];
cz2(end)    = [];
pz          = (cz - cz2);
pz(end)     = 0;
pz(end)     = 1 - sum(pz);

% 1E. Set up all other parameters
% -------------------------------
Svec    = linspace(0,K,K+1);      % vector of the state space
check   = 2;

% 2. Initial conditions
% ---------------------
Vold    = linspace(8000,8300,K+1);  % sets initial guess of value function for 
                                    % each stock leven in each period
V               = Vold;                
j               = round(K/2);
S               = Svec(j);
count           = 0;
invstar         = 500*ones(K+1,1);

% 3. Value Function Iteration
% ---------------------------

while check > .01           % keep going until max dev less than .01%
   for i = 1:T
    [Invtmp, Vtmp]  = ...
            fminbnd(@(inv)negpayoff(inv,esc,p,c,delta,R,alpha,Svec,S,Vold,z,pz),0,1000);
    V(j)            = -Vtmp;
    invstar(j)      = Invtmp;
    j               = max(round(rand()*(K+1)),1);
    S               = Svec(j);  
  end
    dev         = abs((V - Vold)./V)*100;   % calculate the maximum 
                                            % deviation (in percent) 
                                            % between iterations
    check       = max(dev);                 
    Vold        = V;
     
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

      
        