% First step towards programming economics of salmon hatchery/trucking
% problem: programming VFI problem with a single state variable (population
% size and a single choice variable (hatchery output), given a harvest
% rule.
%
% Eventually we need to have population size, and some phenotypic (or 
% genotypic) variables as the sate variables, and hatchery output, trucking
% choices, and the harvest rule as choice variables.
%
% ------------------------------------------------------------------------
% There software includes the following files:
% ------------------------------------------------------------------------
% 1. "salmon.m"         Computes the VFI solution in the form of investment
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
% 2. stock growth       s_(t+1) = s_(t) + r*S_(t)*(1 - S_(t)/K) - harv*S + f(inv)
%                       where s is stock, r is intrinsic growth, K is
%                       carrying capacity, and f(inv) is the hatchery
%                       transformation formula
% 3. hatchery transformation
%                       f(inv) = 5*log(inv)
%                       This was an arbitrary decision, and should be
%                       changed once we have a good estimate.
% ------------------------------------------------------------------------
% September 30, 2014
% Amanda Faig
% ------------------------------------------------------------------------

clear all

% 1A. Set up the economic parameters
% ----------------------------------
p       = 1;            % price of harvest
c       = 0.075;        % cost of harvest
delta   = 1/1.1;        % discount factor
harv    = 0.5;          % harvest control rule: for now just what portion 
                        % of the stock will be harvested
                        

% 1B. Set up the biological parameters
% ------------------------------------
r       = 0.3;          % intrinsic growth rate
K       = 100;          % carrying capacity

% 1C. Set up all other parameters
% -------------------------------
Sgrid   = K;                    % determines how fine the grid is
Svec    = linspace(0,K,Sgrid);  % vector of the state space


% 2. Initial condition
% --------------------
V       = zeros(1, Sgrid);      % sets initial guess of value function for 
                                % each stock leven in each period
                                
% 3. Choose simulation length
% ---------------------------                                
T       = 40;                  


% 4. Value Function Iteration
% ---------------------------

for t = T:-1:1              % loop backwards over time
    for i = 1:length(Svec)  % loop over stocks
        S               = Svec(i);
        [Invtmp, Vtmp]  = ...
            fminbnd(@(inv)negpayoff(inv,harv,p,c,delta,r,K,Svec,S,V),0,1000);
                            % temporary investment decision and value 
                            % function is equal to the levels that maximize 
                            % the payoff function (minimize -1*payoff)
        invstar(i)      = Invtmp;   
                            % investment decisions for each level of stock
        V(i)            = -Vtmp;   
                            % make current value function the starting 
                            % point for the next iteration
    end
    
        %plot the value function and policy functions backwards through time
    subplot(1,2,1)
    plot(Svec,V)
    xlabel('stock')
    ylabel('Value Function value')
    hold on
    
    subplot(1,2,2)
    plot(Svec,invstar)
    xlabel('stock')
    ylabel('optimal investment')
    hold on
    pause(.1)
end

      
        