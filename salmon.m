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
% 2. stock growth       s_(t+1) = 
% 3. hatchery transformation
%                       f(inv) = 5*log(inv)
%                       This was an arbitrary decision, and should be
%                       changed once we have a good estimate.
% 4. harvest            harvest = S_t - escapement*S_t
%                       This should also get updated
% ------------------------------------------------------------------------
% Created: September 30, 2014
% Last Updated: October 6, 2014
% Amanda Faig
% ------------------------------------------------------------------------

clear all

% 1A. Set up the economic parameters
% ----------------------------------
p       = 5;            % price of harvest
c       = 0.075;        % cost of harvest
delta   = 1/1.03;       % discount factor
esc     = 0.4;          % excapement rule: for now just what portion 
                        % of the stock will be harvested
                        

% 1B. Set up the biological parameters
% ------------------------------------
R       = 11;          
alpha   = .1;          
K       = round((R-1)/alpha);          % carrying capacity

% 1C. Set up stochasticity
% ------------------------
z       = [0.3; 0; -0.3];       % the shocks to salmon growth
pz      = [0.3, 0.4, 0.3];      % the probability that each shock will 
                                % occur
                                
        % check that the probability vector makes sense
        % ---------------------------------------------
        if sum(pz) ~= 1                 

            disp('ERROR: improper probability matrix')
            break
        end

% 1C. Set up all other parameters
% -------------------------------
Svec    = linspace(0,K,K);      % vector of the state space
check   = 2;


% 2. Initial condition
% --------------------
Vold    = 2700*ones(1, K);  % sets initial guess of value function for 
                                % each stock leven in each period
                                               

% 3. Value Function Iteration
% ---------------------------

while check > .01           % keep going until max dev less than .01%
    for i = 1:length(Svec)  % loop over stocks
        S               = Svec(i);
        [Invtmp, Vtmp]  = ...
            fminbnd(@(inv)negpayoff(inv,esc,p,c,delta,R,alpha,Svec,S,Vold,z,pz),0,1000);
                            % temporary investment decision and value 
                            % function is equal to the levels that maximize 
                            % the payoff function (minimize -1*payoff)
        invstar(i)      = Invtmp;   
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
    
    colorvec    = [0.1 ,0.5 ,1 ; 1, 0, 0];  % creates colormap for plot
    
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

      
        