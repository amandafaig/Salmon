% First step towards programming economics of salmon hatchery/trucking
% problem: programming VFI problem with a single state variable (population
% size and a single choice variable (hatchery output), given a harvest
% rule, and stochasticity.
%
% Eventually we need to have population size, and some phenotypic (or 
% genotypic) variables as the sate variables, and hatchery output, trucking
% choices, and the harvest rule as choice variables.
%
% ------------------------------------------------------------------------
% There software includes the following files:
% ------------------------------------------------------------------------
% 1. "salmon_VFI.m"         Computes the VFI solution in the form of investment
%                       policy function.
% 2. "negpayoff.m"      Calculates the value function at each given state.
%                       This file also includes the hatchery transformation
%                       (from investment to fish) specification.
% ------------------------------------------------------------------------
% Created: September 30, 2014
% Last Updated: October 6, 2014
% Amanda Faig
% ------------------------------------------------------------------------

clear all
%dbstop in negpayoff

% 1A. Set up the economic parameters
% ----------------------------------
p       = 10;            % price per kg
c       = 0.075;        % cost of harvest
delta   = 1/1.03;       % discount factor
esc     = 0.4;          % excapement rule: for now just what portion 
                        % of the stock will be harvested
                        

% 1B. Set up the biological parameters
% ------------------------------------
R       = 201;          
alpha   = 2;          
K       = round((R-1)/alpha);          % carrying capacity

% 1C. Set up stochasticity
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


        % check that the probability vector makes sense
        % ---------------------------------------------
        if sum(pz) ~= 1               

            disp('ERROR: improper probability matrix')
            break
        end

% 1D. Set up all other parameters
% -------------------------------
Svec    = linspace(0,K,K+1);      % vector of the state space
check   = 2;


% 2. Initial condition
% --------------------
Vold    = linspace(8000,8300,K+1);  % sets initial guess of value function for 
                                    % each stock leven in each period
                                               

% 3. Value Function Iteration
% ---------------------------

while check > .01           % keep going until max dev less than .01%
    for i = 1:length(Svec)  % loop over stocks
        S               = i;
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

      
        