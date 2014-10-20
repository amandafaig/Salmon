function negvalue = negpayoff(inv,esc,p,c,delta,R,alpha,Svec,S,V,z,pz)
% Finds negative payoff from hatchery output 'hatch' when stock is S and
% the value function is V

% Harvest
% -------
harv = harvest(esc,S);

% Today's profit
% ---------------------------
today   = p*harv - c*harv^2 - inv;


% Stock growth formula:
% s_(t+1) = s_(t) + r*S_(t)*(1 - S_(t)/K) - harv*S + f(inv)
% ----------------------------------------------------------------
finv    = hatchery(inv);
SN      = max((1+z)*(R*S/(1+alpha*S) - harv + finv),1);

% The value of having s_(t+1) at the start of the next time period
% ----------------------------------------------------------------
Vnext = interp1(Svec,V,SN,'spline');  % a vector

negvalue = -( today + delta*pz*Vnext');
end

% Hatchery transformation formula
% -------------------------------
function fish = hatchery(inv)
         fish = 5*log(inv/100);
end

% Harvest function
function harv = harvest(esc,S)
         harv = S - esc*S;
end