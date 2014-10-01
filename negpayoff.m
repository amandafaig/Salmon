function negvalue = negpayoff(inv,harv,p,c,delta,r,K,Svec,S,V)
% Finds negative payoff from hatchery output 'hatch' when stock is S and
% the value function is V

% Today's profit
% ---------------------------
today   = p*harv*S - c*(harv*S)^2 - inv;


% Stock growth formula:
% s_(t+1) = s_(t) + r*S_(t)*(1 - S_(t)/K) - harv*S + f(inv)
% ----------------------------------------------------------------
finv    = hatchery(inv);
SN      = S + r*S*(1 - S/K) - harv*S + finv;


% The value of having s_(t+1) at the start of the next time period
% ----------------------------------------------------------------
Vnext = interp1(Svec,V,SN,'pchip'); 

negvalue = -( today + delta*Vnext);
end

% Hatchery transformation formula
% -------------------------------
function fish = hatchery(inv)
            fish = 5*log(inv);
end