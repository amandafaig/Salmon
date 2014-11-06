function negvalue = negpayoff(esc,p,c,delta,R,alpha,Svec,S,V,z,pz)
% Finds negative payoff from hatchery output 'hatch' when stock is S and
% the value function is V

% Harvest
% -------
harv    = S - esc*S;

% Today's profit
% ---------------------------
today   = p*harv - c*harv^2;


% Stock growth formula:
% s_(t+1) = s_(t) + r*S_(t)*(1 - S_(t)/K) - harv 
% ----------------------------------------------------------------
SN      = (1+z)*(R*esc*S/(1+alpha*esc*S));

% The value of having s_(t+1) at the start of the next time period
% ----------------------------------------------------------------
Vnext = interp1(Svec,V,SN,'spline');  % a vector

negvalue = -(today + delta*pz*Vnext');
end


