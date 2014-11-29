function sim_index = simulationpath(T,CDF,start)
sim_index       = zeros(T,1);   % create a simulation path of length T
sim_index(1)    = start;        % seed starting point of simulation path
randpath        = rand(T);      % pick T random draws from uniform dist 
                                % between 0 and 1
for i = 2:T
    r = randpath(i);         
    if r <= CDF(sim_index(i-1),1)   
       sim_index(i) = 1;
                                % if r is smaller than the first value of
                                % the CDF, set tomorrow's index to 1.
    else
       diff            = CDF(sim_index(i-1),:) - r;
       sim_index(i)    = find(diff(1:end-1).*diff(2:end)<0) + 1;
                                % otherwise, find the first value where the
                                % CDF is greater than r by finding where
                                % the negative switches to the positive
                                % (add one because otherwise it gives the
                                % last value where the CDF is less than r)
%%%MS: the above works but is less intuitive.  Try:
%       {~,sim_index(i)] = min(diff(diff>0)); %find the minimum positive diff   
    %%%%AF: using min(diff(diff>0)) will always give sim_index=1 since
    %%%%diff(diff>0) creates a matrix of all the entries where diff>0, and
    %%%%the minimum is of course the first on the left.
     end
end

end