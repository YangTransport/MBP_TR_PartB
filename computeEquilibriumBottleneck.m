function congestion = computeEquilibriumBottleneck(priorityQ,congestion,population)
%% Modifiy the discrete time 

Nstar = length(population.tstar); 
Nt = length(congestion.t);
t = congestion.t;
% tmin_add = min(population.tstar) - 1/congestion.S * sum(population.GroupSize.*priorityQ)/sum(population.GroupSize) ; % hours 
% tmax_add = max(population.tstar) + 1/congestion.S * sum(population.GroupSize.*priorityQ)/sum(population.GroupSize) ;
% dt_add  = 1/(sum(population.GroupSize) * congestion.S); % choose a time interval duration such that at most only one user pass each time

%% Calcualte the priority equlibrium 
SP_alpha = congestion.SP;
SP_alpha = reshape(SP_alpha,[Nt * Nstar,1]);
A = repmat(speye(Nt),1,Nstar); % AX<=capacity. ensures that every time is chosen by at most one person.
bPriority = ones(Nt,1) * congestion.PreemptableS/congestion.S;
Aeq = kron(speye(Nstar),ones(1,Nt));% Aeq*X=Na*ones(Nstar,1) (demand constraint) 
options = optimoptions('linprog','algorithm','interior-point','Display','off');
[arrivalRatesPriority,~,exitflag,~,lambdaPriority]=linprog(SP_alpha,A,bPriority, Aeq, population.GroupSize.* priorityQ,zeros(Nt * Nstar,1),[],options );
if exitflag ~= 1
   options = optimoptions('linprog','algorithm','dual-simplex','Display','off');
   [arrivalRatesPriority,~,exitflag,~,lambdaPriority]=linprog(SP_alpha,A,bPriority, Aeq, population.GroupSize.* priorityQ,zeros(Nt * Nstar,1),[],options );
end
congestion.queueingTimePriority = [];
congestion.queueingTimePriority = lambdaPriority.ineqlin;
congestion.CostPriority = -lambdaPriority.eqlin .* population.alpha ;
arrivalRatesPriority = reshape(arrivalRatesPriority,[Nt,Nstar]);
departureTimes = t - lambdaPriority.ineqlin;
% for i=1:Nstar
%     congestion.eqArrivalsPriority{i}=t(arrivalRatesPriority(:,i)>0);
%     congestion.eqDeparturesPriority{i}=departureTimes(arrivalRatesPriority(:,i)>0);
% end

%% Calculate the nonpriority equlibrium 

% bNonpriority = ones(Nt,1);
% bNonpriority(sum(arrivalRatesPriority,2)>0) =  1- congestion.PreemptableS/congestion.S;
bNonpriority = ones(Nt,1)- sum(arrivalRatesPriority,2);
bNonpriority(bNonpriority<0)=0;
% congestion.busyPeriod =t(sum(arrivalRatesPriority,2)>0) ;
% [arrivalRatesNonpriority,~,exitflag2,~,lambdaNonpriority]=linprog(SP_alpha,A,bNonpriority, Aeq, population.GroupSize.*(1 -priorityQ),zeros(Nt*Nstar,1),[],options);
% if exitflag2 ~= 1
   [arrivalRatesNonpriority,~,exitflag2,~,lambdaNonpriority]=linprog(SP_alpha,A,bNonpriority, Aeq, population.GroupSize.*(1 -priorityQ),zeros(Nt*Nstar,1),[],options);
% end
if exitflag2 ~= 1
    disp('exitflag2')
    disp(exitflag2)
end
congestion.queueingTimeNonpriority = lambdaNonpriority.ineqlin;
congestion.CostNonpriority = -lambdaNonpriority.eqlin .* population.alpha ;
arrivalRatesNonpriority = reshape(arrivalRatesNonpriority,[Nt,Nstar]);
departureTimesNonpriority = t-lambdaNonpriority.ineqlin;
for i = 1:Nstar
    congestion.eqArrivalsNonpriority{i}=t(arrivalRatesNonpriority(:,i)>0);
    congestion.eqDeparturesNonpriority{i}=departureTimesNonpriority(arrivalRatesNonpriority(:,i)>0);
end

end

