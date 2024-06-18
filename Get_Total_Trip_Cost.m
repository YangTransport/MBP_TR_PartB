%% This is a function for calculating the resulting total trip cost given any q

function [TotalTripCost, Trip_Costs_Nonpriority, Trip_Costs_Priority] = Get_Total_Trip_Cost(priorityQ, PreemptRatio, population, congestion )
disp(['Computer equlibirum cost under q: ', num2str(priorityQ)])

priorityQ = priorityQ.';

congestion = computeEquilibriumBottleneck(priorityQ, congestion, population);
disp("the prelimilary solution is ")
disp(congestion.CostPriority)

% disp(['Prelimirary equlibirum costs are : '; num2str(congestion.CostPriority)])

% if max(congestion.CostPriority) > 0
%     [congestion.CostPriority,remained_capacities,new_t, convergence] = refining_costs(congestion.CostPriority, congestion, population, priorityQ,50,[]);
%     count = 1;
%     % if refining process is not converge, we try adding some noises in the costs 
%     while convergence > 1e-1 && count < 10
%         congestion.CostPriority = congestion.CostPriority + rand(length(congestion.CostPriority),1)/1000;
%         [congestion.CostPriority,remained_capacities,new_t, convergence] = refining_costs(congestion.CostPriority, congestion, population, priorityQ,50,[]);
%         disp('convergence after refining') 
%         disp(convergence) 
%         count = count +1;
%     end
%     congestion.accuratet = new_t;
%     [congestion.CostNonpriority,remained_capacities2,new_t2, convergence2] = refining_costs(congestion.CostNonpriority, congestion, population, 1-priorityQ,50,remained_capacities);
% %     if convergence2 > 1e-1
% %         congestion.CostNonpriority = congestion.CostPriority +rand(4,1)/1000;
% %         [congestion.CostNonpriority,remained_capacities2,new_t2, convergence2] = refining_costs(congestion.CostPriority, congestion, population, priorityQ,50,[]);
% %         disp(convergence)
% %     end
% end

TotalTripCost = dot(population.GroupSize .* priorityQ, congestion.CostPriority) + dot(population.GroupSize .* (1- priorityQ),congestion.CostNonpriority);

Trip_Costs_Priority = congestion.CostPriority;

Trip_Costs_Nonpriority = congestion.CostNonpriority;


