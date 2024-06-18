% This function refines the preliminary costs from the LP. An MSA procedure is adopted here.

function [refined_costs, remained_capacities, new_t, convergence] = refining_costs(appox_costs_LP, congestion, population, priority_fraction, iterations, capacities)
Nstar = length(population.tstar); 
refined_costs = appox_costs_LP;

if isempty(capacities) % capacities is the maxium throughput of each intervel (dt)
    congestion.accuratet = congestion.t;
    Nt = length(congestion.accuratet);
    capacities = ones(Nt,1) * congestion.PreemptableS / congestion.S; % number of user per unit time step
else
    Nt = length(congestion.accuratet);
end
SP = zeros(Nt, Nstar);
for i = 1: Nstar
    SP(:,i)= -population.UD{i}(congestion.accuratet)-population.UO{i}(congestion.accuratet);
end

for i = 1:iterations
    stepsize = 1/2;
    [refined_costs, new_capacities, new_max_queueing, new_t, convergence, flag_valid] = tuning_costs(refined_costs, congestion, population, priority_fraction,SP,stepsize,capacities);
    if i == 1
        initial_convergence = convergence;
        disp("convergence——initial")
        disp(initial_convergence)
        if convergence > 1
        disp('bad quality at first')
        end
        temp_new_capacities = new_capacities;
        temp_new_max_queueing = new_max_queueing;
        temp_new_t = new_t;
    end
    if ~flag_valid 
        break
    end
end

disp("convergence——final")
disp(convergence)

if convergence > 1e+0
    disp('The obj function is NOT accurately evaluated')
    disp(-initial_convergence + convergence)
    
else
end
%% calculate the remained capacity for nonpriority users

% if it does not converge, the set it back to the original value
if abs(convergence) > initial_convergence
    disp('it does not converge')
    refined_costs = appox_costs_LP;
    used_capacities = temp_new_capacities;
    used_capacities(temp_new_max_queueing <0)=0;
    remained_capacities = temp_new_capacities * congestion.S / congestion.PreemptableS  - used_capacities;
    new_t = temp_new_t;
else
    used_capacities = new_capacities;
    used_capacities(new_max_queueing <0)=0;
    remained_capacities = new_capacities * congestion.S / congestion.PreemptableS  - used_capacities;
end

end

function Table_para = getTable(refined_costs, congestion, population, priority_fraction,t, SP, capacities)

% 1. Generating Isocost curves under the new appox_costs_LP
Nstar = length(population.tstar);
Nt = length(t);
queueing_time = zeros(Nt, Nstar);
for i = 1: Nstar
    queueing_time(:,i) = refined_costs(i) - SP(:,i);
end
% 2. Geting the arriving order and timing for each group
[max_queueing,idx_max_queueing] = max(queueing_time,[],2);
idx_max_queueing(max_queueing <= 0) = Nstar + 1; % add a virture group Nstar + 1
population.gamma = [population.gamma;0];
population.beta = [population.beta; 0];
population.tstar = [population.tstar, 0];
% find the time point of each interval starts at
diff_idx_max_queueing = idx_max_queueing - circshift(idx_max_queueing,1); % diff(t) > 0 if time t is the start
group_sequence = idx_max_queueing(diff_idx_max_queueing ~= 0); % j1, j2, ..., Nstar + 1
group_sequence = [Nstar + 1; group_sequence]; % Nstar + 1, j1, j2, ...,jn,
idx_time_points = sort(find(diff_idx_max_queueing));
time_points = t(idx_time_points);
point_number = length(time_points);
Table_para = zeros(point_number,5);
capacity_time_points = capacities(diff_idx_max_queueing ~= 0);% the maximum output for threshold points
% disp('capacity_time_points');
% disp(capacity_time_points);
%% maintain a table of parameters. It is based on the observation that, the marginal increase of the cost of group j depends on the parameter of j and the adjacent group of j. 
for i = 1: point_number
    Table_para(i,1) = group_sequence(i);  % the group index before this invervel
    Table_para(i,2) = group_sequence(i+1); % the group index in this invervel
    if population.tstar(group_sequence(i)) < time_points(i) % if the threshold point is after t^*_j
        Table_para(i,3) = - population.gamma(group_sequence(i));
    else
        Table_para(i,3) = population.beta(group_sequence(i));
    end
    if population.tstar(group_sequence(i+1)) < time_points(i)
        Table_para(i,4) = - population.gamma(group_sequence(i+1));
    else
        Table_para(i,4) = population.beta(group_sequence(i+1));
    end
    Table_para(i,5) = capacity_time_points(i) /(t(idx_time_points(i)+1) -time_points(i))/(Table_para(i,4) - Table_para(i,3));
end
end

%% refining the cost by comparing the estimated number of users and the true number of users 
function [tuned_costs, new_capacities, new_max_queueing, new_t, convergence, flag_valid] = tuning_costs(refined_costs, congestion, population, priority_fraction, SP, stepsize, capacities)

% Define the number of groups
Nstar = length(population.tstar);

% Compute the theoretical users for each group
N_theory =  population.GroupSize .* priority_fraction;

% Initialize the time
N_approx = zeros(Nstar,1);

% Initialize the time
t =  congestion.accuratet;
Nt = length(t); 

count = 0;
%% Increase the desity of discrete time points 
while count < 3 
    queueing_time = zeros(Nt, Nstar);
    for i = 1: Nstar
        queueing_time(:,i) = refined_costs(i) - SP(:,i);
    end

    % 2. Geting the arriving order and timing for each group
    [max_queueing, idx_max_queueing] = max(queueing_time,[],2);

    % add a virture group Nstar + 1
    idx_max_queueing(max_queueing <= 0) = Nstar + 1; 

    % find the time point of each interval starts at
    diff_idx_max_queueing = idx_max_queueing - circshift(idx_max_queueing,1); % diff(t) > 0 if time t is the start
    group_sequence = idx_max_queueing(diff_idx_max_queueing ~= 0); % j1, j2, ..., Nstar + 1
    group_sequence = [Nstar + 1; group_sequence]; % Nstar + 1, j1, j2, ...,jn,
    idx_time_points = sort(find(diff_idx_max_queueing));
    time_points = t(idx_time_points);
    
    %% Create new_t by breaking some time invervals into smaller intervals
    accuray2 = 100;
    new_t = t(1:idx_time_points(1) - 2);
    new_capacities = capacities(1:idx_time_points(1)-2);
    diif_newt = new_t - circshift(new_t, 1);
    % if ~isempty(diif_newt(diif_newt == 0))
    %     disp('something')
    %     disp( iter)
    %     disp(idx_time_points)
    %     disp(idx_max_queueing)
    % end
    % add points
    for i = 1:length(time_points) - 1
        idexi = idx_time_points(i);
        if idexi == 1
            disp('group_sequence')
            disp(group_sequence)
            % disp(queueing_time);
        end
        better_t = linspace(t(idexi-1),t(idexi),accuray2 + 1)';
        new_t = [new_t; better_t(1:end-1)];
        diif_newt = new_t - circshift(new_t, 1);
        if ~isempty(diif_newt(diif_newt == 0))
            disp('something')
            disp(idx_time_points)
            disp(idx_max_queueing)
        end        
        new_capacities = [new_capacities;capacities(idexi-1)/accuray2 * ones(accuray2,1)];
        middle_t = t(idexi:idx_time_points(i+1)-1);
        new_t = [new_t; middle_t(1: end-1)];
        diif_newt = new_t - circshift(new_t, 1);
         if ~isempty(diif_newt(diif_newt == 0))
            disp('something')
            disp(idx_time_points)
            disp(idx_max_queueing)
        end 
        middle_capacities = capacities(idexi:idx_time_points(i+1)-1);
        new_capacities = [new_capacities; middle_capacities(1: end-1)];
    end
    idexi = idx_time_points(length(time_points));
    better_t = linspace(t(idexi-1),t(idexi),accuray2+1)';  
    new_t = [new_t; better_t(1:end-1)];
    new_capacities = [new_capacities;capacities(idexi-1)/accuray2 * ones(accuray2,1)];
    middle_t = t(idexi:end);
    new_t = [new_t; middle_t];
    middle_capacities = capacities(idexi:end);
    new_capacities = [new_capacities; middle_capacities];
    capacities = new_capacities;
    diif_newt = new_t - circshift(new_t, 1);
    if ~isempty(diif_newt(diif_newt == 0))
        disp('something')
        disp(idx_time_points)
        disp(idx_max_queueing)
    end
    %% Update t and associated values with new_t 
    t = new_t;
    new_Nt = length(new_t);
    new_SP = zeros(new_Nt, Nstar);
    new_queueing_time = zeros(new_Nt, Nstar);
    for i = 1: Nstar
        new_SP(:,i)= - population.UD{i}(new_t) - population.UO{i}(new_t);
        new_queueing_time(:,i) = refined_costs(i) - new_SP(:,i);
    end
    SP = new_SP;
    Nt = new_Nt;
    count = count + 1;
end

%% Slightly adjust costs to ensure groups with small numbers are also detected 
count = 0;
flag_valid = true;
while any(N_theory > 1e-4 & N_approx <= 0 )
    %% Estimate N_approx 
    [new_max_queueing,new_idx_max_queueing] = max(new_queueing_time,[],2);
    new_idx_max_queueing(new_max_queueing <= 0) = Nstar + 1; % add a virture group Nstar + 1
    N_approx = zeros(Nstar,1);
    for i = 1: Nstar
        N_approx(i) = sum(new_capacities(new_idx_max_queueing == i));
    end
    count = count +1;
    % disp(count)
   if count > 1
        disp('N_approx')
        disp(N_approx)
        disp('refined_costs')
        disp(refined_costs)
        refined_costs(N_theory > 0 & N_approx <= 0) = refined_costs(N_theory > 0 & N_approx <= 0) + 1e-4;
    end
    if count > 20
       flag_valid  = false;
       break 
    end
end


%% Get parameters for the derivative 
tuned_costs = refined_costs;
convergence = sum(abs(N_approx - N_theory));
if flag_valid 
    marginal_shift = zeros(Nstar,1);
    Table_para = getTable(refined_costs, congestion, population, priority_fraction, new_t, new_SP, new_capacities);
    for i = 1: Nstar
        marginal_shift(i)  = sum(Table_para(Table_para(:,1) == i,5)) + sum(Table_para(Table_para(:,2) == i,5)) ;
    end
    
    %% Update the refine_cost
        gap_tau = N_theory - N_approx;
        % make all undetectable gaps to zero
        gap_tau(0 < gap_tau & gap_tau < 1e-4) = 0;
        delta_cost = gap_tau ./marginal_shift;
        delta_cost (isnan(delta_cost )) = 0;
    % Choose the stepsize to adjust the costs, ensuring the costs after adjustment leads to positive number of users for each group  
    while true 
        tuned_costs = refined_costs + stepsize * delta_cost;
    
        N_approx = estimate_N(tuned_costs, new_t, new_capacities, population);
        tem_convergence = sum(abs(N_approx - N_theory));
    %     disp("tem_convergence");
    %      disp(tem_convergence)
        if ~any(N_theory > 1e-4 & N_approx <=0)
            break
        end
        stepsize = 1/10 * stepsize; 
    end
end
end

function N_approx= estimate_N(refined_costs, new_t, new_capacities, population)
 % Initialize N_approx 
Nstar = length(population.tstar);
N_approx = zeros(Nstar,1);

new_Nt = length(new_t);
new_SP = zeros(new_Nt, Nstar);
new_queueing_time = zeros(new_Nt, Nstar);

for i = 1: Nstar
    new_SP(:,i)= - population.UD{i}(new_t) - population.UO{i}(new_t);
    new_queueing_time(:,i) = refined_costs(i) - new_SP(:,i);
end

[new_max_queueing,new_idx_max_queueing] = max(new_queueing_time,[],2);
new_idx_max_queueing(new_max_queueing <= 0) = Nstar + 1; % add a virture group Nstar + 1

for i = 1: Nstar
    N_approx(i) =  sum(new_capacities(new_idx_max_queueing == i));
end
end