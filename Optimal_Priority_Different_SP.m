%% This script is coded for finding the optimal priority fraction $\bm q$ for each group with the SBO-based method 

clear

clc

% Set up the options of surrogateopt. 

options = optimoptions("surrogateopt","MaxFunctionEvaluations",200,'UseParallel',true,'PlotFcn', 'optimplotx' );

% Load the parameters used in Case 3 introduced in Section 3.3

Filename_Parameters = append('./','parameters_case3');

Parameters_Preferences  = load(Filename_Parameters);

% Specify the step size of the preemptive capacity $S^P$  
stepSizeSP = 0.02;

% Initialization 

[NumberGroup,~] = size(Parameters_Preferences.parameters);

resultArrayQ = zeros(1/stepSizeSP, NumberGroup);

resultObj = zeros(1/stepSizeSP, 1);

Total_demand = 6000; 

Mechanism = 'Continuous'; % It depends on if $q_j$ is 'Continuous' or 'Discrete'

rng default % For reproducibility

%% Create the population
tstar = Parameters_Preferences.parameters(:,1);
alpha  = Parameters_Preferences.parameters(:,3);
beta =  Parameters_Preferences.parameters(:,4);
gamma =  Parameters_Preferences.parameters(:,5);
population = generateSParctan(tstar, alpha - beta, alpha + gamma, alpha, 1e100 * ones(length(tstar),1));%
population.GroupSize = Total_demand * Parameters_Preferences.parameters(:,2);
population.alpha = alpha;
population.gamma = gamma;
population.beta = beta;

%% Define discretized time 
Capacity = 1/6; 
PreemptableCapacity = 0;
congestion = generateBottleneck(Capacity,PreemptableCapacity);
Nstar = length(population.tstar); % The number of groups 
tmin = min(population.tstar) - 1/congestion.S; % hours 
tmax = max(population.tstar) + 1/congestion.S;
dt = 1/(sum(population.GroupSize) * congestion.S); % choose a time interval duration such that at most only one user pass each time
t = (tmin:dt:tmax)';
congestion.t = t;
congestion.dt = dt;
Nt = length(t);                                      
SP_alpha = zeros(Nt,Nstar);
for i = 1:Nstar
    SP_alpha(:,i)= (-population.UD{i}(t) - population.UO{i}(t))/population.alpha(i) ;
    SP_real(:,i)= (-population.UD{i}(t) - population.UO{i}(t)) ;
end
congestion.SP = SP_alpha;
congestion.SP_real = SP_real;


%% Find the optimal q given different preempable capacities
k = 0;

for PreemptRatio = 0.8 % SP / S
% for PreemptRatio = 0.1% SP / S

    PreemptableCapacity = PreemptRatio * Capacity;
    congestion.PreemptableS = PreemptableCapacity;
    % Constructing the surrogate optimization problem 
    lb = zeros(1, NumberGroup);
    ub = ones(1, NumberGroup);
    A = ones(1, NumberGroup)/NumberGroup;
    b = PreemptRatio * 1;
    Total_Trip_Cost = @(x)(Get_Total_Trip_Cost_refine(x, PreemptRatio,population, congestion)); 
    switch Mechanism
        case 'Continuous'
             intcon = [];
        case 'Discrete'
             intcon = 1:NumberGroup;
    end
    [sol,fval,exitflag,output,trials] = surrogateopt(Total_Trip_Cost...
        ,lb,ub,intcon,A,b,[],[],options);
    k = k + 1;
    resultArrayQ(k,:) = sol(1,:);
    resultObj(k,:) = fval;
end 
%% Save the results to a .mat file

suffix_elements = 'a' : 'z';
suffix = suffix_elements (randi(26,5,1));
filename = append('./NumericalResultsData/','resultfile',suffix,'.mat');
prompt = "Do you need saving the data? Y/N";
txt = input(prompt,"s");
if isempty(txt)
    txt = 'Y';
end
if txt == 'Y'
   save(filename,'resultArrayQ','resultObj','Filename_Parameters','options','population');
end
