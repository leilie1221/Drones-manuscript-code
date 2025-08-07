function MonteCarlo_BNE_solve_random_allocation
global p
monte_carlo_max = 200;
D_strategy_case        = zeros(monte_carlo_max,4);
A_low_strategy_case    = zeros(monte_carlo_max,4);
A_high_strategy_case   = zeros(monte_carlo_max,4);
efficiency_D_case      = zeros(monte_carlo_max,4);
efficiency_A_case      = zeros(monte_carlo_max,4);
time_record             = zeros(monte_carlo_max,1);
for monte_carlo_num = 1:monte_carlo_max
    tic;
    data_temp = load('benefit_data_matrix_case4.mat');
    benefit_data_matrix = data_temp.benefit_data_matrix;
    high_benefit_matrix = benefit_data_matrix(:,1:4);
    low_benefit_matrix = benefit_data_matrix(:,5:8);
    U = cat(3,...
        high_benefit_matrix,...
        low_benefit_matrix);
    nD = size(high_benefit_matrix, 1);
    nA = size(high_benefit_matrix, 2);
    P_A = [p,1-p];
    P_D = [1,0];
    num_type_A = 2;
    num_type_D = 1;
    num_type = [num_type_A, num_type_D];
    N = 1e2;
    best_D   = 0;
    best_Ah  = 0;
    best_Al  = 0;
    max_D    = -inf;
    max_Ah   = -inf;
    max_Al   = -inf;
    rng('shuffle');
    for k = 1:N
        d_idx = randi(nD);
        ah_idx = randi(nA);
        al_idx = randi(nA);
        F_D = zeros(1,nD);
        F_D(d_idx) = 1;
        F_A = zeros(2,nA);
        F_A(1,ah_idx) = 1;
        F_A(2,al_idx) = 1;
        val_D = 0;
        for t = 1:2
            val_D = val_D + P_A(t) * U{d_idx,ah_idx,t}(1) * p * (t==1) + ...
                             P_A(t) * U{d_idx,al_idx,t}(1) * (1-p) * (t==2);
        end
        val_Ah = 0;
        for d = 1:nD
            val_Ah = val_Ah + P_D(1) * U{d,ah_idx,1}(2) * F_D(d);
        end
        val_Al = 0;
        for d = 1:nD
            val_Al = val_Al + P_D(1) * U{d,al_idx,2}(2) * F_D(d);
        end
        if val_D  > max_D
            max_D  = val_D;
            best_D  = d_idx;
        end
        if val_Ah > max_Ah
            max_Ah = val_Ah;
            best_Ah = ah_idx;
        end
        if val_Al > max_Al
            max_Al = val_Al;
            best_Al = al_idx;
        end
    end
    F_star_D          = zeros(1,nD);
    F_star_A_high     = zeros(1,nA);
    F_star_A_low      = zeros(1,nA);
    F_star_D(best_D)  = 1;
    F_star_A_high(best_Ah) = 1;
    F_star_A_low(best_Al)  = 1;
    F_star_A          = [F_star_A_high; F_star_A_low];
    efficiency_D = zeros(1,nD);
    efficiency_A = zeros(1,nA);
    for k = 1:nD
        efficiency_D(k) = calculate_efficiency_ver2(num_type, P_A, F_star_A, U, k, 'defender');
    end
    for k = 1:nA
        efficiency_A(k) = calculate_efficiency_ver2(num_type, P_D, F_star_D, U, k, 'attacker');
    end
    D_opt_strategy = find(F_star_D);
    A_opt_strategy_high = find(F_star_A(1,:));
    A_opt_strategy_low = find(F_star_A(2,:));
    A_strategies = cell(1, nA);
    D_strategies = cell(1, nD);
    for i = 1:nA
        A_strategies{i} = sprintf('a%d', i);
    end
    for i = 1:nD
        D_strategies{i} = sprintf('d%d', i);
    end
    optimal_A_strategy_high = A_strategies(A_opt_strategy_high);
    optimal_A_strategy_low = A_strategies(A_opt_strategy_low);
    optimal_D_strategy = D_strategies(D_opt_strategy);
    D_strategy_case(monte_carlo_num,:)         = F_star_D;
    A_low_strategy_case(monte_carlo_num,:)     = F_star_A(1,:);
    A_high_strategy_case(monte_carlo_num,:)    = F_star_A(2,:);
    efficiency_D_case(monte_carlo_num,:)       = efficiency_D;
    efficiency_A_case(monte_carlo_num,:)       = efficiency_A;
time_record(monte_carlo_num) = toc;
end
D_strategy_case_average        = sum(D_strategy_case,1)       / monte_carlo_max;
A_low_strategy_case_average    = sum(A_low_strategy_case,1)   / monte_carlo_max;
A_high_strategy_case_average   = sum(A_high_strategy_case,1)  / monte_carlo_max;
efficiency_D_case_average      = sum(efficiency_D_case,1)     / monte_carlo_max;
efficiency_A_case_average      = sum(efficiency_A_case,1)     / monte_carlo_max;
time_record_average             = sum(time_record)              / monte_carlo_max;
disp(['Equilibrium of high-tech A: ',num2str(A_high_strategy_case_average)]);
disp(['Equilibrium of low-tech A: ',num2str(A_low_strategy_case_average)]);
disp(['Equilibrium of D: ',num2str(D_strategy_case_average)]);
disp(['Effectiveness of D: ',num2str(efficiency_D_case_average)]);
function efficiency = calculate_efficiency_ver2(num_type, P, F_star, U, strategy, player)
    efficiency = 0;
    num_type_A = num_type(1);
    num_type_D = num_type(2);
    if strcmp(player,'defender')
        nA = length(F_star);
        player_mode = 1;
        for t_A = 1:num_type_A
            for m_A_i = 1:nA
                efficiency = efficiency + P(t_A) * U{strategy, m_A_i, t_A}(player_mode) * F_star(t_A,m_A_i);
            end
        end
    elseif strcmp(player,'attacker')
        nD = length(F_star);
        player_mode = 2;
        for t_A = 1:num_type_A
            for s_D_j = 1:nD
                efficiency = efficiency + P(num_type_D) * U{s_D_j, strategy, t_A}(player_mode) * F_star(s_D_j);
            end
        end
    end
end

end
