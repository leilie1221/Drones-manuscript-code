function MonteCarlo_BNE_solve_greedy_heuristics
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
    F_star_D = ones(1,nD)/nD;
    F_star_A = ones(2,nA)/nA;
    expected_D_matrix = zeros(nD, nA);
    expected_A_matrix_high = zeros(nD, nA);
    expected_A_matrix_low = zeros(nD, nA);
    expected_A_matrix = zeros(nD, nA);
    for i = 1:nD
        for j = 1:nA
            expected_D_matrix(i, j) = p * high_benefit_matrix{i, j}(1) * F_star_D(i) + (1 - p) * low_benefit_matrix{i, j}(1) * F_star_D(i);
            expected_A_matrix_high(i, j) = p * high_benefit_matrix{i, j}(2) * F_star_A(j);
            expected_A_matrix_low(i, j) = (1 - p) * low_benefit_matrix{i, j}(2) * F_star_A(j);
            expected_A_matrix(i, j) = p * high_benefit_matrix{i, j}(2) * F_star_A(j) + (1 - p) * low_benefit_matrix{i, j}(2) * F_star_A(j);
        end
    end
    dominant_D_strategy = false(nD, 1);
    candidate_D = 1:nD;
    while length(candidate_D) > 1
        dominance_score = zeros(length(candidate_D), 1);
        for i = 1:length(candidate_D)
            for j = 1:length(candidate_D)
                if i ~= j

                    win_count = sum(expected_D_matrix(candidate_D(i), :) >= expected_D_matrix(candidate_D(j), :));
                    dominance_score(i) = dominance_score(i) + win_count;
                end
            end
        end

        [~, min_idx] = min(dominance_score);
        candidate_D(min_idx) = [];
    end
    if ~isempty(candidate_D)
        dominant_D_strategy(candidate_D) = true;
    end
    dominant_A_strategy_high = false(nA, 1);
    candidate_A_high = 1:nA;
    while length(candidate_A_high) > 1
        dominance_score = zeros(length(candidate_A_high), 1);
        for i = 1:length(candidate_A_high)
            for j = 1:length(candidate_A_high)
                if i ~= j

                    win_count = sum(expected_A_matrix_high(:, candidate_A_high(i)) >= expected_A_matrix_high(:, candidate_A_high(j)));
                    dominance_score(i) = dominance_score(i) + win_count;
                end
            end
        end
        [~, min_idx] = min(dominance_score);
        candidate_A_high(min_idx) = [];
    end
    if ~isempty(candidate_A_high)
        dominant_A_strategy_high(candidate_A_high) = true;
    end
    dominant_A_strategy_low = false(nA, 1);
    candidate_A_low = 1:nA;
    while length(candidate_A_low) > 1
        dominance_score = zeros(length(candidate_A_low), 1);
        for i = 1:length(candidate_A_low)
            for j = 1:length(candidate_A_low)
                if i ~= j
                    win_count = sum(expected_A_matrix_low(:, candidate_A_low(i)) >= expected_A_matrix_low(:, candidate_A_low(j)));
                    dominance_score(i) = dominance_score(i) + win_count;
                end
            end
        end
        [~, min_idx] = min(dominance_score);
        candidate_A_low(min_idx) = [];
    end
    if ~isempty(candidate_A_low)
        dominant_A_strategy_low(candidate_A_low) = true;
    end
    dominant_A_strategy = false(nA, 1);
    candidate_A = 1:nA;
    while length(candidate_A) > 1
        dominance_score = zeros(length(candidate_A), 1);
        for i = 1:length(candidate_A)
            for j = 1:length(candidate_A)
                if i ~= j
                    win_count = sum(expected_A_matrix(:, candidate_A(i)) >= expected_A_matrix(:, candidate_A(j)));
                    dominance_score(i) = dominance_score(i) + win_count;
                end
            end
        end
        [~, min_idx] = min(dominance_score);
        candidate_A(min_idx) = [];
    end
    if ~isempty(candidate_A)
        dominant_A_strategy(candidate_A) = true;
    end
    if any(dominant_D_strategy)
        D_opt_strategy = find(dominant_D_strategy, 1);
    else


        A_strategy_index = find(dominant_A_strategy, 1);
        D_subopt_benefit_matrix = expected_A_matrix(:,A_strategy_index);
        [~,D_opt_strategy] = max(D_subopt_benefit_matrix);

        dominant_D_strategy(D_opt_strategy) = 1;
    end
    if any(dominant_A_strategy_high)
        A_opt_strategy_high = find(dominant_A_strategy_high, 1);
    else
        D_strategy_index = find(dominant_D_strategy, 1);
        A_subopt_benefit_matrix_high = expected_A_matrix_high(D_strategy_index,:);
        [~,A_opt_strategy_high] = max(A_subopt_benefit_matrix_high);

        dominant_A_strategy_high(A_opt_strategy_high) = 1;
    end
    if any(dominant_A_strategy_low)
        A_opt_strategy_low = find(dominant_A_strategy_low, 1);
    else
        D_strategy_index = find(dominant_D_strategy, 1);
        A_subopt_benefit_matrix_low = expected_A_matrix_low(D_strategy_index,:);
        [~,A_opt_strategy_low] = max(A_subopt_benefit_matrix_low);
        dominant_A_strategy_low(A_opt_strategy_low) = 1;
    end
    efficiency_D = zeros(1,nD);
    efficiency_A = zeros(1,nA);
    F_star_A = [dominant_A_strategy_high';dominant_A_strategy_low'];
    F_star_D = dominant_D_strategy;
    for k = 1:nD
        efficiency_D(k) = calculate_efficiency_ver2(num_type, P_A, F_star_A, U, k, 'defender');
    end
    for k = 1:nA
        efficiency_A(k) = calculate_efficiency_ver2(num_type, P_D, F_star_D, U, k, 'attacker');
    end
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
