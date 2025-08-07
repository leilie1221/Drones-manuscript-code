function MonteCarlo_BNE_solve_SSABG
global p
monte_carlo_max = 200;
D_strategy_case        = zeros(monte_carlo_max,4);
A_low_strategy_case    = zeros(monte_carlo_max,4);
A_high_strategy_case   = zeros(monte_carlo_max,4);
efficiency_D_case      = zeros(monte_carlo_max,4);
efficiency_A_case      = zeros(monte_carlo_max,4);
time_record             = zeros(monte_carlo_max,1);
for monte_carlo_num = 1:monte_carlo_max
    output_D      = [];
    output_A_high = [];
    output_A_low  = [];
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
    for i = 1:nD
        is_dominant = true;
        for j = 1:nD
            if i ~= j
                for k = 1:nA
                    if expected_D_matrix(i, k) < expected_D_matrix(j, k)
                        is_dominant = false;
                        break;
                    end
                end
                if ~is_dominant
                    break;
                end
            end
        end
        dominant_D_strategy(i) = is_dominant;
    end
    dominant_A_strategy_high = false(nA, 1);
    for i = 1:nA
        is_dominant = true;
        for j = 1:nA
            if i ~= j
                for k = 1:nD
                    if expected_A_matrix_high(k, i) < expected_A_matrix_high(k, j)
                        is_dominant = false;
                        break;
                    end
                end
                if ~is_dominant
                    break;
                end
            end
        end

        dominant_A_strategy_high(i) = is_dominant;
    end
    dominant_A_strategy_low = false(nA, 1);
    for i = 1:nA
        is_dominant = true;
        for j = 1:nA
            if i ~= j
                for k = 1:nD
                    if expected_A_matrix_low(k, i) < expected_A_matrix_low(k, j)
                        is_dominant = false;
                        break;
                    end
                end
                if ~is_dominant
                    break;
                end
            end
        end
        dominant_A_strategy_low(i) = is_dominant;
    end
    dominant_A_strategy = false(nA, 1);
    for i = 1:nA
        is_dominant = true;
        for j = 1:nA
            if i ~= j
                for k = 1:nD
                    if expected_A_matrix(k, i) < expected_A_matrix(k, j)
                        is_dominant = false;
                        break;
                    end
                end
                if ~is_dominant
                    break;
                end
            end
        end
        dominant_A_strategy(i) = is_dominant;
    end
    if ~any(dominant_D_strategy) && ~any(dominant_A_strategy_high) && ~any(dominant_A_strategy_low)
        f_D = ones(nD,1);
        A_D = -expected_D_matrix';
        b_D = -ones(nD,1);
        Aeq_D = [];
        beq_D = [];
        lb_D = zeros(1,nD);
        ub_D = ones(1,nD);
        ctype="UUUU";
        [x_D, fval_D] = glpk(f_D, A_D, b_D, lb_D, ub_D, ctype);
        F_star_D = x_D'/fval_D;
        f_A = -ones(nA,1);
        A_A_high = expected_A_matrix_high;
        A_A_low = expected_A_matrix_low;
        b_A = ones(nA,1);
        Aeq_A = [];
        beq_A = [];
        lb_A = zeros(1,nA);
        ub_A = ones(1,nD);
        [x_A_high, fval_A_high] = glpk(f_A, A_A_high, b_A, lb_A, ub_A, ctype);
        [x_A_low, fval_A_low] = glpk(f_A, A_A_low, b_A, lb_A, ub_A, ctype);
        F_star_A_high = x_A_high'/-fval_A_high;
        F_star_A_low = x_A_low'/-fval_A_low;
        F_star_A = [F_star_A_high;F_star_A_low];
        if isempty(F_star_D)
            [~,D_opt_index] = max(sum(expected_D_matrix,2));
            F_star_D = zeros(1,nD);
            F_star_D(D_opt_index) = 1;
        end
        if isempty(F_star_A_high)
            [~,A_opt_high_index] = max(sum(expected_A_matrix_high,1));
            F_star_A_high = zeros(1,nD);
            F_star_A_high(A_opt_high_index) = 1;
            F_star_A(1,:) = F_star_A_high;
        end
        if isempty(F_star_A_low)
            [~,A_opt_low_index] = max(sum(expected_A_matrix_high,1));
            F_star_A_low = zeros(1,nD);
            F_star_A_low(A_opt_low_index) = 1;
            F_star_A(2,:) = F_star_A_low;
        end
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
    else

        if any(dominant_D_strategy)
            D_opt_strategy = find(dominant_D_strategy, 1);
        else


            A_strategy_index = find(dominant_A_strategy, 1);
            D_subopt_benefit_matrix = expected_A_matrix(:,A_strategy_index);
            [~,D_opt_strategy] = max(D_subopt_benefit_matrix);
            dominant_D_strategy = ones(1,4)/nD;
        end
        if any(dominant_A_strategy_high)
            A_opt_strategy_high = find(dominant_A_strategy_high, 1);
        else
            D_strategy_index = find(dominant_D_strategy, 1);
            A_subopt_benefit_matrix_high = expected_A_matrix_high(D_strategy_index,:);
            [~,A_opt_strategy_high] = max(A_subopt_benefit_matrix_high);
            dominant_A_strategy_high = ones(4,1)/nA;
        end
        if any(dominant_A_strategy_low)
            A_opt_strategy_low = find(dominant_A_strategy_low, 1);
        else
            D_strategy_index = find(dominant_D_strategy, 1);
            A_subopt_benefit_matrix_low = expected_A_matrix_low(D_strategy_index,:);
            [~,A_opt_strategy_low] = max(A_subopt_benefit_matrix_low);
            dominant_A_strategy_low = ones(4,1)/nA;
        end
        efficiency_D = zeros(1,nD);
        efficiency_A = zeros(1,nA);
        F_star_A = [dominant_A_strategy_high';dominant_A_strategy_low'];
        F_star_A_high = dominant_A_strategy_high';
        F_star_A_low = dominant_A_strategy_low';
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
    end
    F_star_D(isnan(F_star_D)) = 0;
    F_star_A_low(isnan(F_star_A_low)) = 0;
    F_star_A_high(isnan(F_star_A_high)) = 0;
    efficiency_D(isnan(efficiency_D)) = 0;
    efficiency_A(isnan(efficiency_A)) = 0;
    if isempty(output_D)
        output_D.iterations = 0;
    end
    if isempty(output_A_high)
        output_A_high.iterations = 0;
    end
    if isempty(output_A_low)
        output_A_low.iterations = 0;
    end
    D_strategy_case(monte_carlo_num,:)         = F_star_D;
    A_low_strategy_case(monte_carlo_num,:)     = F_star_A_low;
    A_high_strategy_case(monte_carlo_num,:)    = F_star_A_high;
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
