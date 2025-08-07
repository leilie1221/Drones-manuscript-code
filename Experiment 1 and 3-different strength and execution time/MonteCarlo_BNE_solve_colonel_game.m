
function MonteCarlo_BNE_solve_colonel_game
global case_flag
monte_carlo_max = 200;
D_strategy_case1        = zeros(monte_carlo_max,4);
D_strategy_case2        = zeros(monte_carlo_max,4);
D_strategy_case3        = zeros(monte_carlo_max,4);
D_strategy_case4        = zeros(monte_carlo_max,4);
A_low_strategy_case1    = zeros(monte_carlo_max,4);
A_low_strategy_case2    = zeros(monte_carlo_max,4);
A_low_strategy_case3    = zeros(monte_carlo_max,4);
A_low_strategy_case4    = zeros(monte_carlo_max,4);
A_high_strategy_case1   = zeros(monte_carlo_max,4);
A_high_strategy_case2   = zeros(monte_carlo_max,4);
A_high_strategy_case3   = zeros(monte_carlo_max,4);
A_high_strategy_case4   = zeros(monte_carlo_max,4);
efficiency_D_case1      = zeros(monte_carlo_max,4);
efficiency_D_case2      = zeros(monte_carlo_max,4);
efficiency_D_case3      = zeros(monte_carlo_max,4);
efficiency_D_case4      = zeros(monte_carlo_max,4);
efficiency_A_case1      = zeros(monte_carlo_max,4);
efficiency_A_case2      = zeros(monte_carlo_max,4);
efficiency_A_case3      = zeros(monte_carlo_max,4);
efficiency_A_case4      = zeros(monte_carlo_max,4);
time_record_case1       = zeros(monte_carlo_max,1);
time_record_case2       = zeros(monte_carlo_max,1);
time_record_case3       = zeros(monte_carlo_max,1);
time_record_case4       = zeros(monte_carlo_max,1);
for monte_carlo_num = 1:monte_carlo_max
    for case_num = 1:4
        switch case_num
            case 1
                case_flag = 'case1';
            case 2
                case_flag = 'case2';
            case 3
                case_flag = 'case3';
            case 4
                case_flag = 'case4';
        end
        switch case_flag
            case 'case1'
                data_temp = load('benefit_data_matrix_case1.mat');
            case 'case2'
                data_temp = load('benefit_data_matrix_case2.mat');
            case 'case3'
                data_temp = load('benefit_data_matrix_case3.mat');
            case 'case4'
                data_temp = load('benefit_data_matrix_case4.mat');
        end
        benefit_data_matrix = data_temp.benefit_data_matrix;
        attri_x_rate = data_temp.attri_x_rate;
        attri_y_rate = data_temp.attri_y_rate;
        D_matrix_high = attri_x_rate(:,2);
        D_matrix_low  = attri_x_rate(:,1);
        A_matrix_high = attri_y_rate(:,2);
        A_matrix_low  = attri_y_rate(:,1);
        D_matrix_high = reshape(D_matrix_high,[4,4]);
        D_matrix_low  = reshape(D_matrix_low,[4,4]);
        A_matrix_high = reshape(A_matrix_high,[4,4]);
        A_matrix_low  = reshape(A_matrix_low,[4,4]);
        win_matrix_high = double(D_matrix_high < A_matrix_high);
        win_matrix_low  = double(D_matrix_low < A_matrix_low);
        p = 0.45;
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
        expected_D_attri_matrix      = p * ( win_matrix_high*2 - 1 ) + (1 - p) * ( win_matrix_low*2 - 1 );
        expected_A_attri_matrix_high = -(p * ( win_matrix_high*2 - 1 ));
        expected_A_attri_matrix_low  = -((1 - p) * ( win_matrix_low*2 - 1 ));
        for i = 1:nD
            for j = 1:nA
                expected_D_matrix(i, j) = p * high_benefit_matrix{i, j}(1) * F_star_D(i) + (1 - p) * low_benefit_matrix{i, j}(1) * F_star_D(i);
                expected_A_matrix_high(i, j) = p * high_benefit_matrix{i, j}(2) * F_star_A(j);
                expected_A_matrix_low(i, j) = (1 - p) * low_benefit_matrix{i, j}(2) * F_star_A(j);
                expected_A_matrix(i, j) = p * high_benefit_matrix{i, j}(2) * F_star_A(j) + (1 - p) * low_benefit_matrix{i, j}(2) * F_star_A(j);
            end
        end
        tic;
        f_D = ones(nD,1);
        A_D = -expected_D_attri_matrix';
        b_D = -ones(nD,1);
        Aeq_D = [];
        beq_D = [];
        lb_D = zeros(1,nD);
        ub_D = ones(1,nD);
        ctype = "UUUU";
        [x_D, fval_D] = glpk(f_D, A_D, b_D, lb_D, ub_D, ctype);
        F_star_D = x_D'/fval_D;
        f_A = -ones(nA,1);
        A_A_high = expected_A_attri_matrix_high;
        A_A_low = expected_A_attri_matrix_low;
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
        if isna(F_star_D)
            [~,D_opt_index] = max(sum(expected_D_matrix,2));
            F_star_D = zeros(1,nD);
            F_star_D(D_opt_index) = 1;
        end
        if isna(F_star_A_high)
            [~,A_opt_high_index] = max(sum(expected_A_matrix_high,1));
            F_star_A_high = zeros(1,nD);
            F_star_A_high(A_opt_high_index) = 1;
            F_star_A(1,:) = F_star_A_high;
        end
        if isna(F_star_A_low)
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
        switch case_flag
            case 'case1'
                D_strategy_case1(monte_carlo_num,:)         = F_star_D;
                A_low_strategy_case1(monte_carlo_num,:)     = F_star_A_low;
                A_high_strategy_case1(monte_carlo_num,:)    = F_star_A_high;
                efficiency_D_case1(monte_carlo_num,:)       = efficiency_D;
                efficiency_A_case1(monte_carlo_num,:)       = efficiency_A;
                time_record_case1(monte_carlo_num)          = toc;
            case 'case2'
                D_strategy_case2(monte_carlo_num,:)         = F_star_D;
                A_low_strategy_case2(monte_carlo_num,:)     = F_star_A_low;
                A_high_strategy_case2(monte_carlo_num,:)    = F_star_A_high;
                efficiency_D_case2(monte_carlo_num,:)       = efficiency_D;
                efficiency_A_case2(monte_carlo_num,:)       = efficiency_A;
                time_record_case2(monte_carlo_num)          = toc;
            case 'case3'
                D_strategy_case3(monte_carlo_num,:)         = F_star_D;
                A_low_strategy_case3(monte_carlo_num,:)     = F_star_A_low;
                A_high_strategy_case3(monte_carlo_num,:)    = F_star_A_high;
                efficiency_D_case3(monte_carlo_num,:)       = efficiency_D;
                efficiency_A_case3(monte_carlo_num,:)       = efficiency_A;
                time_record_case3(monte_carlo_num)          = toc;
            case 'case4'
                D_strategy_case4(monte_carlo_num,:)         = F_star_D;
                A_low_strategy_case4(monte_carlo_num,:)     = F_star_A_low;
                A_high_strategy_case4(monte_carlo_num,:)    = F_star_A_high;
                efficiency_D_case4(monte_carlo_num,:)       = efficiency_D;
                efficiency_A_case4(monte_carlo_num,:)       = efficiency_A;
                time_record_case4(monte_carlo_num)          = toc;
        end
    end
disp(['仿真进度：',num2str(monte_carlo_num),'/',num2str(monte_carlo_max)]);
end
D_strategy_case1_average        = sum(D_strategy_case1,1)       / monte_carlo_max;
A_low_strategy_case1_average    = sum(A_low_strategy_case1,1)   / monte_carlo_max;
A_high_strategy_case1_average   = sum(A_high_strategy_case1,1)  / monte_carlo_max;
efficiency_D_case1_average      = sum(efficiency_D_case1,1)     / monte_carlo_max;
efficiency_A_case1_average      = sum(efficiency_A_case1,1)     / monte_carlo_max;
D_strategy_case2_average        = sum(D_strategy_case2,1)       / monte_carlo_max;
A_low_strategy_case2_average    = sum(A_low_strategy_case2,1)   / monte_carlo_max;
A_high_strategy_case2_average   = sum(A_high_strategy_case2,1)  / monte_carlo_max;
efficiency_D_case2_average      = sum(efficiency_D_case2,1)     / monte_carlo_max;
efficiency_A_case2_average      = sum(efficiency_A_case2,1)     / monte_carlo_max;
D_strategy_case3_average        = sum(D_strategy_case3,1)       / monte_carlo_max;
A_low_strategy_case3_average    = sum(A_low_strategy_case3,1)   / monte_carlo_max;
A_high_strategy_case3_average   = sum(A_high_strategy_case3,1)  / monte_carlo_max;
efficiency_D_case3_average      = sum(efficiency_D_case3,1)     / monte_carlo_max;
efficiency_A_case3_average      = sum(efficiency_A_case3,1)     / monte_carlo_max;
D_strategy_case4_average        = sum(D_strategy_case4,1)       / monte_carlo_max;
A_low_strategy_case4_average    = sum(A_low_strategy_case4,1)   / monte_carlo_max;
A_high_strategy_case4_average   = sum(A_high_strategy_case4,1)  / monte_carlo_max;
efficiency_D_case4_average      = sum(efficiency_D_case4,1)     / monte_carlo_max;
efficiency_A_case4_average      = sum(efficiency_A_case4,1)     / monte_carlo_max;
time_record_case1_average       = sum(time_record_case1)        / monte_carlo_max;
time_record_case2_average       = sum(time_record_case2)        / monte_carlo_max;
time_record_case3_average       = sum(time_record_case3)        / monte_carlo_max;
time_record_case4_average       = sum(time_record_case4)        / monte_carlo_max;
time_record_average             = (time_record_case1_average+time_record_case2_average+...
                                  time_record_case3_average+time_record_case4_average)/4;
disp('---------------------- Case 1 ----------------------')
disp(['Equilibrium of high-tech A: ',num2str(A_high_strategy_case1_average)]);
disp(['Equilibrium of low-tech A: ',num2str(A_low_strategy_case1_average)]);
disp(['Equilibrium of D: ',num2str(D_strategy_case1_average)]);
disp(['Effectiveness of D: ',num2str(efficiency_D_case1_average)]);
disp('---------------------- Case 2 ----------------------')
disp(['Equilibrium of high-tech A: ',num2str(A_high_strategy_case2_average)]);
disp(['Equilibrium of low-tech A: ',num2str(A_low_strategy_case2_average)]);
disp(['Equilibrium of D: ',num2str(D_strategy_case2_average)]);
disp(['Effectiveness of D: ',num2str(efficiency_D_case2_average)]);
disp('---------------------- Case 3 ----------------------')
disp(['Equilibrium of high-tech A: ',num2str(A_high_strategy_case3_average)]);
disp(['Equilibrium of low-tech A: ',num2str(A_low_strategy_case3_average)]);
disp(['Equilibrium of D: ',num2str(D_strategy_case3_average)]);
disp(['Effectiveness of D: ',num2str(efficiency_D_case3_average)]);
disp('---------------------- Case 4 ----------------------')
disp(['Equilibrium of high-tech A: ',num2str(A_high_strategy_case4_average)]);
disp(['Equilibrium of low-tech A: ',num2str(A_low_strategy_case4_average)]);
disp(['Equilibrium of D: ',num2str(D_strategy_case4_average)]);
disp(['Effectiveness of D: ',num2str(efficiency_D_case4_average)]);
disp('----------------------------------------------------')

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
