clc;clear all;close all;

%% Run code
global p
for p = 0.1:0.1:0.9
 fprintf('================== p=%.1f ==========================\n',p);
    disp('================== SSABG ==========================')
    MonteCarlo_BNE_solve_SSABG;

    disp('============== Random allocation ==================')
    MonteCarlo_BNE_solve_random_allocation;

    disp('============== Greedy heruistics ==================')
    MonteCarlo_BNE_solve_greedy_heuristics;

    disp('============= Rule-based assignment ===============')
    MonteCarlo_BNE_solve_ruled_assignment;

    disp('================= Colonel game ====================')
    MonteCarlo_BNE_solve_colonel_game
end
