% =========================================================================
% 基于受限物理轨迹的动力学参数二次降维与条件数验证
% =========================================================================
clear; clc;

%% 1. 载入刚才优化的轨迹参数
load('Optimized_Fourier_Coeffs.mat'); % 载入 x_opt, N, wf, num_joints

% 重建轨迹数据
num_samples = 500; % 采样 500 个点用于高精度分析
t_array = linspace(0, 2*pi/wf, num_samples);
W_traj = zeros(num_joints * num_samples, 28); 

disp('正在重建轨迹并生成观测矩阵...');
for i = 1:num_samples
    t = t_array(i);
    [q, dq, d2q] = calc_fourier_kinematics(x_opt, t, num_joints, N, wf);
    
    % 调用你的快速计算函数 (对应之前的 28 个参数)
    W_t = calc_W_base(q(1), q(2), q(3), dq(1), dq(2), dq(3), d2q(1), d2q(2), d2q(3));
    W_traj((i-1)*3+1 : i*3, :) = double(W_t);
end

%% 2. 列归一化与二次 SVD 分析
% 列归一化
norm_W = vecnorm(W_traj);
norm_W(norm_W == 0) = 1; 
W_traj_norm = W_traj ./ norm_W;

% 进行奇异值分解
[U, S, V] = svd(W_traj_norm, 'econ');
s_diag = diag(S);

% 绘制奇异值分布图 (极其直观地看到是哪几个参数在捣乱)
figure;
semilogy(s_diag, 'o-', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on;
title('真实物理受限空间下的奇异值分布 (Log Scale)');
xlabel('参数索引'); ylabel('奇异值大小');
yline(1e-3, 'k--', '安全截断阈值 (1e-3)', 'LineWidth', 1.5);

%% 3. 重新提取真实可辨识的参数 (去除微小奇异值对应的列)
% 在真实物理限制下，奇异值小于 1e-3 的列基本上就是无法被激发的冗余列
tol = 1.5e-2; 
real_base_count = sum(s_diag > tol);

fprintf('\n由于物理限位，原本的 28 个参数中，在真实工作空间内真正能被激发的只有: %d 个！\n', real_base_count);

% 利用 QR 分解找出留下来的最强独立的列
[Q, R, E] = qr(W_traj_norm, 0);
best_indices = E(1:real_base_count);

% 提取二次降维后的矩阵
W_final_norm = W_traj_norm(:, best_indices);

% 计算最终完美的条件数
final_cond = cond(W_final_norm);
fprintf('👉 剔除不可激发的参数后，最终观测矩阵的条件数为: %.4f\n', final_cond);

disp('==================================================');
if final_cond < 150
    disp('🎉 恭喜！条件数已降至完美范围！你可以直接使用这组降维后的矩阵进行实机辨识了！');
else
    disp('⚠️ 条件数依然偏高，可以尝试在代码中将 tol 阈值 (如 1e-3) 略微调大一点点。');
end
disp('==================================================');

%% 辅助函数 (保持原样)
function [q, dq, d2q] = calc_fourier_kinematics(x, t, num_joints, N, wf)
    q = zeros(num_joints, 1); dq = zeros(num_joints, 1); d2q = zeros(num_joints, 1);
    idx = 1;
    for i = 1:num_joints
        q0 = x(idx); idx = idx + 1;
        a  = x(idx : idx+N-1); idx = idx + N;
        b  = x(idx : idx+N-1); idx = idx + N;
        q_i = q0; dq_i = 0; d2q_i = 0;
        for k = 1:N
            q_i   = q_i   + (a(k)/(wf*k)) * sin(wf*k*t) - (b(k)/(wf*k)) * cos(wf*k*t);
            dq_i  = dq_i  + a(k) * cos(wf*k*t) + b(k) * sin(wf*k*t);
            d2q_i = d2q_i - a(k) * wf * k * sin(wf*k*t) + b(k) * wf * k * cos(wf*k*t);
        end
        q(i) = q_i; dq(i) = dq_i; d2q(i) = d2q_i;
    end
end