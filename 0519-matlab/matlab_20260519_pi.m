% finite_element_pi_extrapolation.m
% 有限元法基本概念示例：正多边形逼近圆计算 pi，展示外推对精度的提升
clear; clc; close all;

%% 1. 参数设置
n_list = 2.^(0:10);          % n = 1,2,4,8,...,1024
h_list = 1 ./ n_list;         % h = 1/n
pi_exact = pi;

%% 2. 计算原始 pi_n 及绝对误差
pi_n = zeros(size(n_list));
error = zeros(size(n_list));
for i = 1:length(n_list)
    n = n_list(i);
    pi_n(i) = n * sin(pi / n);
    error(i) = abs(pi_exact - pi_n(i));
end

%% 3. 显示部分结果表格
fprintf('   n\t\t pi_n\t\t\t 误差\n');
fprintf('--------------------------------------------\n');
for i = 1:length(n_list)
    fprintf('%4d\t %.15f\t %.2e\n', n_list(i), pi_n(i), error(i));
end

%% 4. 收敛阶分析（原始误差）
idx = n_list >= 16;
p = polyfit(log(h_list(idx)), log(error(idx)), 1);
order_orig = p(1);
fprintf('\n原始收敛阶 = %.4f (理论值 2)\n', order_orig);

%% 5. Richardson 外推（n 与 2n 配对）
% 外推公式: pi_rich = (4*pi_{2n} - pi_n)/3
n_rich = n_list(2:end);                % 配对所用的较小 n（实际对应较密的网格）
h_rich = 1 ./ (2 * n_rich);            % 外推后的等效 h（对应 2n 的网格尺寸）
pi_rich = (4 * pi_n(2:end) - pi_n(1:end-1)) / 3;
error_rich = abs(pi_exact - pi_rich);

fprintf('\nRichardson 外推结果（n 与 2n 配对）:\n');
fprintf('   n\t\t 外推值\t\t\t 误差\n');
fprintf('--------------------------------------------\n');
for i = 1:length(n_rich)
    fprintf('%4d\t %.15f\t %.2e\n', n_rich(i), pi_rich(i), error_rich(i));
end

%% 6. Wynn-ε 外推（基于最后几个 pi_n）
seq_short = pi_n(end-3:end);   % n = 128,256,512,1024
pi_wynn = wynn_epsilon(seq_short(:));
error_wynn = abs(pi_exact - pi_wynn);
fprintf('\nWynn-ε 外推（最后4个 pi_n）:\n');
fprintf('外推值 = %.15f, 误差 = %.2e\n', pi_wynn, error_wynn);

%% 7. 绘图：原始误差 + Richardson 外推误差 + Wynn-ε 参考线
figure;
loglog(h_list, error, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
loglog(h_rich, error_rich, 'rs--', 'LineWidth', 1.5, 'MarkerSize', 6);
% 绘制 Wynn-ε 外推的误差水平线（从最小 h 到最大 h）
xlim_min = min(h_list(end), h_rich(end));
xlim_max = max(h_list(1), h_rich(1));
x_hline = [xlim_min, xlim_max];
loglog(x_hline, [error_wynn, error_wynn], 'k--', 'LineWidth', 1.2);
% 可选：在 Wynn-ε 水平线上标注一个点，用于图例
% 在合适位置画一个空点，仅为了图例显示
loglog(NaN, NaN, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Wynn-ε 外推误差');

% 添加拟合直线（原始误差渐近线）
h_fit = logspace(log10(min(h_list(idx))), log10(max(h_list(idx))), 50);
error_fit = exp(polyval(p, log(h_fit)));
loglog(h_fit, error_fit, 'c-.', 'LineWidth', 1.2);

xlabel('h = 1/n (或 h_{rich} = 1/(2n) )');
ylabel('绝对误差');
title('不同逼近方式的误差对比');
legend('原始 π_n', 'Richardson 外推', 'Wynn-ε 外推误差', ...
       sprintf('原始渐近线 (斜率=%.2f)', order_orig), ...
       'Location', 'southeast');
grid on;

%% 8. 最终精度比较
fprintf('\n===== 最终精度比较 =====\n');
fprintf('原始最大 n = %d : 误差 = %.2e\n', n_list(end), error(end));
fprintf('Richardson 外推 (n=%d) : 误差 = %.2e\n', n_rich(end), error_rich(end));
fprintf('Wynn-ε 外推 : 误差 = %.2e\n', error_wynn);
fprintf('精确 π 值 = %.15f\n', pi_exact);

%% ========== 局部函数 ==========
function result = wynn_epsilon(seq)
    % Wynn-epsilon 算法用于序列加速
    seq = seq(:);
    n = length(seq);
    eps_table = nan(n, n);
    eps_table(:,1) = seq;
    % 第二列
    for i = 1:n-1
        diff = eps_table(i+1,1) - eps_table(i,1);
        if diff ~= 0
            eps_table(i,2) = 1 / diff;
        else
            eps_table(i,2) = Inf;
        end
    end
    % 后续列
    for j = 3:n
        for i = 1:n-j+1
            denom = eps_table(i+1,j-2) - eps_table(i,j-2) + ...
                     eps_table(i+1,j-1) - eps_table(i,j-1);
            if denom ~= 0 && ~isinf(denom)
                eps_table(i,j) = 1 / denom;
            else
                eps_table(i,j) = Inf;
            end
        end
    end
    % 提取偶数列对角线元素
    result = [];
    for k = 2:2:n
        if ~isinf(eps_table(1,k)) && ~isnan(eps_table(1,k))
            result = [result; eps_table(1,k)];
        end
    end
    if isempty(result)
        result = seq(end);
    else
        result = result(end);
    end
end