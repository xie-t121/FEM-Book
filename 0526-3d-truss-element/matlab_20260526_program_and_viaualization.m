%% ========================================================================
% 程序作业：三维杆单元刚度矩阵与应力计算
% 学号：20252110002_姓名：解云涛_三维杆单元程序设计作业
% 语言：MATLAB R2023b (R2020a及以上)
% =========================================================================

%% 任务1：公式推导（本节注释中给出完整推导过程）
% 
% 设三维杆单元节点1坐标 (x1,y1,z1)，节点2坐标 (x2,y2,z2)
% 1. 单元长度 L = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)
% 2. 方向余弦：
%    cx = (x2-x1)/L,  cy = (y2-y1)/L,  cz = (z2-z1)/L
% 3. 单元轴向伸长量 Δ = ( -cx, -cy, -cz, cx, cy, cz ) * de
%    其中 de = [u1, v1, w1, u2, v2, w2]^T
% 4. 应变 ε = Δ / L = (1/L) * ( -cx, -cy, -cz, cx, cy, cz ) * de
%    应变-位移矩阵 B = (1/L) * [ -cx, -cy, -cz, cx, cy, cz ]   (1×6)
% 5. 单元刚度矩阵（全局坐标系）：
%    Ke = E*A/L * [  cx^2,    cx*cy,   cx*cz,  -cx^2,  -cx*cy,  -cx*cz;
%                   cx*cy,    cy^2,    cy*cz,  -cx*cy, -cy^2,   -cy*cz;
%                   cx*cz,    cy*cz,   cz^2,   -cx*cz, -cy*cz,  -cz^2;
%                  -cx^2,   -cx*cy,  -cx*cz,   cx^2,    cx*cy,   cx*cz;
%                 -cx*cy,   -cy^2,   -cy*cz,   cx*cy,   cy^2,    cy*cz;
%                 -cx*cz,   -cy*cz,  -cz^2,    cx*cz,   cy*cz,   cz^2 ]
% 6. 应力 σ = E * ε，轴力 N = σ * A
% =========================================================================

%% ====================== 主程序：验证算例 ===============================
fprintf('========== 三维杆单元程序验证 ==========\n\n');

%% 算例1：沿x轴的一维杆单元（任务2 & 任务3的一部分）
fprintf('--- 算例1：沿x轴杆单元 ---\n');
x1 = [0,0,0]; x2 = [2,0,0];
E = 200e9;      % 200 GPa
A = 1.0e-4;     % 0.0001 m^2
de = [0,0,0, 1.0e-3, 0, 0]';   % 节点2 x方向位移1mm

% 调用任务2实现的函数
[L, dir_cos, Ke] = truss3d_element_stiffness(x1, x2, E, A);
[epsilon, sigma, N] = truss3d_element_stress(x1, x2, E, A, de);

% 输出结果
fprintf('单元长度 L = %.6f m (理论 2 m)\n', L);
fprintf('方向余弦 (cx,cy,cz) = (%.1f, %.1f, %.1f)\n', dir_cos(1), dir_cos(2), dir_cos(3));
fprintf('刚度矩阵 Ke (6x6):\n');
disp(Ke);
fprintf('轴向应变 epsilon = %.4e (理论 5e-4)\n', epsilon);
fprintf('轴向应力 sigma = %.2f MPa (理论 100 MPa)\n', sigma/1e6);
fprintf('轴力 N = %.2f N (理论 10000 N)\n\n', N);

% 验证
assert(abs(L-2)<1e-6 && norm(dir_cos-[1,0,0])<1e-6, '算例1失败');
fprintf('算例1验证通过 ✓\n\n');

%% 算例2：空间任意方向杆单元（任务2 & 任务3的一部分）
fprintf('--- 算例2：空间杆单元 (1,2,2)方向 ---\n');
x1 = [0,0,0]; x2 = [1,2,2];
E = 210e9;      % 210 GPa
A = 2.0e-4;     % 0.0002 m^2
de = [0,0,0, 1.0e-3, 2.0e-3, 2.0e-3]';   % 节点2位移沿(1,2,2)

[L, dir_cos, Ke] = truss3d_element_stiffness(x1, x2, E, A);
[epsilon, sigma, N] = truss3d_element_stress(x1, x2, E, A, de);

fprintf('单元长度 L = %.6f m (理论 3 m)\n', L);
fprintf('方向余弦 (cx,cy,cz) = (%.4f, %.4f, %.4f)\n', dir_cos(1), dir_cos(2), dir_cos(3));
fprintf('刚度矩阵 Ke (6x6):\n');
disp(Ke);
fprintf('轴向应变 epsilon = %.4e (理论 1e-3)\n', epsilon);
fprintf('轴向应力 sigma = %.2f MPa (理论 210 MPa)\n', sigma/1e6);
fprintf('轴力 N = %.2f N (理论 42000 N)\n\n', N);

assert(abs(L-3)<1e-6 && norm(dir_cos-[1/3,2/3,2/3])<1e-6, '算例2长度/方向余弦错');
assert(abs(epsilon-1e-3)<1e-8 && abs(sigma-210e6)<1e4 && abs(N-42000)<1, '算例2应力/轴力错');
fprintf('算例2验证通过 ✓\n\n');

%% 任务3：单元刚度矩阵性质验证（使用算例2的Ke）
fprintf('========== 任务3：刚度矩阵性质验证 ==========\n');

% 3.1 对称性检查
sym_error = norm(Ke - Ke', 'fro');
fprintf('1. 对称性检查：||Ke - Ke^T||_F = %.2e (应为0)\n', sym_error);
if sym_error < 1e-10
    fprintf('   → 刚度矩阵对称 ✓\n');
end

% 3.2 奇异性检查（行列式）
det_Ke = det(Ke);
fprintf('2. 行列式 det(Ke) = %.2e (应接近0，奇异)\n', det_Ke);
if abs(det_Ke) < 1e-6
    fprintf('   → 刚度矩阵奇异 ✓\n');
end

% 3.3 半正定性检查（特征值）
eigvals = eig(Ke);
fprintf('3. 特征值：\n');
for i = 1:6
    fprintf('   λ%d = %.4e\n', i, eigvals(i));
end
if all(eigvals >= -1e-10)
    fprintf('   → 所有特征值非负（半正定）✓\n');
end

% 3.4 刚体平移测试（给刚体位移，应产生零内力）
rigid_de = [0.1, 0.2, 0.3, 0.1, 0.2, 0.3]';   % 刚体平移
Fe_rigid = Ke * rigid_de;
fprintf('4. 刚体平移位移 de = [0.1,0.2,0.3,0.1,0.2,0.3]^T\n');
fprintf('   产生的节点内力 Fe = [%.2e, %.2e, %.2e, %.2e, %.2e, %.2e]^T\n', Fe_rigid);
if norm(Fe_rigid) < 1e-10
    fprintf('   → 刚体平移不产生内力，满足平衡 ✓\n');
end
fprintf('\n');

%% 任务4：刚度矩阵物理意义验证（第j列的含义）
fprintf('========== 任务4：刚度矩阵物理意义验证 ==========\n');
% 选第4个自由度（节点2的x方向位移） j = 4
j = 4;
de_j = zeros(6,1);
de_j(j) = 1.0;   % 只有第j个自由度为1，其余为0
Fe_j = Ke * de_j;

fprintf('取 j = %d （节点2的x方向位移），令该自由度位移=1，其余为0\n', j);
fprintf('产生的节点力列阵 Fe = Ke * e_%d = [', j);
fprintf('%.2f, ', Fe_j(1:5));
fprintf('%.2f]^T\n', Fe_j(6));

fprintf('\n解释：\n');
fprintf('- 该节点力列阵正好等于刚度矩阵的第 %d 列。\n', j);
fprintf('- 物理意义：要维持单元在第 %d 个自由度产生单位位移而其他自由度固定，\n', j);
fprintf('  需要在每个自由度上施加的节点力。\n');
fprintf('- 即 Ke(i,%d) 表示：当第 %d 个自由度有单位位移时，在第 i 个自由度上产生的节点力。\n', j, j);
fprintf('- 例如 Ke(1,%d)=%.2f 表示此时在节点1 x方向产生的力。\n', j, Ke(1,j));
fprintf('- 这验证了课件中所说：“当单元的第j个自由度给定位移而其它自由度固定时，\n');
fprintf('  为了使单元平衡而需要在各结点自由度上施加的结点力”。\n');

fprintf('\n========== 程序运行完毕 ==========\n');

%% ========================= 可视化板块 ==================================
% 不删除原有内容，在函数定义之前添加绘图
% 每个任务对应的可视化展示

% ---------- 任务1：公式推导可视化窗口（修复版：纯文本，无LaTeX）----------
figure('Name', '任务1：三维杆单元公式推导过程', 'NumberTitle', 'off', ...
       'Position', [100 100 900 700], 'Color', 'white');
axis off;

% 标题
text(0.05, 0.95, '三维杆单元刚度矩阵与应力计算公式推导', ...
     'FontSize', 14, 'FontWeight', 'bold');

% 1. 长度
text(0.05, 0.88, '1. 节点坐标与单元长度', 'FontSize', 12, 'FontWeight', 'bold');
text(0.05, 0.83, '   节点1: (x1, y1, z1)    节点2: (x2, y2, z2)', 'FontSize', 11);
text(0.05, 0.78, '   长度: L = sqrt( (x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2 )', 'FontSize', 11);

% 2. 方向余弦
text(0.05, 0.71, '2. 方向余弦', 'FontSize', 12, 'FontWeight', 'bold');
text(0.05, 0.66, '   cx = (x2-x1)/L,   cy = (y2-y1)/L,   cz = (z2-z1)/L', 'FontSize', 11);

% 3. 应变
text(0.05, 0.59, '3. 轴向伸长量与应变', 'FontSize', 12, 'FontWeight', 'bold');
text(0.05, 0.54, '   轴向伸长量: Δ = [-cx, -cy, -cz, cx, cy, cz] * d^e', 'FontSize', 11);
text(0.05, 0.49, '   其中 d^e = [u1, v1, w1, u2, v2, w2]^T', 'FontSize', 11);
text(0.05, 0.44, '   应变: ε = Δ / L = (1/L) * [-cx, -cy, -cz, cx, cy, cz] * d^e', 'FontSize', 11);

% 4. 刚度矩阵（分行写）
text(0.05, 0.37, '4. 单元刚度矩阵 (全局坐标系)', 'FontSize', 12, 'FontWeight', 'bold');
text(0.05, 0.32, '   Ke = (E*A/L) * [', 'FontSize', 11);
text(0.05, 0.27, '        cx^2,   cx*cy,   cx*cz,  -cx^2,  -cx*cy,  -cx*cz;', 'FontSize', 10);
text(0.05, 0.23, '        cx*cy,   cy^2,    cy*cz,  -cx*cy, -cy^2,   -cy*cz;', 'FontSize', 10);
text(0.05, 0.19, '        cx*cz,   cy*cz,   cz^2,   -cx*cz, -cy*cz,  -cz^2;', 'FontSize', 10);
text(0.05, 0.15, '       -cx^2,  -cx*cy,  -cx*cz,   cx^2,   cx*cy,   cx*cz;', 'FontSize', 10);
text(0.05, 0.11, '       -cx*cy, -cy^2,   -cy*cz,   cx*cy,  cy^2,    cy*cz;', 'FontSize', 10);
text(0.05, 0.07, '       -cx*cz, -cy*cz,  -cz^2,    cx*cz,  cy*cz,   cz^2 ]', 'FontSize', 10);

% 5. 应力轴力
text(0.55, 0.88, '5. 应力与轴力', 'FontSize', 12, 'FontWeight', 'bold');
text(0.55, 0.83, '   应力: σ = E * ε', 'FontSize', 11);
text(0.55, 0.78, '   轴力: N = σ * A = E*A*ε', 'FontSize', 11);

% 6. 性质
text(0.55, 0.71, '6. 单元刚度矩阵性质', 'FontSize', 12, 'FontWeight', 'bold');
text(0.55, 0.66, '   • 对称性: Ke = Ke^T', 'FontSize', 11);
text(0.55, 0.61, '   • 奇异性: det(Ke) = 0 (存在刚体位移)', 'FontSize', 11);
text(0.55, 0.56, '   • 半正定性: 特征值 λ_i ≥ 0', 'FontSize', 11);
text(0.55, 0.51, '   • 每行(列)元素之和为0', 'FontSize', 11);

% 7. 物理意义
text(0.55, 0.44, '7. 物理意义示例', 'FontSize', 12, 'FontWeight', 'bold');
text(0.55, 0.39, '   第 j 列: 当 d_j = 1 (其余为0) 时', 'FontSize', 11);
text(0.55, 0.34, '   节点力列阵 F^e = Ke(:,j)', 'FontSize', 11);

title('三维杆单元公式推导总览');

% ---------- 原有可视化内容（保持不变）----------
% 创建图形窗口（不覆盖命令行输出）
figure('Name', '三维杆单元任务可视化', 'Position', [50 50 1200 800]);

%% 子图1：算例1与算例2的几何变形（任务1/2）
subplot(2,3,1);
% 算例1：一维杆单元
x1_1 = [0,0,0]; x2_1 = [2,0,0];
de_1 = [0,0,0, 1e-3,0,0]';
[L1, dc1, ~] = truss3d_element_stiffness(x1_1, x2_1, E, A);
x1_def1 = x1_1 + 10*[de_1(1), de_1(2), de_1(3)];  % 放大10倍
x2_def1 = x2_1 + 10*[de_1(4), de_1(5), de_1(6)];
plot3([x1_1(1), x2_1(1)], [x1_1(2), x2_1(2)], [x1_1(3), x2_1(3)], 'b-o', 'LineWidth', 2); hold on;
plot3([x1_def1(1), x2_def1(1)], [x1_def1(2), x2_def1(2)], [x1_def1(3), x2_def1(3)], 'r--o', 'LineWidth', 1.5);
text(x1_1(1), x1_1(2), x1_1(3), ' 1', 'Color', 'b');
text(x2_1(1), x2_1(2), x2_1(3), ' 2', 'Color', 'b');
text(x1_def1(1), x1_def1(2), x1_def1(3), ' 1''', 'Color', 'r');
text(x2_def1(1), x2_def1(2), x2_def1(3), ' 2''', 'Color', 'r');
xlabel('X'); ylabel('Y'); zlabel('Z'); title('算例1：一维杆变形 (放大10倍)');
grid on; axis equal; view(30,30); legend('原始','变形','Location','best');

subplot(2,3,2);
% 算例2：空间杆单元
x1_2 = [0,0,0]; x2_2 = [1,2,2];
de_2 = [0,0,0, 1e-3,2e-3,2e-3]';
[L2, dc2, ~] = truss3d_element_stiffness(x1_2, x2_2, E, A);
scale = 10;
x1_def2 = x1_2 + scale*[de_2(1), de_2(2), de_2(3)];
x2_def2 = x2_2 + scale*[de_2(4), de_2(5), de_2(6)];
plot3([x1_2(1), x2_2(1)], [x1_2(2), x2_2(2)], [x1_2(3), x2_2(3)], 'b-o', 'LineWidth', 2); hold on;
plot3([x1_def2(1), x2_def2(1)], [x1_def2(2), x2_def2(2)], [x1_def2(3), x2_def2(3)], 'r--o', 'LineWidth', 1.5);
% 绘制方向余弦向量（从原点出发）
quiver3(0,0,0, dc2(1)*L2, dc2(2)*L2, dc2(3)*L2, 'g', 'LineWidth', 1.5, 'MaxHeadSize', 0.2);
text(x1_2(1), x1_2(2), x1_2(3), ' 1', 'Color', 'b');
text(x2_2(1), x2_2(2), x2_2(3), ' 2', 'Color', 'b');
text(x1_def2(1), x1_def2(2), x1_def2(3), ' 1''', 'Color', 'r');
text(x2_def2(1), x2_def2(2), x2_def2(3), ' 2''', 'Color', 'r');
text(dc2(1)*L2/2, dc2(2)*L2/2, dc2(3)*L2/2, ' 方向余弦', 'Color', 'g');
xlabel('X'); ylabel('Y'); zlabel('Z'); title(sprintf('算例2：空间杆变形 (放大%d倍)', scale));
grid on; axis equal; view(135,30); legend('原始','变形','方向','Location','best');
hold off;

%% 子图3：任务3 - 刚度矩阵热图与特征值分布
subplot(2,3,3);
% 使用算例2的Ke
imagesc(Ke); colorbar; colormap(jet);
title('刚度矩阵 Ke 热图 (6×6)');
xlabel('列号'); ylabel('行号');
% 在每个单元格显示数值（保留两位有效数字）
for i = 1:6
    for j = 1:6
        text(j, i, sprintf('%.1e', Ke(i,j)), ...
            'HorizontalAlignment', 'center', 'Color', 'w', 'FontSize', 8);
    end
end
axis equal tight;

subplot(2,3,4);
% 特征值柱状图
bar(eigvals, 'FaceColor', [0.2 0.6 0.8]);
xlabel('特征值序号'); ylabel('特征值 (N/m)');
title('刚度矩阵特征值 (半正定性)');
grid on;
% 标注零特征值（刚体模态）
hold on;
for i = 1:length(eigvals)
    if abs(eigvals(i)) < 1e-6
        plot(i, eigvals(i), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    end
end
legend('特征值', '零特征值 (刚体模态)');
hold off;

%% 子图4：任务4 - 单位位移对应的节点力向量
subplot(2,3,5);
% j = 4 时的节点力
bar(Fe_j, 'FaceColor', [0.8 0.2 0.2]);
xlabel('自由度编号 (1~6)'); ylabel('节点力 (N)');
title(sprintf('单位位移 d_%d = 1 时的节点力 (Ke的第%d列)', j, j));
grid on;
% 添加数值标签
for i = 1:6
    text(i, Fe_j(i), sprintf('%.1f', Fe_j(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

%% 子图5：刚体位移验证（任务3的补充）
subplot(2,3,6);
% 刚体平移位移的节点力向量
bar(Fe_rigid, 'FaceColor', [0.2 0.8 0.2]);
xlabel('自由度编号'); ylabel('节点力 (N)');
title('刚体平移位移 (0.1,0.2,0.3) 产生的节点力 (应为零)');
grid on;
ylim([-1e-10 1e-10]);  % 由于数值误差，范围极小
% 标注接近零
for i = 1:6
    if abs(Fe_rigid(i)) < 1e-12
        text(i, Fe_rigid(i), '≈0', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end

sgtitle('三维杆单元作业可视化总览');

% 额外单独窗口：任务3的刚度矩阵热图放大（可选）
figure('Name', '刚度矩阵详细视图');
imagesc(Ke); colorbar; colormap(jet);
title('刚度矩阵 Ke (6×6) - 颜色表示数值大小');
xlabel('列 j'); ylabel('行 i');
for i = 1:6
    for j = 1:6
        text(j, i, sprintf('%.2e', Ke(i,j)), ...
            'HorizontalAlignment', 'center', 'Color', 'k', 'FontSize', 9);
    end
end
axis equal tight;

fprintf('\n可视化绘图已完成。请查看弹出的图形窗口。\n');

%% ====================== 任务2：函数实现 =================================
% 函数1：计算单元长度、方向余弦和全局刚度矩阵
function [L, dir_cos, Ke] = truss3d_element_stiffness(x1, x2, E, A)
    % 输入：x1, x2 - 1×3 坐标向量 [x,y,z]
    %       E    - 弹性模量 (Pa)
    %       A    - 截面积 (m^2)
    % 输出：L       - 单元长度 (m)
    %       dir_cos - 1×3 方向余弦 [cx, cy, cz]
    %       Ke      - 6×6 全局刚度矩阵 (N/m)
    
    dx = x2(1) - x1(1);
    dy = x2(2) - x1(2);
    dz = x2(3) - x1(3);
    L = sqrt(dx^2 + dy^2 + dz^2);
    
    % 退化单元检查
    if L < 1e-12
        error('错误：两个节点重合，单元长度为零，无法继续计算。');
    end
    
    cx = dx / L;  cy = dy / L;  cz = dz / L;
    dir_cos = [cx, cy, cz];
    
    k = E * A / L;
    c2x = cx^2;   c2y = cy^2;   c2z = cz^2;
    cxy = cx*cy;  cxz = cx*cz;  cyz = cy*cz;
    
    Ke = k * [ c2x,  cxy,  cxz, -c2x, -cxy, -cxz;
               cxy,  c2y,  cyz, -cxy, -c2y, -cyz;
               cxz,  cyz,  c2z, -cxz, -cyz, -c2z;
              -c2x, -cxy, -cxz,  c2x,  cxy,  cxz;
              -cxy, -c2y, -cyz,  cxy,  c2y,  cyz;
              -cxz, -cyz, -c2z,  cxz,  cyz,  c2z ];
end

% 函数2：根据节点位移计算应变、应力和轴力
function [epsilon, sigma, N] = truss3d_element_stress(x1, x2, E, A, de)
    % 输入：de - 6×1 节点位移列阵 [u1;v1;w1;u2;v2;w2] (m)
    % 输出：epsilon - 轴向应变 (无量纲)
    %       sigma   - 轴向应力 (Pa)
    %       N       - 轴力 (N)，受拉为正
    
    [L, dir_cos, ~] = truss3d_element_stiffness(x1, x2, E, A);
    cx = dir_cos(1); cy = dir_cos(2); cz = dir_cos(3);
    
    delta = (-cx)*de(1) + (-cy)*de(2) + (-cz)*de(3) + ...
             cx*de(4)  +  cy*de(5)  +  cz*de(6);
    
    epsilon = delta / L;
    sigma = E * epsilon;
    N = sigma * A;
end