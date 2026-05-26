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