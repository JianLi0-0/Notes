Active-Set 方法详解（用于求解二次规划 QP）

下面把 Active-Set 方法从直观思想、数学推导、算法流程、数值实现要点到工程实践技巧一步步讲清楚——要点都尽量配上公式与伪代码，方便你把理论和实现对上号。

⸻

一、直观概念（为什么要用 Active-Set）

Active-Set 方法把不等式约束分成两类：活跃约束（等式） 和 非活跃约束（不影响当前解的约束）。算法在每次迭代里假设一组约束是“活着”的（working set W），把它们当作等式约束去求解一个等式约束的 QP 子问题。根据子问题的解，动态加入或剔除约束，最后收敛到满足 KKT 条件的最优解。

直观上：Active-Set 就像人在爬山找最低点，每次锁住一组“墙”（活跃约束），沿着墙面找最优方向；如果能在该墙面继续下降就移动；如果不能再下降则检验墙的拉格朗日乘子，看是不是应该把某面墙拆掉。

⸻

二、问题与符号

我们考虑标准 QP：
\min_x\ \tfrac12 x^\top H x + g^\top x
\text{s.t. } \; a_i^\top x \le b_i \quad(i\in\mathcal{I}),\quad E x = d \quad(\text{等式})
边界（box）约束也可视为特殊的线性不等式。令工作集 W 为当前认为“活跃”的不等式集合（把它们当作等式 a_j^\top x = b_j 来处理），再把已有的真正等式一起构成 A_W x = b_W。

记当前点为 x^k，目标梯度 \nabla f(x^k)=H x^k+g（记作 q）。

⸻

三、等式约束子问题与 KKT 系统（数学核心）

在当前工作集 W 下，我们要求一个保持活约束可行且能最速下降的步长方向 p。子问题为（线性化）：
\min_p\ \tfrac12 p^\top H p + q^\top p\quad\text{s.t. } A_W p = 0
（因为若要保持 A_W(x^k+p)=b_W，需要 A_W p = 0。）

对子问题写 KKT 条件（未知量是 p 与对应乘子 \lambda）：
\begin{bmatrix} H & A_W^\top \\ A_W & 0 \end{bmatrix}
\begin{bmatrix} p \\ \lambda \end{bmatrix}

\begin{bmatrix} -q \\ 0 \end{bmatrix}.
	•	若得到 p\neq 0，说明存在沿着满足活约束的方向下降的方向。
	•	若 p=0，说明在当前工作集上无法继续降低目标，此时检查乘子 \lambda 的符号来判定最优性（对形式 a^\top x \le b 的约束，最优要求对应乘子 \lambda\ge0；若存在 \lambda_j<0，说明该约束在最优性上不应该被锁定，应当从 W 中移除）。

⸻

四、一步的移动与阻塞约束（step length）

当 p\neq0 时，需计算最大可行步长 \alpha\in(0,1] 使得对所有未被激活的约束仍成立：
对于某个未激活约束 a_i^\top x \le b_i（当前残量 r_i=b_i-a_i^\top x^k>0），如果 a_i^\top p>0（前进会增加 a_i^\top x），则可能会被“撞上”：
\alpha_i = \frac{r_i}{a_i^\top p}\quad(\text{当 } a_i^\top p>0).
取最小的正 \alpha_i（及1）作为步长：\alpha=\min(1,\min_{a_i^\top p>0}\alpha_i)。如果 \alpha<1，说明某个不活跃约束变为活跃，将其加入 W。

⸻

五、完整伪代码（Primal Active-Set）

给定可行起点 x, 构造初始工作集 W（包含所有当前活跃约束）
repeat:
  1. 计算 q = H x + g
  2. 解 KKT 系统: [H A_W^T; A_W 0] [p; λ] = [-q; 0]
  3. if ||p|| > tol:
       # 有下降方向，沿 p 前进到阻塞约束或完全前进
       计算 α = min(1, min_i {(b_i - a_i^T x)/(a_i^T p)} over i with a_i^T p > 0)
       x ← x + α p
       if α < 1:
         将对应阻塞约束加入 W
       continue
     else:
       # p ≈ 0：在当前工作集上已无下降方向
       检查 λ 的符号（对应不等式约束应满足 λ >= 0）
       if 所有 λ 满足可行性:
         return x 为最优解
       else:
         从 W 中移除一个违反符号条件的约束（常选最负的 λ）
         continue
end


⸻

六、数值实现要点（如何高效解上面 KKT）

直接解 KKT 的线性系统是关键且代价最高的一步。常见方案：
	1.	Null-space method（零空间法）
找到 Z 使 A_W Z = 0（列向量构成），将 p = Z p_z。
子问题化为较小维度的线性系统：
(Z^\top H Z) p_z = -Z^\top q.
这个方法避免了求解非定的 KKT 矩阵，适合当活跃约束数目小于变量数时。
	2.	Range-space / Schur complement（伴随空间法）
假设 H 可逆，先解 H 的系统，推导对 \lambda 的方程：
-A_W H^{-1} A_W^\top \lambda = A_W H^{-1} q.
求出 \lambda 后再回代得到 p。当变量维数小、或 H^{-1} 易得时适用。
	3.	直接 KKT 因式分解
对 KKT 矩阵做因式分解（LDLᵀ / 变形 Cholesky），适合不频繁改变工作集或可用高效更新技巧。
	4.	因式分解的增量更新
Active-Set 每次只加入或删除一个约束，利用 rank-1 更新/降级（update/downdate）来更新已有分解，远比每次重建分解高效。qpOASES 就是大量使用这种增量更新以支持 warm start。

⸻

七、最优性判据（KKT 条件回顾）

对于原问题，最优需满足 KKT：
	•	可行性（原始约束满足）
	•	驻点条件： H x^* + g + A^\top \lambda^* = 0
	•	对不等式约束的乘子满足：若 a_i^\top x^* < b_i 则 \lambda_i^=0（互补松弛），若 a_i^\top x^ = b_i 则 \lambda_i^*\ge0（当约束形式是 a_i^\top x \le b_i）。

Active-Set 的检验方式就是用子问题得到的 \lambda 来确认工作集是否能维持最优性。

⸻

八、Phase-I（如何得到初始可行点）

Active-Set 需要一个可行起点。若没有，需要做 Phase-I：
	•	构造人工变量 s\ge0 来允许约束违反，目标尽量减小 s（例如最小化 \sum s_i）。
	•	找到使 s=0 的可行解或判断原问题不可行。
很多实现把 Phase-I 也作为 QP 解决。

⸻

九、数值鲁棒性与工程实践建议
	•	H 的正定性：若 H 严格正定（convex QP），Active-Set 会在有限步内找最优解且唯一。若 H 半正定或退化，需小的正则项（add \epsilon I）或专门处理。
	•	尺度与预处理：对变量/约束进行尺度化可以显著提升数值稳定性。
	•	容差设定：设置合适的可行性、乘子和步长容差避免频繁添加/删除约束导致震荡。
	•	避免循环（cycling）：使用反循环策略或微扰/基于乘子规则避免无限往返。
	•	warm-start：Active-Set 对 warm-start 友好——当 QP 参数随时间小幅变化（MPC 场景），用上次解和工作集初始化可极大加速。
	•	把 bounds 单独处理：实现中通常把箱约束（lb\le x\le ub）作为特殊、快速检查和更新的部分（qpOASES 也是这样），因为处理更简单且更频繁。

⸻

十、优缺点比较（与 Interior-Point / Gradient 法）

优点：
	•	对 warm-start 非常好（在线 MPC 最适合）。
	•	若最终活跃约束很少，迭代次数通常少，性能优。
	•	可直接得到稀疏结构的 KKT（便于解释）。

缺点：
	•	每次迭代需要解线性系统且可能需要多次加入/删除约束，若约束很多或活跃集频繁变化可能开销大。
	•	对高度约束或密集约束的大规模问题，Interior-Point 方法常更稳定和更快（批量处理大量约束）。
	•	可能遇到退化/循环问题，需要工程策略处理。

⸻

十一、常见变体与实现细节
	•	Primal Active-Set（上文所述）：在原始变量空间搜索下降方向。
	•	Dual Active-Set：在对偶空间操作（当变量数远大于约束数，或反过来，也有选择）。
	•	在线 Active-Set（如 qpOASES）：做高效的分解更新、支持快速 warm-start、适合 MPC。
	•	组合策略：先用 Active-Set warm-start，再用 IPM 微调——在某些场景可取两者优点。

⸻

十二、一个简单的二次例子（帮助理解步骤）

问题：最小化 \tfrac12(x_1^2+x_2^2)，约束 x_1+x_2 \ge 1、x_1,x_2\ge0。
	•	初始可行点取 x=(1,0)（满足所有约束），工作集 W 包含 x_2\ge0（是否活跃取决于点）与等式 x_1+x_2\ge1 是否活跃（此点上等式成立所以活跃）。
按上面的子问题求 p，若 p=0 且乘子都满足非负，则最优；否则沿 p 前进并根据阻塞约束调整 W。（这里目标在可行域边界上通常会移动到 (0.5,0.5) 等）

⸻

十三、结论与实践建议（总结）
	•	Active-Set 方法通过在每次迭代上解一个等式约束 QP 子问题，按需增删约束来逼近最优解。
	•	它非常适合问题序列、需要 warm-start 的在线场景（如 MPC）。
	•	实现上要注意 KKT 系统的高效求解（null-space、Schur 或 LDLᵀ）、分解的增量更新、Phase-I 可行性处理、数值正则化与容差策略。
	•	在工程中常把 box bounds 单独优化处理以提高速度（qpOASES 的做法）。

⸻

如果你愿意，我可以：
	•	写一段带注释的 伪代码 + 矩阵运算步骤（便于编码实现），或者
	•	给出 简单 2D 数值例子（逐步计算每次 p,\lambda,\alpha 的值）来演示算法走向，或者
	•	针对你使用的 qpOASES，解释它内部是怎样做增量分解与 warm-start 的（如果你想我可以结合 qpOASES 的接口细节来说明）。

你想要哪一种后续材料？