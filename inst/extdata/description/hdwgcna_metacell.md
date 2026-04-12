
考虑到 metacell 的构建依赖于基于低维空间的 $k$ 近邻图（$k$-nearest neighbors, kNN），参数选择需满足图结构约束，以避免过度重叠或无法构建有效聚合单元。具体而言，我们基于各细胞亚群中最小细胞数 $n_{\min}$ 自适应地确定关键参数。首先，邻域大小 $k$ 设置为：

$$
k \approx \sqrt{n_{\min}}
$$

以在局部平滑与结构分辨率之间取得平衡。同时，为避免在极小或极大样本规模下产生不稳定结构，我们对 $k$ 施加上下界约束：

$$
k = \min\left(\max\left(\lfloor \sqrt{n_{\min}} \rfloor, k_{\min}\right), \left\lfloor \frac{n_{\min}}{2} \right\rfloor, k_{\max}\right)
$$

其中 $k_{\min} = ⟦k_min⟧$，$k_{\max} = ⟦k_max⟧$，并进一步限制 $k \le \left\lfloor \frac{n_{\min}}{2} \right\rfloor$，以防止邻域过度重叠。

其次，通过约束每个细胞可参与的最大共享次数（$max_{shared}$），保证在给定细胞数条件下能够生成足够数量的 metacell。该参数通过满足以下条件进行估计：

$$
\frac{n_{\min} \cdot max_{shared}}{k} \ge M_{\text{target}}
$$

其中 $M_{\text{target}}$ 表示期望构建的 metacell 数量 (设置为 ⟦target_metacells⟧)。

此外，最小细胞数阈值（$min_{cells}$）设置为：

$$
min_{cells} = \lfloor 0.9 \cdot n_{\min} \rfloor
$$

用于排除极小亚群以提高网络稳定性。该参数选择策略从图论约束出发，在保证 metacell 可构建性的同时，引入参数边界以增强方法在不同数据规模下的稳健性，从而最大程度保留生物学变异信息，提高后续加权基因共表达网络分析的鲁棒性与可解释性。


