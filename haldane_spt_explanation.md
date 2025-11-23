# Haldane相为什么是SPT？

## 概述

Haldane相是**对称保护拓扑相（Symmetry Protected Topological Phase, SPT）**的典型范例。本文档将从物理本质上解释这一结论。

---

## 1. 什么是SPT相？

SPT相具有以下**核心特征**：

| 特征 | 描述 |
|------|------|
| **块体有能隙** | 基态与激发态之间有有限能量差 |
| **无局域序参量** | 不能用传统的对称性破缺来刻画 |
| **边界有受保护的低能模式** | 边界态由对称性保护，不能被局域微扰消除 |
| **绝热连续性** | 保持对称性时不能绝热变换到平庸相；破坏对称性后则可以 |

**关键**：SPT相的"拓扑性"完全依赖于对称性的存在。

---

## 2. Haldane相的基本设定

### 2.1 模型

考虑自旋-1反铁磁Heisenberg链：

$$H = J \sum_i \vec{S}_i \cdot \vec{S}_{i+1}, \quad J > 0$$

其中 $\vec{S}_i = (S^x_i, S^y_i, S^z_i)$ 是自旋-1算符。

### 2.2 Haldane猜想（1983）

Haldane预言：
- **整数自旋链**（S=1, 2, ...）：有能隙
- **半整数自旋链**（S=1/2, 3/2, ...）：无能隙（Lieb-Schultz-Mattis定理）

这个猜想已被数值和实验验证，其深刻原因与拓扑有关。

---

## 3. Haldane相的对称性

### 3.1 相关对称性

Haldane相受到以下对称性的保护：

1. **时间反演对称性** $\mathcal{T}$
   $$\mathcal{T}: \vec{S} \to -\vec{S}$$

2. **二面体群 $D_2 \cong \mathbb{Z}_2 \times \mathbb{Z}_2$**
   - $R_x = e^{i\pi \sum_i S^x_i}$：绕x轴旋转π
   - $R_y = e^{i\pi \sum_i S^y_i}$：绕y轴旋转π
   - $R_z = e^{i\pi \sum_i S^z_i}$：绕z轴旋转π

   满足 $R_x R_y = R_z$

3. **键反演（键中心反演）**

任一对称性都足以保护Haldane相的拓扑性质。

### 3.2 为什么需要对称性保护？

如果**破坏所有保护对称性**，Haldane相可以绝热地连接到平庸的直积态。例如：
- 加入单离子各向异性 $D(S^z)^2$ 并增大 $D \to \infty$
- 系统变成每个格点的 $|S^z = 0\rangle$ 直积态（平庸相）

但**保持对称性时**，这条路径被阻断——必须经过相变点。

---

## 4. 边界态与投影表示

### 4.1 自旋-1/2边界态

Haldane相最显著的特征是**边界有效自旋-1/2自由度**。

**物理图像（AKLT态）**：
- 将每个自旋-1分解为两个虚拟自旋-1/2
- 相邻格点的虚拟自旋形成单态（valence bond）
- 边界剩余一个未配对的自旋-1/2

```
●━━●  ●━━●  ●━━●  ●━━●
 ↑        ↑        ↑        ↑      → bulk pairs
↑                              ↑   → edge spin-1/2
```

### 4.2 投影表示（核心概念）

**这是理解SPT的关键！**

考虑 $D_2 = \mathbb{Z}_2 \times \mathbb{Z}_2$ 对称性：
- 块体上：$R_x R_y = R_z$（线性表示）
- **边界上**：$U_x U_y = -U_z$（投影表示）

边界自旋-1/2在 $R_\alpha$ 下变换为：
$$U_\alpha = e^{i\pi \sigma^\alpha/2} = i\sigma^\alpha$$

验证：
$$U_x U_y = (i\sigma^x)(i\sigma^y) = -\sigma^z \cdot i = -U_z$$

**投影表示不能在零维（边界）实现为线性表示！**

这意味着：
1. 边界必须有简并（基态）
2. 简并不能被保持对称性的微扰消除
3. 这就是"拓扑保护"的本质

---

## 5. 弦序参量

### 5.1 定义

Haldane相没有长程序，但有**非局域的弦序**：

$$\mathcal{O}_\alpha(i,j) = \langle S^\alpha_i \exp\left(i\pi \sum_{k=i+1}^{j-1} S^\alpha_k\right) S^\alpha_j \rangle$$

### 5.2 物理意义

- 在Haldane相中：$\lim_{|i-j|\to\infty} \mathcal{O}_\alpha \neq 0$
- 在平庸相中：$\mathcal{O}_\alpha \to 0$

弦序参量探测的是**隐藏的 $\mathbb{Z}_2 \times \mathbb{Z}_2$ 对称性破缺**。

**AKLT态的弦序**：可以精确计算为 $\mathcal{O}_z = 4/9$。

---

## 6. 拓扑纠缠熵

### 6.1 纠缠谱

对于有能隙的1D系统，纠缠谱的**简并结构**编码了拓扑信息。

在Haldane相中，纠缠谱有**2重简并**，对应于：
- 边界投影表示的维度
- $\gamma = \ln(2)$ 的拓扑纠缠熵

### 6.2 与 $\mathbb{Z}_2 \times \mathbb{Z}_2$ SPT的类比

本仓库实现的cluster模型：
$$H = -\sum_i Z_{i-1} X_i Z_{i+1}$$

与Haldane相有相同的 $\mathbb{Z}_2 \times \mathbb{Z}_2$ 保护结构：
- 相同的拓扑纠缠熵 $\gamma = \ln(2)$
- 相同的弦序结构
- 相同的边界投影表示

---

## 7. MPS/张量网络视角

### 7.1 矩阵乘积态分类

1D有能隙相可以用MPS有效描述：
$$|\psi\rangle = \sum_{\{s_i\}} \text{Tr}(A^{s_1} A^{s_2} \cdots A^{s_L}) |s_1 s_2 \cdots s_L\rangle$$

对称性 $g$ 在MPS上的作用：
$$\sum_{s'} U^g_{ss'} A^{s'} = e^{i\phi_g} V_g A^s V_g^{-1}$$

### 7.2 SPT分类

**关键定理**：不同的SPT相由 $V_g$ 的投影表示分类。

对于 $D_2 = \mathbb{Z}_2 \times \mathbb{Z}_2$：
- 线性表示 → 平庸相
- 投影表示 → Haldane相

投影表示由群的**二阶上同调群**分类：
$$H^2(G, U(1))$$

对于 $D_2$：$H^2(D_2, U(1)) = \mathbb{Z}_2$，所以有两个相（平庸 + Haldane）。

---

## 8. 总结：Haldane相是SPT的理由

| 性质 | Haldane相的表现 | SPT特征 |
|------|----------------|---------|
| **能隙** | 有限能隙 $\Delta \sim 0.41J$ | ✓ |
| **局域序** | 无 | ✓ |
| **边界态** | 自旋-1/2自由度 | ✓ |
| **对称性依赖** | 破坏 $D_2$ 或 $\mathcal{T}$ 后变平庸 | ✓ |
| **投影表示** | 边界为对称群的投影表示 | ✓ |
| **弦序** | $\mathcal{O}_\alpha \neq 0$ | ✓ |
| **拓扑熵** | $\gamma = \ln(2)$ | ✓ |

### 核心结论

**Haldane相是SPT因为：**

1. **对称性保护的边界简并**：边界自旋-1/2携带 $D_2$ 的投影表示，无法被局域微扰消除

2. **非平庸的上同调类**：属于 $H^2(\mathbb{Z}_2 \times \mathbb{Z}_2, U(1)) = \mathbb{Z}_2$ 的非平凡元素

3. **不可绝热变形**：保持对称性时，不能连续变换到直积态

---

## 9. 参考文献

1. Haldane, F. D. M. (1983). *Phys. Rev. Lett.* 50, 1153.
2. Affleck, Kennedy, Lieb, Tasaki (AKLT) (1987). *Phys. Rev. Lett.* 59, 799.
3. Gu, Z.-C. & Wen, X.-G. (2009). *Phys. Rev. B* 80, 155131.
4. Pollmann, Berg, Turner, Oshikawa (2012). *Phys. Rev. B* 85, 075125.
5. Chen, Gu, Wen (2011). *Phys. Rev. B* 83, 035107.

---

## 附录：代码验证

本仓库的 `z2z2_spt.py` 实现了 $\mathbb{Z}_2 \times \mathbb{Z}_2$ SPT的数值验证：

```python
# 运行计算
python z2z2_spt.py
```

输出将验证：
- 拓扑纠缠熵 $\gamma = \ln(2) \approx 0.693$
- 纠缠谱2重简并
- 非零弦序参量

这些结果与Haldane相具有相同的拓扑结构。
