# Kramers-Wannier对偶缺陷的Fusion规则推导

## 目标

证明：$\mathcal{N} \times \mathcal{N} = 1 + \psi$

其中 $\mathcal{N}$ 是KW对偶缺陷，$\psi$ 是Majorana费米子。

---

## 1. 预备知识

### 1.1 Ising模型的算符内容

临界Ising模型（c = 1/2 CFT）有三个初级场：

| 场 | 共形维度 $h$ | 量子维度 $d$ | 物理意义 |
|----|-------------|-------------|----------|
| $1$ | 0 | 1 | 恒等算符 |
| $\psi$ | 1/2 | 1 | Majorana费米子 |
| $\sigma$ | 1/16 | $\sqrt{2}$ | 自旋场/序参量 |

### 1.2 KW对偶缺陷 $\mathcal{N}$

对偶缺陷 $\mathcal{N}$ 是一个**拓扑缺陷线**，实现KW对偶变换：
- 量子维度：$d_{\mathcal{N}} = \sqrt{2}$
- 作用：$\mathcal{N}: \sigma \leftrightarrow \mu$（序-无序对偶）

---

## 2. 格点模型推导

### 2.1 对偶变换的定义

在格点上，KW对偶由以下变换定义：

**原始变量 → 对偶变量**：
$$\tau^x_i = \sigma^z_i \sigma^z_{i+1}$$
$$\tau^z_i \tau^z_{i+1} = \sigma^x_{i+1}$$

### 2.2 对偶缺陷算符

在位置 $n$ 插入对偶缺陷，定义缺陷算符：

$$\mathcal{N}_n = \prod_{i \leq n} \sigma^x_i$$

这个算符将 $i \leq n$ 的区域变换到对偶表示。

### 2.3 两个缺陷的融合

考虑在位置 $m$ 和 $n$（$m < n$）插入两个缺陷：

$$\mathcal{N}_m \cdot \mathcal{N}_n = \left(\prod_{i \leq m} \sigma^x_i\right) \left(\prod_{j \leq n} \sigma^x_j\right)$$

由于 $(\sigma^x)^2 = 1$：

$$\mathcal{N}_m \cdot \mathcal{N}_n = \prod_{m < i \leq n} \sigma^x_i$$

### 2.4 识别融合产物

结果 $\prod_{m < i \leq n} \sigma^x_i$ 是什么？

**情况1**：如果 $n - m$ 为**偶数**
- 这是一个偶数个自旋翻转的乘积
- 属于**偶宇称扇区** → 对应 $1$（恒等）

**情况2**：如果 $n - m$ 为**奇数**
- 这是一个奇数个自旋翻转的乘积
- 属于**奇宇称扇区** → 对应 $\psi$（费米子）

因此：
$$\boxed{\mathcal{N} \times \mathcal{N} = 1 + \psi}$$

---

## 3. 弦算符推导

### 3.1 对偶缺陷作为弦算符

更精确地，对偶缺陷可以写成**半无限弦**：

$$\mathcal{N}(x) = \exp\left(i\frac{\pi}{4}\right) \prod_{i < x} \sigma^z_i \cdot \sigma(x)$$

其中 $\sigma(x)$ 是自旋场。

### 3.2 两个缺陷的乘积

$$\mathcal{N}(x_1) \mathcal{N}(x_2) = e^{i\pi/2} \prod_{i < x_1} \sigma^z_i \cdot \sigma(x_1) \cdot \prod_{j < x_2} \sigma^z_j \cdot \sigma(x_2)$$

当 $x_1 < x_2$ 时：

$$= e^{i\pi/2} \prod_{x_1 \leq i < x_2} \sigma^z_i \cdot \sigma(x_1) \sigma(x_2)$$

### 3.3 使用OPE

自旋场的OPE：
$$\sigma(x_1) \sigma(x_2) \sim \frac{1}{|x_1 - x_2|^{1/8}} \left(1 + |x_1 - x_2|^{1/2} \psi + \cdots\right)$$

乘上中间的弦 $\prod_{x_1 \leq i < x_2} \sigma^z_i$：
- 弦为空（$x_1 = x_2$）→ 得到 $1 + \psi$
- 弦非空 → 选择特定的宇称扇区

融合极限（$x_1 \to x_2$）：
$$\mathcal{N} \times \mathcal{N} = 1 + \psi$$

---

## 4. 量子维度验证

### 4.1 Fusion规则与量子维度

对于fusion规则 $a \times b = \sum_c N^c_{ab} \, c$，量子维度满足：

$$d_a \cdot d_b = \sum_c N^c_{ab} \, d_c$$

### 4.2 验证 $\mathcal{N} \times \mathcal{N} = 1 + \psi$

左边：
$$d_{\mathcal{N}} \times d_{\mathcal{N}} = \sqrt{2} \times \sqrt{2} = 2$$

右边：
$$d_1 + d_\psi = 1 + 1 = 2$$

✓ **一致！**

---

## 5. 范畴论视角

### 5.1 Ising范畴

Ising拓扑序由**Ising模范畴**描述，包含三个简单对象：$\{1, \psi, \sigma\}$

Fusion规则：
- $\psi \times \psi = 1$
- $\psi \times \sigma = \sigma$
- $\sigma \times \sigma = 1 + \psi$

### 5.2 对偶缺陷与自旋场

注意到 $\mathcal{N}$ 的fusion规则与 $\sigma$ 相同！

$$\sigma \times \sigma = 1 + \psi$$
$$\mathcal{N} \times \mathcal{N} = 1 + \psi$$

这不是巧合：**对偶缺陷 $\mathcal{N}$ 与自旋场 $\sigma$ 在拓扑层面等价**。

物理原因：
- $\sigma$ 在缺陷线 $\mathcal{N}$ 的端点
- $\mathcal{N}$ 可以看作 $\sigma$ 的"拓扑化"

### 5.3 F-矩阵

fusion空间的结合律由F-矩阵编码。对于 $\mathcal{N} \times \mathcal{N} \times \mathcal{N}$：

$$(\mathcal{N} \times \mathcal{N}) \times \mathcal{N} = (1 + \psi) \times \mathcal{N} = \mathcal{N} + \mathcal{N} = 2\mathcal{N}$$

$$\mathcal{N} \times (\mathcal{N} \times \mathcal{N}) = \mathcal{N} \times (1 + \psi) = \mathcal{N} + \mathcal{N} = 2\mathcal{N}$$

F-矩阵：
$$F^{\mathcal{N}\mathcal{N}\mathcal{N}}_{\mathcal{N}} = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$$

这是Hadamard矩阵，对应Ising任意子的非阿贝尔统计。

---

## 6. 物理解释

### 6.1 宇称与费米子

$\mathcal{N} \times \mathcal{N} = 1 + \psi$ 的物理意义：

两个对偶缺陷融合时，结果依赖于**费米宇称**：
- **偶宇称**：融合到 $1$（真空）
- **奇宇称**：融合到 $\psi$（费米子）

### 6.2 非阿贝尔性

由于有两个fusion通道，对偶缺陷是**非阿贝尔缺陷**：
- 多个缺陷的态空间有简并
- 编织操作给出非平凡的幺正变换
- 这是拓扑量子计算的基础

### 6.3 与Majorana零模的关系

在Majorana表示中：
- $\mathcal{N}$ 对应Majorana零模 $\gamma$
- $\gamma^2 = 1$ 但 $\{\gamma_1, \gamma_2\} = 0$
- 两个Majorana组成一个复费米子：$c = \gamma_1 + i\gamma_2$

Fusion结果：
- $|0\rangle$：无费米子 → $1$
- $c^\dagger|0\rangle$：有费米子 → $\psi$

---

## 7. 总结

### 推导要点

1. **格点构造**：$\mathcal{N}_m \mathcal{N}_n = \prod_{m<i\leq n} \sigma^x_i$，宇称决定fusion通道

2. **弦算符**：使用 $\sigma \times \sigma$ 的OPE加上弦贡献

3. **量子维度**：$\sqrt{2} \times \sqrt{2} = 1 + 1$ ✓

4. **范畴论**：$\mathcal{N}$ 与 $\sigma$ 有相同的fusion规则，同属Ising范畴

### 最终结果

$$\boxed{\mathcal{N} \times \mathcal{N} = 1 + \psi}$$

- $1$：偶宇称扇区（真空）
- $\psi$：奇宇称扇区（Majorana费米子）

这个fusion规则是Ising拓扑序和非阿贝尔任意子的核心特征。

---

## 参考文献

1. Oshikawa, M. & Affleck, I. (1997). *Phys. Rev. B* 55, 8323.
2. Petkova, V. B. & Zuber, J.-B. (2001). *Nucl. Phys. B* 603, 449.
3. Fröhlich, J. et al. (2010). *Commun. Math. Phys.* 294, 353.
4. Aasen, D., Mong, R. S. K. & Fendley, P. (2016). *J. Phys. A* 49, 354001.
