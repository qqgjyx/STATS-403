# Problem 5. (30pt)

The Maximum Mean Discrepancy (MMD) between two domain samples $X_s, X_t \subset \mathbb{R}^D$ is defined to be

$$\text{MMD}(X_s, X_t) = \left\| \frac{1}{n_s} \sum_{i=1}^{n_s} \phi(x^s_i) - \frac{1}{n_t} \sum_{i=1}^{n_t} \phi(x^t_i) \right\|_{\mathcal{H}},$$

where $\mathcal{H}$ is a real-valued reproducing kernel Hilbert space (RKHS) with a kernel $k$ such that $k(x, x^{\prime}) = \langle \phi(x), \phi(x^{\prime}) \rangle_{\mathcal{H}}$.

Define the $(n_s + n_t) \times (n_s + n_t)$ matrix $K$ to be:

$$K = \begin{bmatrix}
K_{s,s} & K_{s,t} \\
K^\top_{s,t} & K_{t,t}
\end{bmatrix}$$

in which the $(i, j)$-th entry of $K_{u,v}$ is given by $k(x^u_i, x^v_j)$, where $u, v \in \{s, t\}$.

Also define a matrix $L$ of the same shape:

$$L = \begin{bmatrix}
L_{s,s} & L_{s,t} \\
L^\top_{s,t} & L_{t,t}
\end{bmatrix}$$

where:
- $L_{s,s}$: all entries $\frac{1}{n_s^2}$
- $L_{t,t}$: all entries $\frac{1}{n_t^2}$
- $L_{s,t}$: all entries $-\frac{1}{n_s n_t}$

1. Prove: $\text{MMD}^2(X_s, X_t) = \text{tr}(KL)$
2. In general, for two distributions $P$ and $Q$, define $\text{MMD}(P, Q) = \| \mathbb{E}_{x \sim P}[\phi(x)] - \mathbb{E}_{x \sim Q}[\phi(x)] \|_{\mathcal{H}}$.
   - (a) Prove: For any $f \in \mathcal{H}$, $\|f\|_{\mathcal{H}} = \sup_{\|g\|_{\mathcal{H}} \leq 1} \langle g, f \rangle_{\mathcal{H}}$.
   - (b) Prove:
$$\text{MMD}(P, Q) = \sup_{\|f\|_{\mathcal{H}} \leq 1} \mathbb{E}_{x \sim P}[f(x)] - \mathbb{E}_{x \sim Q}[f(x)]$$

---

Let

$$\mu_s=\frac{1}{n_s}\sum_{i=1}^{n_s}\phi(x^s_i),\qquad
\mu_t=\frac{1}{n_t}\sum_{j=1}^{n_t}\phi(x^t_j).$$

Then

$$\operatorname{MMD}^2(X_s,X_t)=\|\mu_s-\mu_t\|_{\mathcal{H}}^2
=\langle\mu_s-\mu_t,\mu_s-\mu_t\rangle_{\mathcal{H}}.$$

Step 1: Expand the inner product

$$
\begin{aligned}
\langle\mu_s,\mu_s\rangle_{\mathcal{H}} &=
\frac{1}{n_s^2}\sum_{i=1}^{n_s}\sum_{j=1}^{n_s}
\langle\phi(x^s_i),\phi(x^s_j)\rangle_{\mathcal{H}}
\;=\;\frac{1}{n_s^2}\sum_{i,j}k(x^s_i,x^s_j),\\[2mm]
\langle\mu_t,\mu_t\rangle_{\mathcal{H}} &=\frac{1}{n_t^2}\sum_{i,j}k(x^t_i,x^t_j),\\[2mm]
\langle\mu_s,\mu_t\rangle_{\mathcal{H}} &=\frac{1}{n_s n_t}\sum_{i,j}k(x^s_i,x^t_j).
\end{aligned}
$$

Hence

$$
\operatorname{MMD}^2(X_s,X_t)
=\frac{1}{n_s^2}\sum_{i,j}k(x^s_i,x^s_j) \;+\;
\frac{1}{n_t^2}\sum_{i,j}k(x^t_i,x^t_j) \;-\;
\frac{2}{n_s n_t}\sum_{i,j}k(x^s_i,x^t_j).
\tag{★}
$$

Step 2: Write (★) in matrix form

Stack the $n_s+n_t$ samples in one order $(x_1,\dots,x_{n_s+n_t})$.
Define $K$ as in the statement; define $L$ so that its $(p,q)$-entry is

$$L_{pq}=
\begin{cases}
\frac{1}{n_s^2} & p,q\leq n_s\\
\frac{1}{n_t^2} & p,q>n_s\\
-\frac{1}{n_s n_t} & \text{otherwise}.
\end{cases}$$

Then each double sum in (★) is a block‑wise Frobenius inner product $K_{u,v}\!:\!L_{u,v}$.
Collecting the three terms yields

$$\operatorname{MMD}^2(X_s,X_t)=\sum_{p,q}K_{pq}L_{pq}
\;=\;\langle K,L\rangle_F
\;=\;\operatorname{tr}(K^\top L)
\;=\;\operatorname{tr}(KL),$$

because $K$ is symmetric.

---

2  Population MMD as a dual norm

(a) Dual‑norm identity in a Hilbert space

For any $f,g\in\mathcal{H}$ the Cauchy–Schwarz inequality gives
$|\langle g,f\rangle|\le\|g\|\,\|f\|$.
Taking the supremum over all $g$ with $\|g\|\le1$:

$\sup_{\|g\|\le1}\langle g,f\rangle
\le\sup_{\|g\|\le1}\|g\|\,\|f\|=\|f\|$.

Equality is attained by choosing $g=f/\|f\|$ (or any unit‑vector multiple of $f$).
Therefore

$$
\boxed{\;\|f\|_{\mathcal{H}}=
\sup_{\|g\|_{\mathcal{H}}\le1}\langle g,f\rangle_{\mathcal{H}}\;}
\tag{1}
$$

---

(b) Variational form of MMD

Let $\mu_P=\mathbb{E}_{x\sim P}[\phi(x)]$ and $\mu_Q=\mathbb{E}_{x\sim Q}[\phi(x)]$.
Then

$\operatorname{MMD}(P,Q)
=\|\mu_P-\mu_Q\|_{\mathcal{H}}$.

Apply (1) with $f=\mu_P-\mu_Q$:

$\operatorname{MMD}(P,Q)
=\sup_{\|g\|\le1}\langle g,\mu_P-\mu_Q\rangle
=\sup_{\|g\|\le1}
\Bigl( \mathbb{E}_{x\sim P}[g(x)]-\mathbb{E}_{x\sim Q}[g(x)] \Bigr)$.

Rename $g$ to $f$ and we have the desired result:

$$\boxed{\;
\operatorname{MMD}(P,Q)=
\sup_{\|f\|_{\mathcal{H}}\le1}\;
\mathbb{E}_{x\sim P}[f(x)]
-\mathbb{E}_{x\sim Q}[f(x)]
\;}$$


$$
\begin{aligned}
\min_{\pi}\mathbb E_{\pi}[d]
&=1-\sum_x\min\{p(x),q(x)\}\\
&=\tfrac12\sum_x\bigl(p(x)+q(x)-2\min\{p(x),q(x)\}\bigr)\\
&=\tfrac12\sum_x|p(x)-q(x)|\\
&=\tfrac12\,\operatorname{TV}(D,D{\prime}).\\[-3pt]
\end{aligned}
$$