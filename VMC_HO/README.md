# Simulating Harmonic Oscillator using Variational Monte Carlo Method

The Variational Monte Carlo (VMC) method is used to estimate the ground state energy of a quantum system using the **variational principle**, **Monte Carlo integration**, and the **Metropolis algorithm**.

---

## 🔑 Key Equations

The variational energy is given by:

$$
E_T = \int \rho(R)\, E_L(R)\, dR
$$

where:

- Probability density:

$$
\rho(R) = \frac{\Psi_T^2(R)}{\int \Psi_T^2(R)\, dR}
$$

- Local energy:

$$
E_L(R) = \frac{\hat{H} \Psi_T(R)}{\Psi_T(R)}
$$

Monte Carlo estimate:

$$
E_T \approx \frac{1}{M} \sum_{i=1}^{M} E_L(R_i)
$$

---

## ⚙️ Algorithm Steps

### 1. Choose Trial Wavefunction

Select a trial wavefunction with variational parameters:

$$
\Psi_T(R; \alpha)
$$

Example (1D harmonic oscillator):

$$
\Psi_T(x) = e^{-\beta x^2}
$$

---

### 2. Initialize

- Choose initial configuration: $R_0$
- Set step size: $h$
- Set number of Monte Carlo steps: $N_{\text{total}}$
- Set thermalization steps: $N_{\text{therm}}$

---

### 3. Metropolis Sampling

For each Monte Carlo step:

1. Propose a new configuration:
    $$
    R' = R + \delta, \quad \delta \in [-h/2, h/2]
    $$

2. Compute acceptance ratio:
    $$
    A = \frac{|\Psi_T(R')|^2}{|\Psi_T(R)|^2}
    $$

3. Accept or reject:
   - If $A \ge 1$ → accept  
   - Else accept with probability $A$

4. Update configuration:
   - If accepted: $R \leftarrow R'$
   - Else: keep current $R$

**Note**: _In a MC step when we reject a movement/step, we set the current position to the position that was at the starting of the step. After rejection, we do need to calculate the energy of the current position. that should not be skipped over._

---

### 4. Thermalization

- Discard the first $N_{\text{therm}}$ samples  
- Ensures independence from initial conditions  

---

### 5. Compute Local Energy

For each sampled configuration $R_i$:

$$
E_L(R_i) = \frac{\hat{H} \Psi_T(R_i)}{\Psi_T(R_i)}
$$

Example (1D harmonic oscillator):

$$
E_L(x) = \beta + (0.5 - 2\beta^2)x^2
$$

---

### 6. Estimate Energy

$$
E_T = \frac{1}{M} \sum_{i=1}^{M} E_L(R_i)
$$

---

### 7. Compute Variance

$$
\sigma^2 = \frac{1}{M} \sum_{i=1}^{M} E_L^2(R_i) - E_T^2
$$

---

### 8. Optimize Variational Parameters

- Vary parameters (e.g. $\beta$)
- Minimize $E_T$
- Optimal parameters give best approximation to ground state

---

## 🧩 Compact Algorithm (Pseudocode)

```text
1. Initialize R
2. For i = 1 → N_total:
    2.1 Propose R'
    2.2 Compute A = |Ψ_T(R')|^2 / |Ψ_T(R)|^2
    2.3 Accept/reject using Metropolis rule
    2.4 If i > N_therm:
        2.4.1 Compute E_L(R)
        2.4.2 Accumulate energy
3. Compute average energy E_T
4. Compute variance σ²
5. Optimize parameters α```
