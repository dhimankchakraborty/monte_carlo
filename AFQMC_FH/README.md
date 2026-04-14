# Auxiliary-Field Quantum Monte Carlo (AFQMC)

## Slater Determinant
A **Slater determinant** is the simplest many-body wavefunction used to describe a system of **fermions** (e.g., electrons), ensuring the required **antisymmetry** under particle exchange. In AFQMC, slater determinants are used as a simplest representation of many body wave function.

---

### Definition
A Slater determinant is constructed from (M) single-particle orbitals:
$$|\phi\rangle = \hat{\phi}_1^\dagger \hat{\phi}_2^\dagger \cdots \hat{\phi}_M^\dagger |0\rangle  $$
Each orbital is expanded in a chosen single-particle basis:
$$ \hat{\phi}_m^\dagger = \sum_i \phi_{i,m} , c_i^\dagger  $$

---

### Matrix Representation
The Slater determinant is represented by a matrix:
$$
\Phi \in \mathbb{C}^{N \times M}  
$$

- $N$: number of basis states (e.g., lattice sites)
- $M$: number of fermions
- Each **column** = one single-particle orbital
$$
\Phi =  
\left(\begin{array}{cc}
\phi_{1,1} &\phi_{1,2} &\cdots \\  
\phi_{2,1} &\phi_{2,2} &\cdots \\  
\vdots & \vdots & \ddots  
\end{array}\right)
$$
---
### Key Properties
- Columns (orbitals) are typically **orthonormal**:  
$$ 
    \Phi^\dagger \Phi = I$$
- Fermions are **indistinguishable** → columns do NOT correspond to specific particles
- The wavefunction in real space is:  
$$
    \Psi(R) = \det[\phi_m(r_i)]  
$$
- Overlap between two Slater determinants:  
$$ \langle \phi | \phi' \rangle = \det(\Phi^\dagger \Phi')  $$
---
### Important Note
- If the system contains different species of spin, then we need different slater determinants for each spin species and the final many body state is tensor product of all those slater determinants.
- **Initialization**
	- There are many possible way to initialize a Slater Determinant wave function
	- One possible way is to use **energy Eigen vectors** of the **non-interacting** counterpart of the interacting system that we are trying to simulate
	- In this formalism, each column will be one Eigen vector, starting with the ground state at the leftmost column and so on.
---
### Summary

- Slater determinant = antisymmetric fermionic wavefunction
- Represented as an $N \times M$ orbital matrix
- Foundation of mean-field methods and AFQMC
- Enables efficient many-body simulations via linear algebra