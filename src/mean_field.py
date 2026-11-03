"""Mean-field and Bogoliubov foundation; see derivations/mean_field_bogoliubov.md."""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import root


@dataclass
class Chain:
    N: int = 40
    J: float = 1.0
    U: float = -2e-4
    V: float = 0.0
    kappa: float = 1.0
    phi: float = np.pi / 3
    delta: float = 0.5
    drive: float = 38.0
    boundary: str = "obc"
    profile: str = "homogeneous"
    width: float = 5.0

    def __post_init__(self):
        if self.N < 3 or self.boundary not in ("obc", "pbc"):
            raise ValueError("Require N >= 3 and boundary obc or pbc")
        if self.kappa < 0 or self.profile not in ("homogeneous", "damped"):
            raise ValueError("Invalid loss rate or profile")
        if self.profile == "damped" and (self.boundary != "obc" or self.width <= 0):
            raise ValueError("Damped profiles require OBC and width > 0")
        self.C = np.zeros((self.N, self.N))
        self.T = np.zeros((self.N, self.N), dtype=complex)
        bonds = [(j, j + 1) for j in range(self.N - 1)]
        if self.boundary == "pbc":
            bonds.append((self.N - 1, 0))
        for j, k in bonds:
            self.C[j, k] = self.C[k, j] = 1
            self.T[k, j] = self.J * np.exp(1j * self.phi)
            self.T[j, k] = self.T[k, j].conjugate()
        s = np.ones(self.N)
        if self.profile == "damped":
            d = np.minimum(np.arange(1, self.N + 1), np.arange(self.N, 0, -1))
            s = np.tanh(d / self.width)
        self.detuning = self.delta * s**2
        self.pump = self.drive * s

    def force(self, alpha):
        """f in i d(alpha)/dt = f."""
        n = np.abs(alpha)**2
        return self.pump + self.T @ alpha + (
            self.detuning + self.U - 0.5j * self.kappa
            + 2 * self.U * n + self.V * (self.C @ n)
        ) * alpha

    def blocks(self, alpha):
        n = np.abs(alpha)**2
        A = self.T + np.diag(self.detuning + self.U - 0.5j * self.kappa
                            + 4 * self.U * n + self.V * (self.C @ n))
        A = A + self.V * self.C * np.outer(alpha, alpha.conjugate())
        B = np.diag(2 * self.U * alpha**2) + self.V * self.C * np.outer(alpha, alpha)
        return A, B

    def bogoliubov(self, alpha):
        A, B = self.blocks(alpha)
        return np.block([[A, B], [-B.conjugate(), -A.conjugate()]])

    def stationary(self, initial, residual_tol=1e-9):
        """Local root solve. Failure does not establish absence of a steady state."""
        def unpack(x):
            return x[:self.N] + 1j * x[self.N:]

        def fun(x):
            f = self.force(unpack(x))
            return np.r_[f.real, f.imag]

        def jac(x):
            A, B = self.blocks(unpack(x))
            return np.block([[(A + B).real, -(A - B).imag],
                             [(A + B).imag, (A - B).real]])

        initial = np.asarray(initial, dtype=complex)
        if initial.shape != (self.N,):
            raise ValueError("initial must have shape (N,)")
        result = root(fun, np.r_[initial.real, initial.imag], jac=jac, tol=1e-10)
        alpha = unpack(result.x)
        residual = float(np.max(np.abs(self.force(alpha))))
        return {"alpha": alpha, "residual": residual,
                "converged": bool(np.isfinite(residual) and residual < residual_tol),
                "solver_success": bool(result.success), "message": result.message}
