"""Run with: python src/validate_mean_field.py."""
import numpy as np
from scipy.optimize import linear_sum_assignment
from mean_field import Chain


def main():
    rng = np.random.default_rng(20261103)
    for boundary, profile in [("obc", "homogeneous"), ("obc", "damped"),
                              ("pbc", "homogeneous")]:
        model = Chain(N=8, U=-0.03, V=-0.02, boundary=boundary, profile=profile)
        a = rng.normal(size=8) + 1j * rng.normal(size=8)
        d = rng.normal(size=8) + 1j * rng.normal(size=8)
        A, B = model.blocks(a)
        h = 1e-5
        finite = (model.force(a + h*d) - model.force(a - h*d)) / (2*h)
        error = np.max(np.abs(finite - A@d - B@d.conjugate()))
        assert error < 1e-8, error
        np.testing.assert_allclose(B, B.T, atol=1e-12)
        np.testing.assert_allclose(A + 0.5j*model.kappa*np.eye(8),
                                   (A + 0.5j*model.kappa*np.eye(8)).conj().T)
        print(boundary, profile, "linearization error", error)

    for V in (0.0, -2e-4):
        m = Chain(N=12, V=V, boundary="pbc")
        c = m.delta + m.U + 2*m.J*np.cos(m.phi)
        g = 2*(m.U + m.V)
        roots = np.roots([g*g, 2*c*g, c*c + (m.kappa/2)**2, -m.drive**2])
        for n in roots:
            if abs(n.imag) > 1e-7 or n.real < 0:
                continue
            n = n.real
            a = -m.drive/(c + g*n - 0.5j*m.kappa)
            state = np.full(m.N, a)
            assert np.max(np.abs(m.force(state))) < 1e-8
            result = m.stationary(state * 1.001)
            assert result["converged"], result
            real_space = np.linalg.eigvals(m.bogoliubov(state))
            momentum = []
            for q in 2*np.pi*np.arange(m.N)/m.N:
                def aq(k):
                    return (m.delta + m.U + 4*m.U*n + 2*m.V*n
                            + 2*m.J*np.cos(k-m.phi) + 2*m.V*n*np.cos(k)
                            - 0.5j*m.kappa)
                bq = 2*a*a*(m.U + m.V*np.cos(q))
                momentum.extend(np.linalg.eigvals([[aq(q), bq],
                                                    [-bq.conjugate(), -aq(-q).conjugate()]]))
            distances = np.abs(real_space[:, None] - np.array(momentum)[None, :])
            rows, cols = linear_sum_assignment(distances)
            error = np.max(distances[rows, cols])
            assert error < 1e-8, error
            print("PBC V=", V, "density=", n, "spectrum error=", error,
                  "growth rate=", np.max(real_space.imag))

    m = Chain(N=8, U=0, V=0, drive=1)
    exact = np.linalg.solve(m.T + np.diag(m.detuning-0.5j*m.kappa), -m.pump)
    result = m.stationary(np.zeros(m.N, dtype=complex))
    assert result["converged"]
    np.testing.assert_allclose(result["alpha"], exact, atol=1e-10)
    print("Linear driven-chain benchmark passed; all checks passed.")


if __name__ == "__main__":
    main()
