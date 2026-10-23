"""Least-squares solvers for the fits. A problem supplies evaluate(x) -> (residual, Jacobian, ...) and fun(x)."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np


def damped_gauss_newton(prob, x, lower, upper, max_nfev=1000, ftol=1e-10, gtol=1e-6, lazy_jacobian=False):
    """Gauss-Newton with SVD damping of each trial step, in Jacobian-column-scaled coordinates, within bounds
    (the core fit). With `lazy_jacobian` a trial costs only its residuals; the Jacobian is formed once accepted.

    Stops when the scaled projected gradient falls below gtol (status 1), after eight accepted steps of relative
    gain below ftol (2), or when no damped step improves (3)."""
    x = x.copy()
    r, J, _ = prob.evaluate(x)
    lam, nfev, njev, accepted, stagnant, rejected = 1e-5, 1, 1, 0, 0, 0
    c = float(r@r)
    status, message = 0, "Maximum evaluations reached"
    while nfev < max_nfev:
        gradient = J.T@r
        active = ((x-lower) < 1e-9) & (gradient > 0) | ((upper-x) < 1e-9) & (gradient < 0)
        norms = np.maximum(np.sqrt(np.sum(J*J, axis=0)), 1e-15)
        free = ~active & (norms > 1e-14)
        optimality = float(np.max(abs(gradient[free]/norms[free]))) if np.any(free) else 0.
        if optimality < gtol:
            status, message = 1, "Scaled projected gradient converged"
            break
        U, s, Vt = np.linalg.svd(J[:, free]/norms[free], full_matrices=False)
        ub = U.T@r
        old, improved, backtracking = c, False, False
        while nfev < max_nfev and lam < 1e12:
            step = np.zeros_like(x)
            step[free] = -(Vt.T@(s/(s*s+lam)*ub))/norms[free]
            trial = np.clip(x+step, lower, upper)
            if lazy_jacobian:
                rn, Jn = prob.fun(trial), None
            else:
                rn, Jn, _ = prob.evaluate(trial)
                njev += 1
            cn = float(rn@rn)
            nfev += 1
            if cn < c:
                if lazy_jacobian:
                    rn, Jn, _ = prob.evaluate(trial)
                    njev += 1
                    cn = float(rn@rn)        # the full evaluation decides acceptance
                    if cn >= c:
                        rejected += 1
                        lam *= 10
                        backtracking = True
                        continue
                predicted = float(-(gradient@(trial-x))*2-np.linalg.norm(J@(trial-x))**2)
                rho = (c-cn)/predicted if predicted > 0 else 0.
                x, r, J, c = trial, rn, Jn, cn
                accepted += 1
                improved = True
                lam = max(lam*(.2 if rho > .7 else .7 if rho > .1 else 2.), 1e-16)
                break
            rejected += 1
            lam *= 10
            backtracking = True
        if not improved:
            status, message = 3, "No improving damped step found"
            break
        stagnant = stagnant+1 if old-c < ftol*max(old, 1.) else 0
        if stagnant >= 8:
            status, message = 2, "Eight accepted steps had negligible relative improvement"
            break
    gradient = J.T@r
    active = ((x-lower) < 1e-9) & (gradient > 0) | ((upper-x) < 1e-9) & (gradient < 0)
    norms = np.maximum(np.linalg.norm(J, axis=0), 1e-15)
    free = ~active & (norms > 1e-14)
    optimality = float(np.max(abs(gradient[free]/norms[free]))) if np.any(free) else 0.
    return SimpleNamespace(x=x, success=status in (1, 2), status=status, message=message, nfev=nfev, njev=njev,
                           cost=c, optimality=optimality, accepted_steps=accepted, rejected_steps=rejected,
                           active_mask=np.where(x-lower < 1e-8, -1, np.where(upper-x < 1e-8, 1, 0)))


def marquardt(prob, x, max_nfev, scale, tol=1e-9):
    """Levenberg-Marquardt, diagonally scaled, started heavily damped (the continuation).

    The damping floor |g_i|/scale_i keeps a parameter the data barely see from moving more than ~scale_i/lambda per
    step. Every evaluation returns its Jacobian, so each trial is a full iterate; rejected trials raise the damping.
    """
    r, J = prob.evaluate(x)
    cost, lam, nfev, small = float(r@r), 1., 1, 0
    while nfev < max_nfev:
        A, g = J.T@J, J.T@r
        D = np.maximum(np.maximum(np.diag(A), abs(g)/scale), 1e-12)
        step = np.linalg.solve(A+lam*np.diag(D), -g)
        r2, J2 = prob.evaluate(x+step)
        nfev += 1
        cost2 = float(r2@r2)
        if np.isfinite(cost2) and cost2 < cost:
            small = small+1 if (cost-cost2)/cost < tol else 0
            x, r, J, cost, lam = x+step, r2, J2, cost2, max(lam/3, 1e-9)
            if small >= 2:
                return dict(x=x, cost=cost, nfev=nfev, converged=True)
        else:
            lam *= 4
            if lam > 1e12:
                return dict(x=x, cost=cost, nfev=nfev, converged=True)
    return dict(x=x, cost=cost, nfev=nfev, converged=False)


def gauss_newton(prob, x, max_nfev, tol=1e-9):
    """Undamped Gauss-Newton steps on the column-scaled Jacobian, halved until the cost falls (the continuation's
    finish). Near the minimum Marquardt's damping floor lets each step take only a few per cent of the way along
    weakly determined directions; the linear model is accurate there. Converged when a full step would gain less
    than tol of the cost."""
    r, J = prob.evaluate(x)
    cost, nfev = float(r@r), 1
    while nfev < max_nfev:
        c = np.linalg.norm(J, axis=0)
        c[c == 0] = 1.
        step = -np.linalg.lstsq(J/c, r, rcond=None)[0]/c
        if cost-float(np.sum((r+J@step)**2)) < tol*cost:
            return dict(x=x, cost=cost, nfev=nfev, converged=True)
        for _ in range(8):
            r2, J2 = prob.evaluate(x+step)
            nfev += 1
            if np.isfinite(r2@r2) and float(r2@r2) < cost:
                break
            step /= 2
        else:
            return dict(x=x, cost=cost, nfev=nfev, converged=False)
        x, r, J, cost = x+step, r2, J2, float(r2@r2)
    return dict(x=x, cost=cost, nfev=nfev, converged=False)
