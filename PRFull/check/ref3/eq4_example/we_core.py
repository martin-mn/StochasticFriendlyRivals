"""Core for the worked-example checks (independent code).
State s = 4a+b, a = last outcome, b = outcome before, focal view; outcome o = 2*fD + qD
(CC=0, CD=1, DC=2, DD=3). Co-player sees CD<->DC exchanged.
"""
import numpy as np

SW = [0, 2, 1, 3]
P_EX = [1, 0.5, 1, 0.5, 0, 0, 1, 0.5, 1, 0.25, 1, 0.25, 0.25, 0, 0.25, 0]


def coplayer_state(s):
    a, b = divmod(s, 4)
    return 4 * SW[a] + SW[b]


def all_binary():
    codes = np.arange(65536)
    return ((codes[:, None] >> np.arange(16)[None, :]) & 1).astype(float)


def transition_batch(p, Q, eps):
    """p: (16,) focal; Q: (B,16) co-players. Returns (B,16,16)."""
    B = Q.shape[0]
    p = np.asarray(p, float)
    pc = eps + (1 - 2 * eps) * p                      # focal coop prob per state
    qidx = [coplayer_state(s) for s in range(16)]
    qc = eps + (1 - 2 * eps) * Q[:, qidx]             # co-player coop prob per (focal) state
    T = np.zeros((B, 16, 16))
    for s in range(16):
        a = s // 4
        for fD in (0, 1):
            pf = pc[s] if fD == 0 else 1 - pc[s]
            for qD in (0, 1):
                pq = qc[:, s] if qD == 0 else 1 - qc[:, s]
                o = 2 * fD + qD
                T[:, s, 4 * o + a] += pf * pq
    return T


def gth_batch(T):
    """Stationary distributions by GTH (no subtractions). T (B,n,n)."""
    T = T.copy()
    B, n, _ = T.shape
    for k in range(n - 1, 0, -1):
        s = T[:, k, :k].sum(axis=1)                   # (B,)
        T[:, :k, :k] += T[:, :k, k][:, :, None] * T[:, k, :k][:, None, :] / s[:, None, None]
        T[:, k, :k] /= s[:, None]   # store normalised row for back substitution
    pi = np.zeros((B, n))
    pi[:, 0] = 1.0
    for k in range(1, n):
        pi[:, k] = (pi[:, :k] * T[:, :k, k]).sum(axis=1)
    # careful: back substitution uses column k of the reduced matrices; recompute properly below
    return pi


def gth(T):
    """Single-chain GTH, returns stationary distribution. Standard formulation."""
    A = np.array(T, float).copy()
    n = A.shape[0]
    for k in range(n - 1, 0, -1):
        s = A[k, :k].sum()
        A[:k, k] /= s
        A[:k, :k] += np.outer(A[:k, k], A[k, :k])
    pi = np.zeros(n)
    pi[0] = 1.0
    for k in range(1, n):
        pi[k] = pi[:k] @ A[:k, k]
    return pi / pi.sum()


def gth_b(T):
    """Batched standard GTH."""
    A = T.copy()
    B, n, _ = A.shape
    for k in range(n - 1, 0, -1):
        s = A[:, k, :k].sum(axis=1)
        A[:, :k, k] /= s[:, None]
        A[:, :k, :k] += A[:, :k, k][:, :, None] * A[:, k, :k][:, None, :]
    pi = np.zeros((B, n))
    pi[:, 0] = 1.0
    for k in range(1, n):
        pi[:, k] = np.einsum('bi,bi->b', pi[:, :k], A[:, :k, k])
    return pi / pi.sum(axis=1, keepdims=True)


def margins(p, Q, eps=1e-12):
    T = transition_batch(p, Q, eps)
    pi = gth_b(T)
    w = pi.reshape(-1, 4, 4).sum(axis=2)   # outcome of last round
    return w[:, 2] - w[:, 1], w
