import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ─── Physical Constants (SI) ────────────────────────────────────────────────
hbar = 1.054571817e-34       # J·s
m_e  = 9.10938356e-31        # kg (electron mass)
# ─── Simulation Parameters ──────────────────────────────────────────────────
nx = 200                     # number of spatial points
L  = 1.0e-6                  # box length in meters (1 μm)
dx = L / nx                  # spatial step (m)
dt = 1.0e-18                 # time step (s)
tl = 2000                    # total number of time steps

# Derived grid
x = np.linspace(-L/2, L/2 - dx, nx)

# ─── Initial Gaussian Wavepacket ───────────────────────────────────────────
x0    = -L/4                 # initial center (m)
sigma = 1.0e-7               # width (m)
p0    = 1.0e-14              # initial momentum (kg·m/s)
k0    = p0 / hbar            # wave number (1/m)

# Construct initial psi
psi = np.exp(-0.5*((x - x0)/sigma)**2) * np.exp(1j*k0*x)
# Normalize
norm = np.sqrt(np.sum(np.abs(psi)**2) * dx)
psi /= norm

# ─── Crank–Nicolson Coefficients ───────────────────────────────────────────
alpha = 1j * hbar * dt / (2.0 * m_e * dx**2)

# A tri‐diagonal (1 + 2α on diag, −α on off‐diag)
a_main = np.full(nx, 1.0 + 2.0*alpha, dtype=complex)
a_off  = np.full(nx-1, -alpha,      dtype=complex)

# B tri‐diagonal (1 − 2α on diag, +α on off‐diag)
b_main = np.full(nx, 1.0 - 2.0*alpha, dtype=complex)
b_off  = np.full(nx-1,  alpha,       dtype=complex)

# Hard‐wall boundaries: only cut the very first/last linkimport numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ─── Physical Constants (SI) ────────────────────────────────────────────────
hbar = 1.054571817e-34      # J·s
m_e  = 9.10938356e-31       # kg

# ─── Simulation Parameters ──────────────────────────────────────────────────
nx = 500                    # spatial grid points
L  = 1e-6                   # box length (m)
dx = L/nx                   # spatial step (m)

dt = 1e-15                  # time step (s) – big enough to see motion
tl = 10000                   # number of steps

# Derived grid
x = np.linspace(-L/2, L/2-dx, nx)

# ─── Initial Gaussian Wavepacket ───────────────────────────────────────────
x0    = -L/4                # center (m)
sigma = 1e-7                # width (m)
p0    = 5e-24               # momentum (kg·m/s) – choose so v = p/m ~5e6 m/s
k0    = p0 / hbar           # wave number (1/m)

# initial psi
psi = np.exp(-0.5*((x-x0)/sigma)**2) * np.exp(1j*k0*x)
psi /= np.sqrt(np.sum(np.abs(psi)**2)*dx)

# ─── Crank–Nicolson Coefficients ───────────────────────────────────────────
alpha = 1j*hbar*dt/(4*m_e*dx**2)

# Build A (I + i*H*dt/2) & B (I - i*H*dt/2) as tri−diagonals
a_main = np.full(nx, 1+2*alpha, dtype=complex)
a_off  = np.full(nx-1, -alpha,    dtype=complex)
b_main = np.full(nx, 1-2*alpha, dtype=complex)
b_off  = np.full(nx-1,  alpha,  dtype=complex)

# Hard‐wall BC: only cut the very first & last off‐diagonal link
a_main[0] = a_main[-1] = 1.0
b_main[0] = b_main[-1] = 1.0
a_off[0]    = 0.0
a_off[-1]   = 0.0
b_off[0]    = 0.0
b_off[-1]   = 0.0

def thomas_solve(a_diag, a_off, rhs):
    """Thomas algorithm to solve tridiagonal A·x = rhs."""
    n = len(a_diag)
    diag = a_diag.copy()
    off  = a_off.copy()
    b    = rhs.copy()
    # forward
    for i in range(1, n):
        w = off[i-1]/diag[i-1]
        diag[i] -= w*off[i-1]
        b[i]    -= w*b[i-1]
    # back
    x_new = np.zeros(n, dtype=complex)
    x_new[-1] = b[-1]/diag[-1]
    for i in range(n-2, -1, -1):
        x_new[i] = (b[i] - off[i]*x_new[i+1]) / diag[i]
    return x_new

# ─── Time Evolution ────────────────────────────────────────────────────────
psi_history = np.zeros((tl, nx), dtype=complex)
for t in range(tl):
    # Build RHS = B·psi
    rhs = np.zeros(nx, dtype=complex)
    rhs[1:-1] = (
        b_off[:-1] * psi[:-2] +
        b_main[1:-1] * psi[1:-1] +
        b_off[1:]   * psi[2:]
    )
    # enforce boundary values
    rhs[0], rhs[-1] = psi[0], psi[-1]
    # solve for next psi
    psi = thomas_solve(a_main, a_off, rhs)
    psi_history[t] = psi

# ─── Animate density + Re & Im on twin axes ────────────────────────────────
# ─── Animate density + Re & Im on twin axes ────────────────────────────────
fig, ax1 = plt.subplots(figsize=(8,4))
ax2 = ax1.twinx()

line_d, = ax1.plot(x*1e6, np.abs(psi_history[0])**2, 'b-', lw=2, label=r'$|\psi|^2$')
line_r, = ax2.plot(x*1e6, psi_history[0].real,  'r--', label='Re[ψ]')
line_i, = ax2.plot(x*1e6, psi_history[0].imag,  'g:',  label='Im[ψ]')

ax1.set_xlabel('x (μm)')
ax1.set_ylabel('Density')
ax2.set_ylabel('Wavefunction')

ax1.set_xlim(x.min()*1e6, x.max()*1e6)
ax1.set_ylim(0, (np.abs(psi_history)**2).max()*1.1)
ax2.set_ylim(psi_history.real.min()*1.1, psi_history.real.max()*1.1)

ax1.set_title('Moving Electron Wavepacket in 1 µm Box')
ax1.legend(loc='upper left')
ax2.legend(loc='upper right')

def update(frame):
    ψ = psi_history[frame]
    line_d.set_ydata(np.abs(ψ)**2)
    line_r.set_ydata(ψ.real)
    line_i.set_ydata(ψ.imag)
    return line_d, line_r, line_i

# Only every 5th time‐step:
frame_steps = range(0, tl, 10)

ani = FuncAnimation(
    fig, update,
    frames=frame_steps,     # instead of range(tl)
    interval=20, blit=True
)
#ani.save('wavepacket_SI_fixed.gif', writer='pillow', fps=30)
plt.show()

print('Saved: wavepacket_SI_fixed.gif')

for arr in (a_main, b_main):
    arr[0] = arr[-1] = 1.0
for arr in (a_off, b_off):
    arr[0] = arr[-1] = 0.0

def thomas_solve(a_diag, a_off, rhs):
    """Thomas algorithm for A x = rhs."""
    n = len(a_diag)
    diag = a_diag.copy()
    off  = a_off.copy()
    b    = rhs.copy()
    # forward elimination
    for i in range(1, n):
        w = off[i-1] / diag[i-1]
        diag[i] -= w * off[i-1]
        b[i]    -= w * b[i-1]
    # back substitution
    x_new = np.zeros(n, dtype=complex)
    x_new[-1] = b[-1] / diag[-1]
    for i in range(n-2, -1, -1):
        x_new[i] = (b[i] - off[i] * x_new[i+1]) / diag[i]
    return x_new

# ─── Time–Evolution Storage ────────────────────────────────────────────────
psi_history = np.zeros((tl, nx), dtype=complex)

for t in range(tl):
    # Build RHS = B * psi
    rhs = np.zeros(nx, dtype=complex)
    rhs[1:-1] = (
        b_off[:-1] * psi[:-2] +
        b_main[1:-1] * psi[1:-1] +
        b_off[1:]   * psi[2:]
    )
    rhs[0], rhs[-1] = psi[0], psi[-1]
    # Solve A ψ_new = RHS
    psi = thomas_solve(a_main, a_off, rhs)
    psi_history[t] = psi

# ─── Animation (Density + Re & Im Parts) ──────────────────────────────────
fig, ax = plt.subplots()
line_d, = ax.plot(x*1e6, np.abs(psi_history[0])**2, lw=2, label=r'$|\psi|^2$')
line_r, = ax.plot(x*1e6, psi_history[0].real,   ls='--', label='Re[$\\psi$]')
line_i, = ax.plot(x*1e6, psi_history[0].imag,   ls=':',  label='Im[$\\psi$]')
ax.set_xlim((x.min()*1e6, x.max()*1e6))
ax.set_xlabel('x (μm)')
ax.legend()
ax.set_title('Electron Wavepacket in 1 μm Box')

def update(frame):
    ψ = psi_history[frame]
    line_d.set_ydata(np.abs(ψ)**2)
    line_r.set_ydata(ψ.real)
    line_i.set_ydata(ψ.imag)
    return line_d, line_r, line_i

ani = FuncAnimation(fig, update, frames=tl, interval=20, blit=True)
#ani.save('wavepacket_SI.gif', writer='pillow', fps=30)
plt.show()

print('Saved SI‐unit animation as wavepacket_SI.gif')
