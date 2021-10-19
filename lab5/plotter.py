import matplotlib.pyplot as plt
import numpy as np
import subprocess


def analytic_solution(x, t):
    return np.sin(t) * np.cos(x)


def analytic_solution_fixed_time(t, h, count=None):
    res = np.array([analytic_solution(x, t) for x in np.arange(0, np.pi / 2, h)])
    if count and count != len(res):
        return res[-1:]
    return res


def analytic_solution_fixed_x(x, tau, count=None):
    res = np.array([analytic_solution(x, t) for t in np.arange(0, 10, tau)])
    if count and count != len(res):
        return res[-1:]
    return res


def error(analytic, numeric):
    return np.sum((analytic - numeric) ** 2) ** 0.5


def get_fixed_x_solution(solution, idx):
    res = []
    for k, v in solution.items():
        res.append(v[idx])
    return res


def run_prog(h: float, tau: float):
    return subprocess.run(['cmake-build-debug/lab5', str(h), str(tau)]).stdout


def read_file(file: str, tau: float) -> dict:
    res = dict()
    count = 0
    with open(file, 'r+') as f:
        line = f.readline()
        while line:
            res[tau * count] = np.array(list(map(float, line.split())))
            line = f.readline()
            count += 1
    return res


run_prog(0.2, 0.02)
impl_solution = read_file('impl.txt', 0.02)
expl_solution = read_file('expl.txt', 0.02)
exim_solution = read_file('exim.txt', 0.02)

x = np.arange(0, np.pi / 2, 0.2)

fig = plt.figure()
fig.suptitle('Solution at t = 0.5')
plt.plot(x, impl_solution[0.5], label='Implicit method')
plt.plot(x, expl_solution[0.5], label='Explicit method')
plt.plot(x, exim_solution[0.5], label='Crank-Nicolson method')
plt.plot(x, analytic_solution_fixed_time(0.5, 0.2), label='Analytic solution')
plt.legend()

fig = plt.figure()
fig.suptitle('Solution at x = 0.6')

t = np.arange(0, 10, 0.02)
plt.plot(t, get_fixed_x_solution(impl_solution, -5), label='Implicit method')
plt.plot(t, get_fixed_x_solution(expl_solution, -5), label='Explicit method')
plt.plot(t, get_fixed_x_solution(exim_solution, -5), label='Crank-Nicolson method')
plt.plot(t, analytic_solution_fixed_x(0.6, 0.02), label='Analytic solution')
plt.legend()

fig = plt.figure()

fig.suptitle('Error - time dependency, x = 0.6')

plt.plot(t, [error(analytic_solution_fixed_time(t, len(impl_solution[t])), impl_solution[t]) for t in
             np.arange(0, 10, 0.02)],
         label='Implicit method error')
plt.plot(t, [error(analytic_solution_fixed_time(t, len(impl_solution[t])), expl_solution[t]) for t in
             np.arange(0, 10, 0.02)],
         label='Explicit method error')
plt.plot(t, [error(analytic_solution_fixed_time(t, len(impl_solution[t])), exim_solution[t]) for t in
             np.arange(0, 10, 0.02)],
         label='Crank-Nicolson method error')

plt.legend()

fig = plt.figure()

fig.suptitle('Error - tau dependency, h = 0.64, t ≈ 5')

tau = np.arange(0.01, 0.21, 0.01)

tau_errors_expl = []
tau_errors_impl = []
tau_errors_exim = []

for i in tau:
    run_prog(0.64, i)
    impl_solution = read_file('impl.txt', i)
    expl_solution = read_file('expl.txt', i)
    exim_solution = read_file('exim.txt', i)

    t = 5
    current_t = 0
    current_diff = np.inf
    for k in impl_solution.keys():
        if abs(k - t) < current_diff:
            current_diff = abs(k - t)
            current_t = k

    tau_errors_expl.append(
        error(analytic_solution_fixed_time(current_t, 0.64, len(expl_solution[current_t])), expl_solution[current_t]))
    tau_errors_impl.append(
        error(analytic_solution_fixed_x(current_t, 0.64, len(expl_solution[current_t])), impl_solution[current_t]))
    tau_errors_exim.append(
        error(analytic_solution_fixed_x(current_t, 0.64, len(expl_solution[current_t])), exim_solution[current_t]))

plt.plot(tau, tau_errors_impl, label='Implicit method')
plt.plot(tau, tau_errors_expl, label='Explicit method')
plt.plot(tau, tau_errors_exim, label='Crank-Nicolson method')
plt.legend()

fig = plt.figure()

fig.suptitle('Error - h dependency, tau = 0.02, t = 5')

h_values = np.arange(0.2, 0.64, 0.02)

h_errors_expl = []
h_errors_impl = []
h_errors_exim = []

for i in h_values:
    run_prog(i, 1)
    impl_solution = read_file('impl.txt', 1)
    expl_solution = read_file('expl.txt', 1)
    exim_solution = read_file('exim.txt', 1)

    h_errors_expl.append(error(analytic_solution_fixed_time(5, len(expl_solution[5])), expl_solution[5]))
    h_errors_impl.append(error(analytic_solution_fixed_time(5, len(expl_solution[5])), impl_solution[5]))
    h_errors_exim.append(error(analytic_solution_fixed_time(5, len(expl_solution[5])), exim_solution[5]))

plt.plot(h_values, h_errors_impl, label='Implicit method')
plt.plot(h_values, h_errors_expl, label='Explicit method')
plt.plot(h_values, h_errors_exim, label='Crank-Nicolson method')
plt.legend()

plt.show()
