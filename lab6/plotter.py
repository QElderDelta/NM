import matplotlib.pyplot as plt
import numpy as np
import subprocess


def analytic_solution(x, t):
    return np.exp(-x - t) * np.cos(x) * np.cos(2 * t)


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
    return np.sum((analytic - numeric) ** 2) ** 0.5 / len(analytic)

def error_x(analytic, numeric):
    return (analytic - numeric) ** 2 / len(numeric)


def get_fixed_x_solution(solution, idx):
    res = []
    for k, v in solution.items():
        res.append(v[idx])
    return res


def run_prog(h: float, tau: float, t_max: float):
    input = str(h) + ' ' + str(tau) + ' ' + str(t_max)
    return subprocess.run(['cmake-build-debug/lab6'], shell=True, input=input.encode('utf-8')).stdout


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


h = 0.3
tau = 0.2

run_prog(h, tau, 10)
impl_solution = read_file('impl.txt', tau)
expl_solution = read_file('expl.txt', tau)

x = np.arange(0, np.pi / 2, h)

t = 1

fig = plt.figure()
fig.suptitle('Solution at t = {}'.format(t))
plt.plot(x, impl_solution[t], label='Implicit method')
plt.plot(x, expl_solution[t], label='Explicit method')
plt.plot(x, analytic_solution_fixed_time(t, h), label='Analytic solution')
plt.legend()

t = 5

fig = plt.figure()
fig.suptitle('Solution at t = {}'.format(t))
plt.plot(x, impl_solution[t], label='Implicit method')
plt.plot(x, expl_solution[t], label='Explicit method')
plt.plot(x, analytic_solution_fixed_time(t, h), label='Analytic solution')
plt.legend()

t = 7

fig = plt.figure()
fig.suptitle('Solution at t = {}'.format(t))
plt.plot(x, impl_solution[t], label='Implicit method')
plt.plot(x, expl_solution[t], label='Explicit method')
plt.plot(x, analytic_solution_fixed_time(t, h), label='Analytic solution')
plt.legend()

tau = 0.01

run_prog(h, tau, 10)
impl_solution = read_file('impl.txt', tau)
expl_solution = read_file('expl.txt', tau)

x = np.arange(0, np.pi / 2, h)

t = 1

fig = plt.figure()
fig.suptitle('Solution at t = {}'.format(t))
plt.plot(x, impl_solution[t], label='Implicit method')
plt.plot(x, expl_solution[t], label='Explicit method')
plt.plot(x, analytic_solution_fixed_time(t, h), label='Analytic solution')
plt.legend()

t = 5

fig = plt.figure()
fig.suptitle('Solution at t = {}'.format(t))
plt.plot(x, impl_solution[t], label='Implicit method')
plt.plot(x, expl_solution[t], label='Explicit method')
plt.plot(x, analytic_solution_fixed_time(t, h), label='Analytic solution')
plt.legend()

t = 7

fig = plt.figure()
fig.suptitle('Solution at t = {}'.format(t))
plt.plot(x, impl_solution[t], label='Implicit method')
plt.plot(x, expl_solution[t], label='Explicit method')
plt.plot(x, analytic_solution_fixed_time(t, h), label='Analytic solution')
plt.legend()

fig = plt.figure()
fig.suptitle('Solution at x = 0.6')

t = np.arange(0, 10, tau)
plt.plot(t, get_fixed_x_solution(impl_solution, 2), label='Implicit method')
plt.plot(t, get_fixed_x_solution(expl_solution, 2), label='Explicit method')
plt.plot(t, analytic_solution_fixed_x(0.6, tau), label='Analytic solution')
plt.legend()

fig = plt.figure()

fig.suptitle('Error - time dependency, x = 0.6')

plt.plot(t, [error(analytic_solution_fixed_time(t, len(impl_solution[t])), impl_solution[t]) for t in
             np.arange(0, 10, tau)],
         label='Implicit method error')
plt.plot(t, [error(analytic_solution_fixed_time(t, len(impl_solution[t])), expl_solution[t]) for t in
             np.arange(0, 10, tau)],
         label='Explicit method error')

plt.legend()

# fig = plt.figure()
#
# fig.suptitle('Error, t = 5')
#
# t = 5
#
# plt.plot(x, error_x(analytic_solution_fixed_time(t, h), impl_solution[t]), label='impl')


# fig = plt.figure()
#
# fig.suptitle('Error - tau dependency, h = 0.64, t ≈ 5')
#
# tau = np.arange(0.1, 2, 0.2)
#
# tau_errors_expl = []
# tau_errors_impl = []
#
# for i in tau:
#     run_prog(0.3, i, 10)
#     impl_solution = read_file('impl.txt', i)
#     expl_solution = read_file('expl.txt', i)
#
#     t = 6
#     current_t = 0
#     current_diff = np.inf
#     for k in impl_solution.keys():
#         if abs(k - t) < current_diff:
#             current_diff = abs(k - t)
#             current_t = k
#
#     tau_errors_expl.append(
#         error(analytic_solution_fixed_time(current_t, 0.3, len(expl_solution[current_t])), expl_solution[current_t]))
#     tau_errors_impl.append(
#         error(analytic_solution_fixed_x(current_t, 0.3, len(expl_solution[current_t])), impl_solution[current_t]))
#
# plt.plot(tau, tau_errors_impl, label='Implicit method')
# plt.plot(tau, tau_errors_expl, label='Explicit method')
# plt.legend()

fig = plt.figure()

fig.suptitle('Error - h dependency, tau = 0.5, t = 5')

h_values = np.arange(0.2, 0.64, 0.02)

h_errors_expl = []
h_errors_impl = []

for i in h_values:
    run_prog(i, 0.5, 6)
    impl_solution = read_file('impl.txt', 0.5)
    expl_solution = read_file('expl.txt', 0.5)

    h_errors_expl.append(error(analytic_solution_fixed_time(5, len(expl_solution[5])), expl_solution[5]))
    h_errors_impl.append(error(analytic_solution_fixed_time(5, len(expl_solution[5])), impl_solution[5]))

plt.plot(h_values, h_errors_impl, label='Implicit method')
plt.plot(h_values, h_errors_expl, label='Explicit method')
plt.legend()

plt.show()
