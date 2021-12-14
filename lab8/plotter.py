import matplotlib.pyplot as plt
import numpy as np
import subprocess


def analytic_solution(x, y, t):
    return x * y * np.cos(t)


def analytic_solution_fixed_yt(hx, y, t, count=None):
    res = np.array([analytic_solution(x, y, t) for x in np.arange(0, 1 + hx, hx)])
    if count and count != len(res):
        return res[-1:]
    return res


def analytic_solution_fixed_t(hx, hy, t, count=None):
    return np.array([[x * y * np.cos(t) for x in np.arange(0, 1 + hx, hx)] for y in np.arange(0, 1 + hy, hy)])


def error(analytic, numeric):
    return np.sum((analytic - numeric) ** 2) ** 0.5 / analytic.size


def get_fixed_x_solution(solution, idx):
    res = []
    for k, v in solution.items():
        res.append(v[idx])
    return res


def run_prog(hx: float, hy: float, tau: float, t_end: float):
    input = str(hx) + ' ' + str(hy) + ' ' + str(tau) + ' ' + str(t_end)
    return subprocess.run(['cmake-build-debug/lab8'], shell=True, input=input.encode('utf-8')).stdout


def read_file(file: str, hy: float, tau: float) -> dict:
    res = dict()
    curr_t = 0
    count = 0
    curr_grid = []
    with open(file, 'r+') as f:
        line = f.readline()
        while line:
            if count == int((1 / hy) + 1):
                count = 0
                res[curr_t] = curr_grid
                curr_t += tau
                curr_grid = []
            curr_grid.append(list(map(float, line.split())))
            line = f.readline()
            count += 1
    return res


def get_closest_to(sol, y):
    diff = 10
    key = -1
    for k in sol.keys():
        if abs(y - k) < diff:
            key = k
            diff = abs(y - k)
    return key


def read_eps_file(file):
    res = []
    with open(file, 'r') as f:
        line = f.readline()
        while line:
            res += map(float, line.split())
            line = f.readline()
    return np.array(res)


hx = 0.1
hy = 0.1
tau = 0.1
t_max = 5

run_prog(hx, hy, tau, t_max)
var_solution = read_file('var.txt', hy, tau)
step_solution = read_file('frac.txt', hy, tau)

x = np.arange(0, 1 + hx, hx)

t = get_closest_to(var_solution, 1)
y = 9

fig = plt.figure()
fig.suptitle(f'Solution at t = {t:5.4f}, y = {y * hy:5.4f}')
plt.plot(x, var_solution[t][y], label='Varying directions method')
plt.plot(x, step_solution[t][y], label='Fractional steps method')
plt.plot(x, analytic_solution_fixed_yt(hx, y * hy, t), label='Analytic solution')
plt.legend()

t = get_closest_to(var_solution, 2)
y = 9

fig = plt.figure()
fig.suptitle(f'Solution at t = {t:5.4f}, y = {y * hy:5.4f}')
plt.plot(x, var_solution[t][y], label='Varying directions method')
plt.plot(x, step_solution[t][y], label='Fractional steps method')
plt.plot(x, analytic_solution_fixed_yt(hx, y * hy, t), label='Analytic solution')
plt.legend()

t = get_closest_to(var_solution, 3)
y = 9

fig = plt.figure()
fig.suptitle(f'Solution at t = {t:5.4f}, y = {y * hy:5.4f}')
plt.plot(x, var_solution[t][y], label='Varying directions method')
plt.plot(x, step_solution[t][y], label='Fractional steps method')
plt.plot(x, analytic_solution_fixed_yt(hx, y * hy, t), label='Analytic solution')
plt.legend()

t = get_closest_to(var_solution, 1)
y = 7

fig = plt.figure()
fig.suptitle(f'Solution at t = {t:5.4f}, y = {y * hy:5.4f}')
plt.plot(x, var_solution[t][y], label='Varying directions method')
plt.plot(x, step_solution[t][y], label='Fractional steps method')
plt.plot(x, analytic_solution_fixed_yt(hx, y * hy, t), label='Analytic solution')
plt.legend()

t = get_closest_to(var_solution, 1)
y = 8

fig = plt.figure()
fig.suptitle(f'Solution at t = {t:5.4f}, y = {y * hy:5.4f}')
plt.plot(x, var_solution[t][y], label='Varying directions method')
plt.plot(x, step_solution[t][y], label='Fractional steps method')
plt.plot(x, analytic_solution_fixed_yt(hx, y * hy, t), label='Analytic solution')
plt.legend()

t = get_closest_to(var_solution, 1)
y = 9

fig = plt.figure()
fig.suptitle(f'Solution at t = {t:5.4f}, y = {y * hy:5.4f}')
plt.plot(x, var_solution[t][y], label='Varying directions method')
plt.plot(x, step_solution[t][y], label='Fractional steps method')
plt.plot(x, analytic_solution_fixed_yt(hx, y * hy, t), label='Analytic solution')
plt.legend()

fig = plt.figure()
fig.suptitle(f'Error-time dependency, tau={tau:5.4f}')
plt.plot([k for k in var_solution.keys()],
         [error(analytic_solution_fixed_t(hx, hy, t), var_solution[t]) for t in var_solution.keys()],
         label='Varying directions method error')
plt.plot([k for k in step_solution.keys()],
         [error(analytic_solution_fixed_t(hx, hy, t), step_solution[t]) for t in step_solution.keys()],
         label='Fractional steps method error')
plt.legend()

plt.show()
