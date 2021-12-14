import matplotlib.pyplot as plt
import numpy as np
import subprocess


def analytic_solution(x, y):
    return np.cos(x) * np.cos(y)


def analytic_solution_fixed_y(y, hx, count=None):
    res = np.array([analytic_solution(x, y) for x in np.arange(0, np.pi / 2, hx)])
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


def run_prog(hx: float, hy: float, eps: float):
    input = str(hx) + ' ' + str(hy) + ' ' + str(eps)
    return subprocess.run(['cmake-build-debug/lab7'], shell=True, input=input.encode('utf-8')).stdout


def read_file(file: str, hy: float) -> dict:
    res = dict()
    count = 0
    with open(file, 'r+') as f:
        line = f.readline()
        while line:
            res[hy * count] = np.array(list(map(float, line.split())))
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


hx = 0.01
hy = 0.01
eps = 0.000001

run_prog(hx, hy, eps)
lieb_solution = read_file('lieb.txt', hy)
seid_solution = read_file('seid.txt', hy)
relax_solution = read_file('relax.txt', hy)

x = np.arange(0, np.pi / 2, hx)

y = get_closest_to(lieb_solution, 0.3)

fig = plt.figure()
fig.suptitle(f'Solution at y = {y}')
plt.plot(x, lieb_solution[y], label='Liebmann method')
plt.plot(x, seid_solution[y], label='Seidel method')
plt.plot(x, relax_solution[y], label='Relaxation method, tau = 0.5')
plt.plot(x, analytic_solution_fixed_y(y, hx), label='Analytic solution')
plt.legend()

y = get_closest_to(lieb_solution, 0.6)

fig = plt.figure()
fig.suptitle(f'Solution at y = {y}')
plt.plot(x, lieb_solution[y], label='Liebmann method')
plt.plot(x, seid_solution[y], label='Seidel method')
plt.plot(x, relax_solution[y], label='Relaxation method, tau = 0.5')
plt.plot(x, analytic_solution_fixed_y(y, hx), label='Analytic solution')
plt.legend()

y = get_closest_to(lieb_solution, 0.9)

fig = plt.figure()
fig.suptitle(f'Solution at y = {y}')
plt.plot(x, lieb_solution[y], label='Liebmann method')
plt.plot(x, seid_solution[y], label='Seidel method')
plt.plot(x, relax_solution[y], label='Relaxation method, tau = 0.5')
plt.plot(x, analytic_solution_fixed_y(y, hx), label='Analytic solution')
plt.legend()

fig = plt.figure()
fig.suptitle(f'Eps - iteration number dependency, target eps = {eps}')
lieb_eps = read_eps_file('lieb_eps.txt')
seid_eps = read_eps_file('seid_eps.txt')
relax_eps = read_eps_file('relax_eps.txt')
plt.plot(lieb_eps, label='Liebmann method')
plt.plot(seid_eps, label='Seidel method')
plt.plot(relax_eps, label='Relaxation method, tau = 0.5')

plt.yscale('log')

plt.legend()

plt.show()
