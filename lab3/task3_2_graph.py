import matplotlib.pyplot as plt

def func(a, b, c, d, x, x_i):
    return a + b * (x - x_i) + c * (x - x_i) ** 2 + d * (x - x_i) ** 3

x = []
y = []

with open("task3_2_plot.txt") as f:
    line = f.readline()
    while line:
        values = line.split()
        assert len(values) == 6
        left, right = float(values[0]), float(values[1])
        a, b, c, d = float(values[2]), float(values[3]), float(values[4]), float(values[5])
        i = left
        while i <= right:
            x.append(i)
            y.append(func(a, b, c, d, i, left))
            i += 0.01
        line = f.readline()    


plt.plot(x, y)
plt.show()
