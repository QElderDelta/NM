import matplotlib.pyplot as plt

x = []
y = []
first_order_curve_x = []
first_order_curve_y = []
first_order_curve_coefs = []
second_order_curve_x = []
second_order_curve_y = []
second_order_curve_coefs = []

def func(coefs, point):
    result = 0.0
    for i in range(len(coefs)):
        result += coefs[i] * point ** i
    return result


with open("task3_3.txt") as f:
    line = f.readline()
    count = int(line)
    for i in range(count):
        line = f.readline()
        x.append(float(line))
    for i in range(count):
        line = f.readline()
        y.append(float(line))

with open("task3_3_plot.txt") as f:
    line = f.readline()
    coefs = line.split()
    for coef in coefs:
        first_order_curve_coefs.append(float(coef))
    line = f.readline()
    coefs = line.split()    
    for coef in coefs:
        second_order_curve_coefs.append(float(coef))

plt.plot(x, y, 'o', color='black')
i = x[0]
while i <= x[len(x) - 1]:
    first_order_curve_x.append(i)
    second_order_curve_x.append(i)
    first_order_curve_y.append(func(first_order_curve_coefs, i))
    second_order_curve_y.append(func(second_order_curve_coefs, i))
    i += 0.01
plt.plot(first_order_curve_x, first_order_curve_y, linestyle='solid', label='First order curve')
plt.plot(second_order_curve_x, second_order_curve_y, linestyle='dashed', label='Second order curve')
plt.legend()
plt.show()