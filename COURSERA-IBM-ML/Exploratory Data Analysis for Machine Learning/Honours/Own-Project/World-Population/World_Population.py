import matplotlib.pyplot as plt

a = [10 ,20 ,30, 40]
b = [1, 2, 3, 4]

fig, ax = plt.subplots(1, 1)
ax.scatter(a, b)
fig.savefig("sample.png")