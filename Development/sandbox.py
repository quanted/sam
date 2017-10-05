import numpy as np

n_scenarios = 3
n_vars = 4
n_dates = 5

a = np.arange(n_scenarios * n_vars * n_dates).reshape((n_scenarios, n_vars, n_dates))

b = a.reshape((n_scenarios, n_vars * n_dates))

c = b.reshape((n_scenarios, n_vars, n_dates))

print(a)
print(b)
print(c)