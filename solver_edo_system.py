import numpy as np

def rk4_system(f, y0, t0, t_end, h):
    n = int((t_end - t0) / h)
    t = np.linspace(t0, t_end, n+1)
    y = np.zeros((n+1, len(y0)))
    y[0] = y0

    for i in range(n):
        k1 = h * f(t[i], y[i])
        k2 = h * f(t[i] + h/2, y[i] + k1/2)
        k3 = h * f(t[i] + h/2, y[i] + k2/2)
        k4 = h * f(t[i] + h, y[i] + k3)
        
        y[i+1] = y[i] + (k1 + 2*k2 + 2*k3 + k4) / 6

    return t, y


### Ejemplo solver sistema edos:

def system(t, y):
    y1, y2 = y
    dydt = np.array([y2, -y1])
    return dydt

# Parámetros
y0 = [0, 1]
t0 = 0
t_end = 2 * np.pi
h = 0.1

# Resolviendo el sistema de EDOs
t, y = rk4_system(system, y0, t0, t_end, h)

# Mostrando los resultados
for i in range(len(t)):
    print(f"t = {t[i]:.2f}, y1 = {y[i,0]:.4f}, y2 = {y[i,1]:.4f}")
