## Universidad del Valle de Guatemala
## Modelación y Simulación
## Laboratorio #3
## Algoritmos Rungue-Kutta para EDO y sistemas de EDO

import numpy as np

def rk4_single(f, y0, t0, t_end, h):
    n = int((t_end - t0) / h)
    t = np.linspace(t0, t_end, n+1)
    y = np.zeros(n+1)
    y[0] = y0

    for i in range(n):
        k1 = h * f(t[i], y[i])
        k2 = h * f(t[i] + h/2, y[i] + k1/2)
        k3 = h * f(t[i] + h/2, y[i] + k2/2)
        k4 = h * f(t[i] + h, y[i] + k3)
        
        y[i+1] = y[i] + (k1 + 2*k2 + 2*k3 + k4) / 6

    return t, y


### Ejemplo solucionando solo una EDO
def dydt(t, y):
    return -2 * y

# Parámetros
y0 = 1
t0 = 0
t_end = 2
h = 0.1

# Resolviendo la EDO
t, y = rk4_single(dydt, y0, t0, t_end, h)

# Mostrando los resultados
for i in range(len(t)):
    print(f"t = {t[i]:.2f}, y = {y[i]:.4f}")


