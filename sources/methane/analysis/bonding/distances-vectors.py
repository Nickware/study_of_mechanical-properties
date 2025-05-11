import numpy as np
from numpy.linalg import norm

# Posiciones (C primero, luego H1-H4)
C = np.array([-0.176681, 0.142355, 0.642642])
H1 = np.array([-0.641291, 0.574156, 0.0672497])
H2 = np.array([-0.696224, 0.857166, 0.293665])
H3 = np.array([-0.14708, 0.721901, -0.588574])
H4 = np.array([-0.180149, 1.11447, -0.237455])

# Vectores CH
v1 = H1 - C
v2 = H2 - C
v3 = H3 - C
v4 = H4 - C

# Cálculo de ángulos
def angle(v, w):
    return np.degrees(np.arccos(np.dot(v,w)/(norm(v)*norm(w))))

print(f"Ángulo H1-C-H2: {angle(v1,v2):.2f}°")
print(f"Ángulo H1-C-H3: {angle(v1,v3):.2f}°")
print(f"Ángulo H1-C-H4: {angle(v1,v4):.2f}°")
print(f"Ángulo H2-C-H3: {angle(v2,v3):.2f}°")
print(f"Ángulo H2-C-H4: {angle(v2,v4):.2f}°")
print(f"Ángulo H3-C-H4: {angle(v3,v4):.2f}°")
