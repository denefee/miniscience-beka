from dolfin import *
from mshr import *

# Create mesh and define function space
domain = Circle(Point(0, 0), 1)
mesh = generate_mesh(domain, 64)
# mesh = UnitSquareMesh(32, 32)
V = FunctionSpace(mesh, "Lagrange", 1)

# Define Dirichlet boundary (x = 0 or x = 1)
def boundary(x):
    return (x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS)and(x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS)
# tol = 1E-14
# def boundary(x):
#     return abs(x[0]) < tol or abs(x[1]) < tol \
#         or abs(x[0] - 1) < tol or abs(x[1] - 1) < tol

# Define boundary condition
u0 = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)
# u0 = Constant(0.0)
bc = DirichletBC(V, u0, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
# f = Expression("sin(40*x[0]) + sin(40*x[1])", degree=2)
f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)c
# g = Expression("sin(5*x[0])", degree=2)
a = inner(grad(u), grad(v))*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bc)

# Save solution in VTK format
file = File("poisson/poisson.pvd")
file << u

# Plot solution
import matplotlib.pyplot as plt
plot(u)
plt.show()