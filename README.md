# MES Solver - Documentation and Mathematical Description

Finite Element Method (FEM) solver for **transient heat conduction** analysis in 2D. Supports quadrilateral meshes (DC2D4), multiple materials, convective boundary conditions (external and internal), and Gauss integration.

---

## Table of Contents

1. [Physical Problem Overview](#1-physical-problem-overview)
2. [Weak Formulation and FEM Discretization](#2-weak-formulation-and-fem-discretization)
3. [Reference Element and Shape Functions](#3-reference-element-and-shape-functions)
4. [Numerical Integration - Gauss Quadrature](#4-numerical-integration---gauss-quadrature)
5. [Mapping Jacobian](#5-mapping-jacobian)
6. [Conductivity Matrix H](#6-conductivity-matrix-h)
7. [Capacity Matrix C](#7-capacity-matrix-c)
8. [Convective Boundary Conditions - Hbc and P](#8-convective-boundary-conditions---hbc-and-p)
9. [Time Scheme - Backward Euler Method](#9-time-scheme---backward-euler-method)
10. [Global System and Solution](#10-global-system-and-solution)
11. [Input File Format](#11-input-file-format)
12. [Classes and Structures Description](#12-classes-and-structures-description)

---

## 1. Physical Problem Overview

The equation being solved is the **transient heat conduction equation** (Fourier–Kirchhoff equation):

$$\rho c_p \frac{\partial T}{\partial t} = \nabla \cdot (k \nabla T) + Q$$

where:

| Symbol | Meaning | Unit |
|--------|---------|------|
| $T$ | temperature | °C or K |
| $t$ | time | s |
| $\rho$ | material density | kg/m³ |
| $c_p$ | specific heat | J/(kg·K) |
| $k$ | thermal conductivity | W/(m·K) |
| $Q$ | heat source (not used in solver) | W/m³ |

On the boundaries, the **Newton convective boundary condition** (type III) applies:

$$-k \frac{\partial T}{\partial n} = \alpha (T - T_\infty)$$

where $\alpha$ [W/(m²·K)] is the convection coefficient and $T_\infty$ is the ambient temperature.

---

## 2. Weak Formulation and FEM Discretization

Multiplying the equation by a test function $v$ and integrating by parts (Green's theorem), we obtain the **weak formulation**:

$$\int_\Omega \rho c_p \frac{\partial T}{\partial t} v \, d\Omega + \int_\Omega k \nabla T \cdot \nabla v \, d\Omega + \int_{\Gamma_\alpha} \alpha T v \, d\Gamma = \int_{\Gamma_\alpha} \alpha T_\infty v \, d\Gamma$$

After FEM approximation ($T \approx \mathbf{N}^T \mathbf{T}_e$, where $\mathbf{N}$ is the shape function vector) and global assembly, the problem reduces to a system of ordinary differential equations:

$$\mathbf{C} \dot{\mathbf{T}} + (\mathbf{H} + \mathbf{H}_{bc}) \mathbf{T} = \mathbf{P}$$

---

## 3. Reference Element and Shape Functions

The solver uses **four-node quadrilateral elements** (bilinear). Each real element is mapped to a reference element in the coordinate system $(\xi, \eta) \in [-1,1]^2$.

Shape functions on the reference element:

$$N_1 = \tfrac{1}{4}(1-\xi)(1-\eta), \quad N_2 = \tfrac{1}{4}(1+\xi)(1-\eta)$$

$$N_3 = \tfrac{1}{4}(1+\xi)(1+\eta), \quad N_4 = \tfrac{1}{4}(1-\xi)(1+\eta)$$

Their derivatives with respect to $\xi$ (denoted `dN_dE`):

$$\frac{\partial N_1}{\partial \xi} = -\tfrac{1}{4}(1-\eta), \quad \frac{\partial N_2}{\partial \xi} = +\tfrac{1}{4}(1-\eta), \quad \ldots$$

Derivatives with respect to $\eta$ (denoted `dN_dN`):

$$\frac{\partial N_1}{\partial \eta} = -\tfrac{1}{4}(1-\xi), \quad \frac{\partial N_2}{\partial \eta} = -\tfrac{1}{4}(1+\xi), \quad \ldots$$

Real coordinates of a point in the element:

$$x = \sum_{i=1}^4 N_i \, x_i, \qquad y = \sum_{i=1}^4 N_i \, y_i$$

---

## 4. Numerical Integration - Gauss Quadrature

Integrals over elements are computed using **Gauss–Legendre quadrature**:

$$\int_{-1}^{1} f(\xi) \, d\xi \approx \sum_{i=1}^{N} w_i \, f(\xi_i)$$

The solver supports schemes from $N=1$ to $N=5$ points. For example, for $N=2$:

$$\xi_{1,2} = \pm \frac{1}{\sqrt{3}}, \qquad w_{1,2} = 1$$

For $N=3$:

$$\xi = \left\lbrace -\sqrt{\tfrac{3}{5}},\ 0,\ \sqrt{\tfrac{3}{5}} \right\rbrace, \qquad w = \left\lbrace \tfrac{5}{9},\ \tfrac{8}{9},\ \tfrac{5}{9} \right\rbrace$$

In 2D, integration points are formed as a **Cartesian product** of 1D sets - for $N\times N$ integration points:

$$\iint f(\xi,\eta) \, d\xi \, d\eta \approx \sum_{j=1}^N \sum_{i=1}^N w_i w_j \, f(\xi_i, \eta_j)$$

Default scheme in the solver: `npc = 4` (16 integration points per element).

---

## 5. Mapping Jacobian

The transformation from the local system $(\xi,\eta)$ to the global system $(x,y)$ is described by the **Jacobian matrix**:

$$\mathbf{J} = \begin{bmatrix} \dfrac{\partial x}{\partial \xi} & \dfrac{\partial y}{\partial \xi} \\[8pt] \dfrac{\partial x}{\partial \eta} & \dfrac{\partial y}{\partial \eta} \end{bmatrix} = \begin{bmatrix} \sum_i \dfrac{\partial N_i}{\partial \xi} x_i & \sum_i \dfrac{\partial N_i}{\partial \xi} y_i \\[6pt] \sum_i \dfrac{\partial N_i}{\partial \eta} x_i & \sum_i \dfrac{\partial N_i}{\partial \eta} y_i \end{bmatrix}$$

Jacobian determinant:

$$\det \mathbf{J} = J_{11} J_{22} - J_{12} J_{21}$$

Inverse matrix:

$$\mathbf{J}^{-1} = \frac{1}{\det \mathbf{J}} \begin{bmatrix} J_{22} & -J_{12} \\ -J_{21} & J_{11} \end{bmatrix}$$

Global derivatives of shape functions (required for matrix H):

$$\frac{\partial N_i}{\partial x} = J^{-1}_{11} \frac{\partial N_i}{\partial \xi} + J^{-1}_{12} \frac{\partial N_i}{\partial \eta}$$

$$\frac{\partial N_i}{\partial y} = J^{-1}_{21} \frac{\partial N_i}{\partial \xi} + J^{-1}_{22} \frac{\partial N_i}{\partial \eta}$$

---

## 6. Conductivity Matrix H

**Local conductivity matrix** for an element:

$$H_{ij}^{(e)} = \int_{\Omega^{(e)}} k \left( \frac{\partial N_i}{\partial x} \frac{\partial N_j}{\partial x} + \frac{\partial N_i}{\partial y} \frac{\partial N_j}{\partial y} \right) d\Omega$$

After substituting the Gauss quadrature in the $(\xi,\eta)$ system:

$$H_{ij}^{(e)} \approx k \sum_{p=1}^{N_{pc}^2} w_{p_x} w_{p_y} \left( \frac{\partial N_i}{\partial x}\bigg|_p \frac{\partial N_j}{\partial x}\bigg|_p + \frac{\partial N_i}{\partial y}\bigg|_p \frac{\partial N_j}{\partial y}\bigg|_p \right) \det\mathbf{J}_p$$

---

## 7. Capacity Matrix C

The solver first computes the **consistent capacity matrix** (consistent mass matrix):

$$C_{ij}^{\text{cons}} = \rho c_p \int_{\Omega^{(e)}} N_i N_j \, d\Omega \approx \rho c_p \sum_p w_{p_x} w_{p_y} N_i(\xi_p, \eta_p) N_j(\xi_p, \eta_p) \det\mathbf{J}_p$$

and then applies the **lumped capacity matrix** (lumped mass matrix) by row summation:

$$C_{ii}^{\text{lump}} = \sum_{j=1}^{4} C_{ij}^{\text{cons}}, \qquad C_{ij}^{\text{lump}} = 0 \text{ for } i \neq j$$

The lumped matrix is diagonal - it simplifies solving and improves numerical stability for large time steps.

---

## 8. Convective Boundary Conditions - Hbc and P

For element edges with convection (BC = 1 or BC = 2), a **boundary matrix** is added:

$$H_{bc,ij}^{(e)} = \int_{\Gamma_\alpha} \alpha N_i N_j \, d\Gamma \approx \alpha \sum_{p=1}^{N_{pc}} w_p N_i(s_p) N_j(s_p) \det\mathbf{J}_{\text{edge}}$$

where $\det\mathbf{J}_{\text{edge}} = \frac{L_\text{edge}}{2}$ ($L$ is the edge length), and $s_p$ are Gauss points along the edge.

The solver handles four element edges: bottom (0–1), right (1–2), top (2–3), and left (3–0), with shape functions interpolated along one variable:

$$N_1^{\text{edge}} = \tfrac{1}{2}(1-s), \qquad N_2^{\text{edge}} = \tfrac{1}{2}(1+s)$$

**Convective load vector**:

$$P_i^{(e)} = \int_{\Gamma_\alpha} \alpha T_\infty N_i \, d\Gamma \approx \alpha T_\infty \sum_{p} w_p N_i(s_p) \det\mathbf{J}_{\text{edge}}$$

---

## 9. Time Scheme - Backward Euler Method

The system of equations:

$$\mathbf{C} \dot{\mathbf{T}} + (\mathbf{H} + \mathbf{H}_{bc}) \mathbf{T} = \mathbf{P}$$

is discretized in time using the **backward Euler method** (implicit):

$$\frac{\mathbf{C}}{\Delta t} (\mathbf{T}^{n+1} - \mathbf{T}^n) + (\mathbf{H} + \mathbf{H}_{bc}) \mathbf{T}^{n+1} = \mathbf{P}$$

After rearranging:

$$\underbrace{\left(\frac{\mathbf{C}}{\Delta t} + \mathbf{H} + \mathbf{H}_{bc}\right)}_{\mathbf{H}_{\text{eff}}} \mathbf{T}^{n+1} = \underbrace{\mathbf{P} + \frac{\mathbf{C}}{\Delta t} \mathbf{T}^n}_{\mathbf{P}_{\text{eff}}}$$

The backward Euler scheme is **unconditionally stable** (any $\Delta t$), though only first-order accurate.

---

## 10. Global System and Solution

Local matrices are assembled into global ones using the direct assembly method:

$$H_{gi, gj}^{\text{glob}} \mathrel{+}= H_{ij}^{(e)}, \quad \text{where } g_i = \text{ID}^{(e)}_i - 1$$

The linear system $\mathbf{H}_\text{eff} \mathbf{T} = \mathbf{P}_\text{eff}$ is solved using **Gaussian elimination with partial pivoting**:

```
for i = 0..N-1:
    find row maxRow with max |A[k][i]|
    swap rows i ↔ maxRow
    forward elimination (zeroing column i)
back substitution
```

---

## 11. Input File Format

A text file with sections separated by keywords:

```
SimulationTime     [s]         # total simulation time
SimulationStepTime [s]         # time step
AlfaOut            [W/m²·K]    # external convection coefficient
AlfaIn             [W/m²·K]    # internal convection coefficient
TotOut             [°C]        # external ambient temperature
TotIn              [°C]        # internal ambient temperature
InitialTemp        [°C]        # initial temperature
Nodes number       [-]
Elements number    [-]

*Node
  ID, x, y

*Element, type=DC2D4
  ID, N1, N2, N3, N4          # nodes in CCW order

*Materials
  ID, name, k[W/mK], rho[kg/m³], cp[J/kgK]

*ElementMat
  elemID, matID               # material assignment to element

*BC                            # boundary nodes (single alpha type)
  n1, n2, ...

*BC_OUT                        # external boundary nodes (AlfaOut, TotOut)
  n1, n2, ...

*BC_IN                         # internal boundary nodes (AlfaIn, TotIn)
  n1, n2, ...
```

**Note:** Element nodes must be given in **counter-clockwise (CCW) order**. The solver automatically detects and corrects wrong orientation.

---

## 12. Classes and Structures Description

### `GaussQuadrature`

Stores Gauss–Legendre quadrature points and weights for $N \in \{1,...,5\}$.

| Field | Type | Description |
|-------|------|-------------|
| `N` | `int` | number of points |
| `xi` | `vector<double>` | point coordinates ($\xi_i$) |
| `w` | `vector<double>` | weights ($w_i$) |

---

### `Node`

A single mesh node.

| Field | Type | Description |
|-------|------|-------------|
| `id` | `int` | node number (from 1) |
| `x`, `y` | `double` | global coordinates |
| `BC` | `int` | boundary condition type: 0=none, 1=OUT, 2=IN |

---

### `ElemUniv`

Reference element - pre-calculation of shape function values and their derivatives at all integration points.

| Field | Type | Description |
|-------|------|-------------|
| `npc` | `int` | number of Gauss points in one direction |
| `gauss` | `GaussQuadrature` | integration scheme |
| `N[p][i]` | `vector<vector<double>>` | $N_i$ at point $p$ |
| `dN_dE[p][i]` | `vector<vector<double>>` | $\partial N_i / \partial \xi$ at point $p$ |
| `dN_dN[p][i]` | `vector<vector<double>>` | $\partial N_i / \partial \eta$ at point $p$ |

Integration points are numbered row-by-row: `p = j * npc + i` for (`i`-th in $\xi$, `j`-th in $\eta$).

---

### `Surface`

Stores shape functions along the four element edges, used for integrating boundary conditions.

| Field | Type | Description |
|-------|------|-------------|
| `npc` | `int` | number of Gauss points along the edge |
| `N[e][p][i]` | 3D `vector` | value of $N_i$ on edge `e` at point `p` |

Edges (local node indices):

| `e` | Edge | Nodes |
|-----|------|-------|
| 0 | bottom | 0–1 |
| 1 | right | 1–2 |
| 2 | top | 2–3 |
| 3 | left | 3–0 |

---

### `Jakobian`

Computes and stores the mapping Jacobian for a single integration point.

| Field | Type | Description |
|-------|------|-------------|
| `J[2][2]` | `double[2][2]` | Jacobian matrix $\mathbf{J}$ |
| `Jodwr[2][2]` | `double[2][2]` | inverse matrix $\mathbf{J}^{-1}$ |
| `detJ` | `double` | determinant $\det\mathbf{J}$ |
| `dN_dx[4]`, `dN_dy[4]` | `double[4]` | global derivatives $\partial N_i/\partial x$, $\partial N_i/\partial y$ |

Constructor `Jakobian(elem, grid, eUniv, p)` automatically computes all values above for point $p$.

---

### `Material`

Material parameters for one layer.

| Field | Type | Description |
|-------|------|-------------|
| `id` | `int` | identifier |
| `name` | `string` | name (e.g. `"plaster"`, `"brick"`) |
| `conductivity` | `double` | $k$ [W/(m·K)] |
| `density` | `double` | $\rho$ [kg/m³] |
| `specificHeat` | `double` | $c_p$ [J/(kg·K)] |

---

### `Element`

Four-node finite element. Performs all local FEM arithmetic.

| Field | Type | Description |
|-------|------|-------------|
| `id` | `int` | element number |
| `ID[4]` | `int[4]` | node numbers (global, from 1) |
| `jakobiany` | `vector<Jakobian>` | Jacobians for each integration point |
| `H[4][4]` | `vector<vector<double>>` | local conductivity matrix |
| `C[4][4]` | `vector<vector<double>>` | local capacity matrix (lumped) |
| `Hbc[4][4]` | `vector<vector<double>>` | local convection boundary matrix |
| `P[4]` | `vector<double>` | local convective load vector |
| `materialId` | `int` | ID of assigned material |

Methods:

| Method | Description |
|--------|-------------|
| `obliczH(eUniv, grid, mat)` | Builds $\mathbf{H}^{(e)}$ and $\mathbf{C}^{(e)}$ for the given material |
| `obliczHbc(surface, grid, alfaOut, alfaIn)` | Builds $\mathbf{H}_{bc}^{(e)}$ for active boundary edges |
| `obliczP(surface, grid, alfaOut, alfaIn, TotOut, TotIn)` | Builds $\mathbf{P}^{(e)}$ for active boundary edges |

---

### `Grid`

Container for the entire computational mesh.

| Field | Type | Description |
|-------|------|-------------|
| `nN` | `int` | number of nodes |
| `nE` | `int` | number of elements |
| `nodes` | `vector<Node>` | list of nodes |
| `elements` | `vector<Element>` | list of elements |
| `boundaryNodes` | `vector<int>` | boundary nodes |

---

### `GlobalData`

Global simulation parameters read from file.

| Field | Type | Description |
|-------|------|-------------|
| `SimulationTime` | `double` | end time [s] |
| `SimulationStepTime` | `double` | time step $\Delta t$ [s] |
| `AlfaOut` | `double` | external $\alpha$ [W/(m²·K)] |
| `AlfaIn` | `double` | internal $\alpha$ [W/(m²·K)] |
| `TotOut` | `double` | external $T_\infty$ [°C] |
| `TotIn` | `double` | internal $T_\infty$ [°C] |
| `InitialTemp` | `double` | initial temperature [°C] |
| `nN`, `nE` | `int` | number of nodes and elements |
| `materials` | `vector<Material>` | list of materials |
| `materialIdToIndex` | `map<int,size_t>` | ID → index mapping in vector |
| `npc` | `int` | number of Gauss points (default 4) |

---

### Helper Functions

| Function | Description |
|----------|-------------|
| `loadFromFile(filename, globalData, grid)` | Input file parser; handles all sections `*Node`, `*Element`, `*BC`, `*Materials`, etc. |
| `solveLinearSystem(A, b)` | Gaussian elimination with partial pivoting |
| `isCCW(x, y)` | Checks CCW orientation of a quadrilateral (signed area algorithm); if CW - swaps nodes 1↔3 |
| `gauss1D(f, G)` | 1D integration of function $f$ using Gauss scheme `G` |
| `gauss2D(f, G)` | 2D integration of function $f(\xi,\eta)$ using scheme `G`×`G` |
| `printMatrix(M, title)` | Prints matrix to `cout` |
| `printDerivativeTable(deriv, label)` | Prints shape function derivative table |

---

## Computation Flow Diagram

```
loadFromFile()
    │
    ├── read GlobalData (parameters, materials)
    └── read Grid (nodes, elements, BC)

for each element:
    ├── compute Jacobians (npc² points)
    ├── obliczH()    →  H[4×4], C[4×4]
    ├── obliczHbc()  →  Hbc[4×4]
    └── obliczP()    →  P[4]
    │
    └── assemble into global matrices:
        Hglobal, Cglobal, Hbcglobal, Pglobal

time loop (steps = SimulationTime / SimulationStepTime):
    │
    ├── H_eff = H + Hbc + C/dt
    ├── P_eff = P + (C/dt)·T^n
    └── T^{n+1} = solveLinearSystem(H_eff, P_eff)

print Tmin/Tmax per material per step
```

---

## Example Input Files

### Small test (4×4 nodes, 9 elements, single material)

```
SimulationTime 500
SimulationStepTime 50
Conductivity 25
Alfa 300
Tot 1200
InitialTemp 100
Nodes number 16
Elements number 9
...
*BC
1, 2, 3, 4, 5, 8, 9, 12, 13, 14, 15, 16
```

### Multi-material wall (plaster + styrofoam + brick)

```
SimulationTime 86400
SimulationStepTime 3600
AlfaIn 8
AlfaOut 25
TotOut 0
TotIn 23
*Materials
1, plaster,    1.0,  1800, 1000
2, styrofoam,  0.040,  40, 1460
3, brick,      0.6,  1600,  840
*BC_OUT
1, 12, 23, ...
*BC_IN
11, 22, 33, ...
```

---

## Compilation Requirements

```bash
g++ -O2 -std=c++17 -o mes_solver main.cpp
./mes_solver
```

The input file is selected in `main()` via the `Pliki` vector index. Change `Pliki[5]` to the desired file.
