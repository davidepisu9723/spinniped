# Spinniped

<p align="center">
  <img src="https://raw.githubusercontent.com/davidepisu9723/spinniped/master/logo/spinniped-logo.png" width="240" alt="Spinniped logo">
</p>

Spinniped is a small finite-element toolbox for deterministic and stochastic
rotor dynamics. Its input API follows a NASTRAN-like structure: grids,
coordinate systems, materials, properties, and elements are immutable Python
dataclasses whose fields contain ordinary Python values. NumPy arrays appear
only after a definition is built into a numerical model.

The project is under active development and is intended for preliminary design,
research, and education.

## Installation

Spinniped requires Python 3.13 or newer, NumPy, and SciPy.

Until a release is available on PyPI, install a local clone in editable mode:

```bash
git clone https://github.com/davidepisu9723/spinniped.git
cd spinniped
python -m pip install -e .
```

Install the test dependency with:

```bash
python -m pip install -e ".[test]"
```

## Model-definition API

A model definition is data, not a collection of stateful finite-element
objects. Records refer to one another by integer IDs:

- `Grid` defines a point and its input coordinate system.
- `CoordinateSystem` defines a local Cartesian frame.
- `Material` defines density, Young's modulus, and Poisson's ratio.
- `ShaftProperty`, `BearingProperty`, and `DiskProperty` define reusable
  physical data.
- `ShaftElement`, `BearingElement`, and `DiskElement` define connectivity and
  reference a compatible property.
- `ModelDefinition` holds the record lists.

`LumpedMassElement` and `LumpedMassProperty` are aliases for the corresponding
disk records.

This example defines a three-grid shaft with a disk and two bearings:

```python
from spinniped import (
    BearingElement,
    BearingProperty,
    DiskElement,
    DiskProperty,
    Grid,
    Material,
    ModelDefinition,
    RandomDistribution,
    ShaftElement,
    ShaftProperty,
)

definition = ModelDefinition(
    grids=[
        Grid(1, z=0.0),
        Grid(2, z=0.5),
        Grid(3, z=1.0),
    ],
    materials=[
        Material(
            id=1,
            density=7850.0,
            young_modulus=(1,),
            poisson_ratio=0.3,
        ),
    ],
    properties=[
        ShaftProperty(
            id=1,
            material=1,
            outer_diameter=(2,),
            theory="timoshenko",
        ),
        BearingProperty(
            id=2,
            kxx=1e8,
            kyy=1e8,
            cxx=100.0,
            cyy=100.0,
        ),
        DiskProperty(
            id=3,
            mass=2.0,
            diametral_inertia=0.01,
            polar_inertia=0.02,
        ),
    ],
    elements=[
        ShaftElement(1, grid_a=1, grid_b=2, property=1),
        ShaftElement(2, grid_a=2, grid_b=3, property=1),
        BearingElement(3, grid=1, property=2),
        BearingElement(4, grid=3, property=2),
        DiskElement(5, grid=2, property=3),
    ],
    distributions=[
        RandomDistribution(
            id=1,
            name="steel Young modulus",
            distribution="normal",
            parameters={"mean": 210e9, "stdv": 5e9},
        ),
        RandomDistribution(
            id=2,
            name="shaft diameter",
            distribution="uniform",
            parameters={"low": 0.019, "high": 0.021},
        ),
    ],
    spin_axis=(0.0, 0.0, 1.0),
)
```

### Record reference

All dimensions and coefficients must use one coherent unit system; examples
and theory notes use SI.

| Record | Fields |
|---|---|
| `Grid` | `id`; coordinates `x`, `y`, `z` (default `0.0`); `coordinate_system` (default global ID `0`) |
| `CoordinateSystem` | `id`; `origin`; direction `x_axis`; direction `xy_plane` |
| `Material` | `id`; `density`; `young_modulus`; `poisson_ratio` |
| `ShaftProperty` | `id`; material reference `material`; `outer_diameter`; `inner_diameter` (default `0.0`); `theory`; `rotary_inertia`; `damping` |
| `BearingProperty` | `id`; stiffnesses `kxx`, `kyy`, `kzz`, `kxy`, `kyx`; damping coefficients `cxx`, `cyy`, `czz`, `cxy`, `cyx` |
| `DiskProperty` | `id`; `mass`; `diametral_inertia`; `polar_inertia`; `damping` |
| `ShaftElement` | `id`; end references `grid_a`, `grid_b`; shaft `property` reference |
| `BearingElement` | `id`; first `grid`; bearing `property`; optional `grid_b`; property `coordinate_system` |
| `DiskElement` | `id`; nodal `grid`; disk `property`; property `coordinate_system` |
| `RandomDistribution` | unique `id`; unique `name`; `distribution`; family-specific `parameters` |

`ShaftProperty.theory` accepts `"timoshenko"` (the default) and `"euler"`.
Bending rotary inertia is enabled by default. Shaft and disk `damping` values
are coefficients in the mass-proportional relation
$\mathbf{C}=\alpha\mathbf{M}$. Bearing coefficients default to zero and occupy
the translational x-y-z block.

`ModelDefinition` groups the `grids`, `coordinate_systems`, `materials`,
`properties`, `elements`, and `distributions` lists; each defaults to an empty
list. Its `spin_axis` field is the global spin direction and defaults to
`(0.0, 0.0, 1.0)`.

IDs must be unique positive integers within each record category. A shaft
element must reference a `ShaftProperty`, a bearing element a
`BearingProperty`, and a disk element a `DiskProperty`. A shaft property also
references a material ID. The builder validates these relationships.

### Coordinate systems

`Grid.x`, `Grid.y`, and `Grid.z` are interpreted in the coordinate system
selected by `Grid.coordinate_system`. System `0` is reserved for the global
Cartesian system and must not be redefined. A custom system is defined by a
global origin, a direction for its positive x-axis, and another direction in
its x-y plane. The two direction vectors need not be normalized, but they must
be nonzero and non-collinear:

```python
from spinniped import CoordinateSystem, Grid

frame = CoordinateSystem(
    id=10,
    origin=(1.0, 0.0, 0.0),
    x_axis=(0.0, 1.0, 0.0),
    xy_plane=(0.0, 0.0, 1.0),
)

point = Grid(20, x=0.5, coordinate_system=10)
```

Shaft matrices are formulated in an element-local frame whose z-axis runs
from `grid_a` to `grid_b`, then transformed into the global frame during
assembly. Shaft elements therefore need not be parallel to the global z-axis.

`BearingElement.coordinate_system` rotates bearing coefficients from the
selected local frame into the global frame. With `grid_b=None` (the default),
the bearing connects `grid` to ground. Supplying `grid_b` creates a relative
two-grid bearing with the block pattern

$$
\begin{bmatrix}
\mathbf{A}_b & -\mathbf{A}_b\\
-\mathbf{A}_b & \mathbf{A}_b
\end{bmatrix},
$$

for both bearing stiffness and damping. `DiskElement.coordinate_system`
similarly orients disk inertia and its local spin axis. These element
coordinate-system fields default to global system `0`. A grid's own
`coordinate_system` locates the grid; it does not implicitly orient an attached
bearing or disk.

### Spin-axis convention

`ModelDefinition.spin_axis` is a nonzero direction expressed in global
coordinates. The builder normalizes it, so it defines orientation rather than
speed. For each shaft, the unit-speed gyroscopic matrix is multiplied by the
projection of this direction onto the shaft's local axis. Reversing
`grid_a`/`grid_b` therefore does not reverse the physical global gyroscopic
effect. A shaft perpendicular to the spin axis has zero gyroscopic
contribution in the current formulation.

For a disk, the same projection uses the local z-axis selected by
`DiskElement.coordinate_system`. The scalar `speed` or each Campbell `speeds`
value later supplies the angular-speed magnitude in rad/s.

## Random parameters

Random distributions are registered once in ``ModelDefinition.distributions``.
Each distribution has a unique integer ``id``, a unique descriptive ``name``,
a family, and its parameters:

```python
RandomDistribution(
    id=17,
    name="bearing stiffness",
    distribution="normal",
    parameters={"mean": 1.0e8, "stdv": 1.0e6},
)
```

Any numeric model parameter can reference this scalar distribution with a
one-item tuple. Reusing the reference reuses exactly the same draw:

```python
BearingProperty(id=1, kxx=(17,), kyy=(17,))
```

Correlated variables use one multivariate distribution and select components
with ``(distribution_id, component_index)``:

```python
RandomDistribution(
    id=18,
    name="correlated bearing stiffness",
    distribution="multivariate_normal",
    parameters={
        "mean": [1.0e8, 1.2e8],
        "stdv": [1.0e7, 2.0e7],
        "correlation": [[1.0, 0.8], [0.8, 1.0]],
    },
)

BearingProperty(id=1, kxx=(18, 0), kyy=(18, 1))
```

The supported families are ``"normal"``, ``"uniform"``, and
``"multivariate_normal"``. Schemas are strict, all parameters must be finite,
standard deviations must be nonnegative, and correlation matrices must be
symmetric and positive semidefinite with a unit diagonal. Distribution IDs and
names must both be unique.

Resolved records and numerical inputs are validated for type, shape, finite
values, references, and physical bounds. For example, shaft length, density,
and Young's modulus must be positive; `outer_diameter > inner_diameter >= 0`;
Poisson's ratio must lie strictly between `-1` and `0.5`; and disk mass and
inertias cannot be negative. A stochastic draw that violates a physical bound
raises an error identifying its sample and field path; the builder does not
silently clip or resample it.

### Deterministic build

The default build is deterministic. Normal distributions resolve to their
means and uniform distributions to their interval midpoints. Exactly one set of
global matrices is assembled:

```python
from spinniped import ModelBuilder

model = ModelBuilder().build(definition)

# Equivalent, explicit form:
model = ModelBuilder().build(definition, stochastic=False)
```

### Stochastic build

A stochastic build draws the requested number of Monte Carlo realizations and
assembles one matrix set per realization:

```python
ensemble = ModelBuilder().build(
    definition,
    stochastic=True,
    samples=1_000,
    seed=42,
)
```

The seed initializes NumPy's random generator, making the sampled definitions
and assembled matrices reproducible. `samples` must be a positive integer.
The convenience function `build_model(definition, **options)` has the same
behavior as `ModelBuilder().build(...)`.

## Built models and matrix conventions

`ModelBuilder.build` returns a `BuiltModel`. Deterministic and stochastic builds
use the same layout:

| Attribute | Shape or type | Meaning |
|---|---|---|
| `K`, `stiffness` | `(samples, ndof, ndof)` | Global stiffness matrices |
| `M`, `mass` | `(samples, ndof, ndof)` | Global mass matrices |
| `C`, `damping` | `(samples, ndof, ndof)` | Global viscous damping matrices |
| `G`, `gyroscopic` | `(samples, ndof, ndof)` | Global unit-speed gyroscopic matrices |
| `coordinates` | `(samples, grids, 3)` | Resolved global grid coordinates |
| `definitions` | tuple | Fully resolved deterministic definitions |
| `grid_ids` | tuple | User grid IDs in assembly order |
| `samples` | integer | Number of realizations |
| `ndof` | integer | Number of global degrees of freedom |
| `stochastic` | Boolean | Whether random sampling was used |
| `seed` | integer or `None` | Requested random seed |

Even a deterministic model has a leading sample dimension of length one. For
example, its stiffness matrix is `model.K[0]`, not `model.K`.

Every grid has six degrees of freedom in this order:

```text
x, y, z, tx, ty, tz
```

If a grid has compact assembly index `p`, its global DOFs are `6*p` through
`6*p + 5`. User grid IDs may be non-consecutive; use `model.grid_ids` to see
their assembly order.

The matrices implement the speed-dependent equation

$$
\mathbf{M}\ddot{\mathbf{q}}+
(\mathbf{C}+\Omega\mathbf{G})\dot{\mathbf{q}}+
\mathbf{K}\mathbf{q}=\mathbf{f}.
$$

`G` is a unit-speed matrix and is multiplied by angular speed `Omega` in the
Campbell, frequency-response, and time-response analyses. Speeds are expressed
in radians per second.

Matrix contributions currently are:

| Element | Stiffness | Mass | Damping | Gyroscopic |
|---|---|---|---|---|
| Shaft | Elastic shaft matrix | Consistent shaft matrix | Mass-proportional | Consistent lateral rotary-inertia matrix, projected onto `spin_axis` |
| Bearing | Local translational coefficients, to ground or between two grids | Zero | Local translational coefficients, to ground or between two grids | Zero |
| Disk/lumped mass | Zero | Oriented concentrated mass and inertia | Mass-proportional | Oriented polar inertia, projected onto `spin_axis` |

The low-level local-matrix functions live in `spinniped.stiffness`,
`spinniped.mass`, `spinniped.damping`, and `spinniped.gyroscopic`. Normal model
construction should go through `ModelBuilder`, which applies coordinate
transformations and global assembly consistently.

## Solving a model

`Solver.solve` is the single analysis entry point. Supply the built model to the
constructor:

```python
from spinniped import Solver

fixed = [2, 5]  # Remove axial and torsional rigid motions.
solver = Solver(model)
result = solver.solve("modal", fixed_dofs=fixed, modes=8)
```

or to an individual call:

```python
solver = Solver()
result = solver.solve(
    "modal",
    model=model,
    fixed_dofs=fixed,
    modes=8,
)
```

`fixed_dofs` contains global DOF indices. The solver removes duplicate indices,
validates their range, and solves on the remaining DOFs. Every result is a
dictionary containing `analysis`, `free_dofs`, `samples`, and `stochastic` in
addition to analysis-specific values. Solver inputs must be finite and have the
documented shapes; unknown options are rejected. Results are regular NumPy
arrays, never ragged object arrays. If an analysis cannot return a consistent
mode count, it raises an error and asks for fewer or a fixed number of modes.

### Modal analysis

The modal solver handles the real, undamped generalized eigenproblem

$$
\mathbf{K}\boldsymbol{\phi}=
\lambda\mathbf{M}\boldsymbol{\phi}.
$$

```python
result = solver.solve(
    "modal",
    fixed_dofs=fixed,
    modes=8,
    track_modes=True,
)
```

Result values include:

- `eigenvalues`: `(samples, modes)` values of $\omega^2$ in $\mathrm{s}^{-2}$;
- `frequencies`: `(samples, modes)` frequencies in Hz;
- `eigenvectors`: `(samples, free_dofs, modes)` reduced mode shapes.

Eigenvalues within `zero_tolerance` of zero are treated as rigid modes and
removed. The default tolerance is scaled from the eigenspectrum; it may be
overridden with a finite, nonnegative `zero_tolerance`. An eigenvalue below the
negative tolerance raises an instability or insufficient-constraint error
rather than being silently discarded.

The real modal route requires symmetric stiffness and mass matrices; use the
Campbell route when cross-coupled bearing coefficients make stiffness
non-symmetric. Damping and gyroscopic terms are not included in modal analysis.
The requested `modes` count is strict: requesting more positive modes than are
available raises an error.

For stochastic models, `track_modes=True` (the default) matches every later
sample to sample zero by maximizing the Modal Assurance Criterion (MAC), then
aligns eigenvector signs. Set `track_modes=False` to retain independent
ascending-eigenvalue order in every sample.

### Campbell analysis

```python
result = solver.solve(
    "campbell",
    fixed_dofs=fixed,
    speeds=[0.0, 100.0, 200.0],  # rad/s
    modes=8,
    track_modes=True,
)
```

The solver evaluates the damped first-order state-space eigenproblem at each
speed. Its analysis-specific results are:

- `speeds`: the requested angular speeds in rad/s, shape `(speeds,)`;
- `speeds_hz`: those speeds divided by $2\pi$, shape `(speeds,)`;
- `eigenvalues`: complex state eigenvalues, shape
  `(samples, speeds, modes)`;
- `frequencies`: absolute imaginary parts in Hz, shape
  `(samples, speeds, modes)`;
- `eigenvectors`: complex state vectors, shape
  `(samples, speeds, 2*free_dofs, modes)`.

With `track_modes=True`, roots after the first speed are assigned to the
previous speed by displacement-vector MAC and phase-aligned, forming continuous
modal branches through crossings. Stochastic Campbell results are then matched
at each speed to sample zero. Set `track_modes=False` to sort each speed and
sample independently by positive imaginary part. The requested `modes` count is
strict at every speed.

### Frequency response

Input frequencies are in Hz. The force may contain either all global DOFs or
only the free DOFs:

```python
force = [0.0] * model.ndof
force[6] = 1.0

result = solver.solve(
    "frequency_response",
    fixed_dofs=fixed,
    frequencies=[10.0, 20.0, 30.0],
    force=force,
    speed=100.0,  # rad/s
)
```

`result["response"]` is complex and has shape
`(samples, frequencies, free_dofs)`.

### Time response

Times are in seconds. `force` may be a constant vector or a callable that
accepts time and returns a full or reduced force vector:

```python
result = solver.solve(
    "time_response",
    fixed_dofs=fixed,
    times=[0.0, 0.001, 0.002],
    force=force,
    speed=100.0,
)
```

An optional `initial_state` contains the free-DOF displacements followed by
the free-DOF velocities. The returned `displacement` and `velocity` arrays have
shape `(samples, times, free_dofs)`.

## Examples

Run the included examples from the repository root:

```bash
python examples/01_simple_shaft_modal.py
python examples/02_simple_shaft_modal_with_bearings.py
```

The theoretical conventions and implemented shaft formulation are documented
in [theory/00_notation.md](theory/00_notation.md) and
[theory/01_shaft_element.md](theory/01_shaft_element.md).

## Testing

Run the complete suite with:

```bash
python -m pytest
```

See [tests/README.md](tests/README.md) for the test organization, benchmark
assumptions, and guidance for adding verification cases.

## Current limitations

- Global matrices are dense.
- Sampling currently uses plain Monte Carlo rather than variance-reduction
  methods such as Latin hypercube sampling.
- Correlation is supported within a multivariate normal distribution, but not
  between separate distribution records.
- Analyses return numerical dictionaries; plotting and result-object APIs are
  not currently provided.

## Contributing

Keep dependencies minimal and accompany changes to formulations, records,
assembly, or solvers with focused tests. Numerical tolerances should be based
on physical or analytical expectations rather than adjusted solely to suppress
a failing test.

Spinniped is released under the MIT License.
