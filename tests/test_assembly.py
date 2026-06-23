"""Tests for global matrix assembly and boundary-condition handling."""

import numpy as np
import pytest

from spinniped import element, node, rotor, solver


def _assemble(elements, node_ids):
    # This private test helper creates the smallest complete Spinniped workflow:
    # mesh -> rotor -> solver -> global matrix assembly.
    meshgrid = node.Grid()

    # User node IDs may be consecutive or deliberately irregular, depending on
    # the assembly behaviour under test. Coordinates are spaced along z only
    # to produce a simple valid mesh.
    for i, nid in enumerate(node_ids):
        meshgrid.add_node(id=nid, z=float(i))

    # Store the mesh and supplied elements inside one rotor model.
    r = rotor.Rotor()
    r.add_meshgrid(meshgrid)
    for e in elements:
        r.add_element(e)

    # Assemble the local element matrices into K_global and M_global.
    s = solver.Solver()
    s.add_rotor(r)
    s.assemble_global_matrices()
    return s


def test_one_element_global_matrices_equal_local_matrices(shaft_element):
    # With exactly one two-node element, the global and local DOF orderings are
    # identical. Assembly should therefore reproduce the local matrices.
    s = _assemble([shaft_element], [1, 2])

    assert np.allclose(s.K_global, shaft_element.K)
    assert np.allclose(s.M_global, shaft_element.M)


def test_two_elements_accumulate_at_shared_node(shaft_properties):
    # Build a three-node chain: node 1 -- e1 -- node 2 -- e2 -- node 3.
    p = shaft_properties
    e1 = element.ShaftElement(1, p.E, p.nu, p.rho, p.d, 0.5, 1, 2)
    e2 = element.ShaftElement(2, p.E, p.nu, p.rho, p.d, 0.5, 2, 3)
    s = _assemble([e1, e2], [1, 2, 3])

    # Three nodes with six DOFs each produce an 18 x 18 global matrix.
    assert s.K_global.shape == (18, 18)

    # The middle node belongs to both elements. Its global 6 x 6 diagonal block
    # must equal the sum of e1's second-node block and e2's first-node block.
    assert np.allclose(
        s.K_global[6:12, 6:12], e1.K[6:12, 6:12] + e2.K[:6, :6]
    )

    # Nodes 1 and 3 share no element, so assembly must not create a direct
    # stiffness coupling between their DOFs.
    assert np.allclose(s.K_global[:6, 12:18], 0.0)


def test_nonconsecutive_node_ids_are_mapped_correctly(shaft_properties):
    # User IDs need not be zero-based or consecutive. The solver must map IDs
    # 10 and 70 to its internal node positions before assembling the element.
    p = shaft_properties
    e = element.ShaftElement(1, p.E, p.nu, p.rho, p.d, p.L, 10, 70)
    s = _assemble([e], [10, 70])

    # With one element, correct ID mapping again makes global matrices equal to
    # the local matrices.
    assert np.allclose(s.K_global, e.K)
    assert np.allclose(s.M_global, e.M)


def test_bearing_matrix_is_inserted_in_node_xy_block():
    # Distinct values make transposition, permutation, and wrong-index errors
    # immediately visible. A symmetric matrix would hide some such mistakes.
    bearing = element.BallBearingElement(
        eid=1,
        n1=30,
        kxx=11.0,
        kxy=12.0,
        kyx=13.0,
        kyy=14.0,
    )

    # Node ID 30 is the second internal node, whose six global DOFs begin at
    # index 6. The bearing acts on its x and y translations: indices [6, 7].
    s = _assemble([bearing], [10, 30, 70])

    # This is the exact local bearing stiffness expected in K_global[6:8, 6:8].
    K_b = np.array([[11.0, 12.0], [13.0, 14.0]])
    assert np.array_equal(s.K_global[6:8, 6:8], K_b)

    # No shaft elements are present, so these four bearing entries must be the
    # only nonzero terms in the entire global stiffness matrix.
    assert np.count_nonzero(s.K_global) == 4


def test_constraints_reduce_matrices_and_do_not_duplicate(shaft_element):
    # A two-node shaft has 12 global DOFs before constraints are applied.
    s = _assemble([shaft_element], [1, 2])

    # DOF 1 is intentionally supplied more than once and in both calls. The
    # solver should treat fixed DOFs as a set rather than removing extra rows.
    s.apply_boundary_conditions([0, 1, 1])
    s.apply_boundary_conditions([1, 2])

    # The canonical stored list must contain three sorted, unique indices.
    assert s.fixed_dofs == [0, 1, 2]

    # Removing three constrained DOFs from twelve leaves nine free DOFs.
    assert s.K_ff.shape == (9, 9)
    assert s.M_ff.shape == (9, 9)


def test_element_referencing_unknown_node_is_rejected(shaft_properties):
    # The mesh contains IDs 1 and 2, while the element incorrectly references
    # node 99 as its second endpoint.
    p = shaft_properties
    e = element.ShaftElement(1, p.E, p.nu, p.rho, p.d, p.L, 1, 99)

    # pytest passes this test only if assembly raises the requested exception.
    # Matching the message also checks that the error identifies the bad node.
    with pytest.raises(ValueError, match="Node ID 99"):
        _assemble([e], [1, 2])
