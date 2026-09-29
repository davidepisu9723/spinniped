# Chapter 12 -- Development notes

## Table of contents

- [Implementation rule](#implementation-rule)
- [Contribution rule](#contribution-rule)

## Implementation rule

Every new local matrix must declare its active coordinate frame and DOF order.
Two-grid shaft matrices must ultimately be expressed in the 12-DOF ordering in
Section 4, and one-grid matrices in the six-DOF ordering in Section 3, before
the builder assembles them.

## Contribution rule

Keep dependencies minimal and accompany changes to formulations, records,
assembly, or solvers with focused tests. Numerical tolerances must follow
physical or analytical expectations rather than being adjusted only to suppress
a failure.

The detailed verification organization is maintained in
[Testing and verification](../../tests/README.md).
