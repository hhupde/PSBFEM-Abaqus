# PSBFEM-Abaqus
Title: PSBFEM-Abaqus: development of User Element Subroutine (UEL) for polygonal scaled boundary finite element method in Abaqus

The polygonal scaled boundary finite element method (PSBFEM) is a novel method integrating the standard scaled
boundary finite element method (SBFEM) and the polygonal mesh technique. This work discusses developing a PSBFEM
framework within the commercial finite element software Abaqus. The PSBFEM is implemented by the User Element
Subroutine (UEL) feature of the software. The details on the main procedures to interact with Abaqus, defining the UEL
element, and solving the stiffness matrix by the eigenvalue decomposition are present. Moreover, we also develop the
preprocessing module and the postprocessing module using the Python script to generate meshes automatically and visualize
results. Several benchmark problems from two-dimensional linear elastostatics are solved to validate the proposed
implementation. The results show that PSBFEM-UEL has significantly better than FEM convergence and accuracy rate with
mesh refinement. The implementation of PSBFEM-UEL can conveniently use arbitrary polygon elements by the polygon/
quadtree discretizations in the Abaqus.

## Citation
If you use PSBFEM-Abaqus for academic research, you are encouraged to cite the following paper:

```
@article{ye2021psbfem,
  title={PSBFEM-Abaqus: Development of User Element Subroutine (UEL) for Polygonal Scaled Boundary Finite Element Method in Abaqus},
  author={Ye, Nan and Su, Chao and Yang, Yang},
  journal={Mathematical Problems in Engineering},
  volume={2021},
  year={2021},
  publisher={Hindawi}
}
```
