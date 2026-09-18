# OxCMData: Automated Data Processing and Mesh Generation for Contour Method Residual Stress Evaluation

[![MATLAB](https://img.shields.io/badge/Platform-MATLAB%20R2020b+-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![FEniCS](https://img.shields.io/badge/Solver-FEniCS%20%7C%20OxCM-orange.svg)](https://fenicsproject.org/)
[![Status](https://img.shields.io/badge/Status-Active%20Development-success.svg)](#active-development--roadmap)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**OxCMData** is an open-source, automated MATLAB computational framework designed to bridge the gap between raw surface profilometry and implicit finite element solvers (such as [OxCM / FEniCS](https://doi.org/10.1007/s00366-024-01959-3)). 

Historically, contour method data analysis and mesh generation have relied heavily on proprietary CAD tracing and commercial finite element pre-processors (Abaqus, ANSYS). **OxCMData** completely eliminates the need for commercial intermediation by providing an end-to-end, automated pipeline: from unaligned, noisy point clouds to solver-ready, Jacobian-validated 3D volumetric meshes exported natively to FEniCS XML format.

---

## 📌 Key Features

- **Automated Artefact Rejection:** Multi-iteration perimeter boundary erosion to purge electrical discharge machining (EDM) wire entry/exit roll-off, combined with median absolute deviation (MAD) statistical outlier filtering.
- **Parallelized Hierarchical ICP:** Multi-core point-to-point Iterative Closest Point (ICP) registration featuring synthetic planar pre-alignment, Delaunay/k-d tree spatial searches, and quaternion-based trajectory extrapolation.
- **Tunable Bivariate B-Spline Smoothing:** B-spline tensor-product surface approximation scaled dynamically by a single normalized smoothing parameter $S \in [0, 1]$ acting as a spatial low-pass filter.
- **CAD-Free $\alpha$-Shape Boundary Detection:** Autonomous detection of arbitrary, non-rectangular cross-sectional footprints using scaled $\alpha$-shape formulations ($1.2 \times r_\alpha$).
- **Structured 3D Volumetric Extrusion:** Generates parallel structured $Z$-strata extruded to user-defined depths.
- **Simplex & Hexahedral Mesh Generation:** Conforming prism-to-tetrahedron subdivision for standard legacy FEniCS solvers, alongside tri-to-quad barycentric subdivision for hexahedral meshes.
- **Strict Jacobian Positivity Enforcement:** Automated audit and correction of inverted elements to guarantee valid element topologies for finite element analysis.

---

## 🔄 Workflow Pipeline Overview
