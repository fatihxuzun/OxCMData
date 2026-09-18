# OxCMData: Automated Data Processing and Mesh Generation for Contour Method Residual Stress Evaluation

[![MATLAB](https://img.shields.io/badge/Platform-MATLAB%20R2020b+-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![FEniCS](https://img.shields.io/badge/Solver-FEniCS%20%7C%20OxCM-orange.svg)](https://fenicsproject.org/)
[![Status](https://img.shields.io/badge/Status-Active%20Development-success.svg)](#active-development--roadmap)
[![License: Academic & Non-Commercial](https://img.shields.io/badge/License-Academic%20%26%20Non--Commercial-blue.svg)](LICENSE)

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

```text
[ Plane-A & Plane-B Raw TXT ]
             │
             ▼
┌───────────────────────────────┐
│ STEP 1: OxCM_profilometryData │ ──► Clean boundary, Align (ICP), Interpolate,
└───────────────────────────────┘     Average & B-Spline Smooth  --> [ myData.txt ]
             │
             ▼
┌─────────────────────────────────┐
│ STEP 2: OxCM_rotateProcessedData│ ──► (Optional) Rotate in XY to align principal
└─────────────────────────────────┘     stress axes & re-center  --> [ myData.txt ]
             │
             ▼
┌─────────────────────────────────┐
│ STEP 3: OxCM_extrudedPerimeterMesh│ ─► Alpha-shape boundary, 2D remesh, 3D layer
└─────────────────────────────────┘     extrusion, Jacobian check --> [ myMesh.xml ]
             │
             ▼
┌───────────────────────────────┐
│    STEP 4: OxCM_plotXMLMesh   │ ──► Inspect boundary faces & verify 3D mesh
└───────────────────────────────┘
             │
             ▼
   [ OxCM / FEniCS Solver ]
```

---

## 💻 Prerequisites & Toolboxes

The pipeline requires **MATLAB (R2020b or later recommended)** with the following standard toolboxes:

- **Parallel Computing Toolbox** *(for parallel ICP registration)*
- **Curve Fitting Toolbox** *(for bivariate B-spline approximation `spap2` / `fnval`)*
- **Statistics and Machine Learning Toolbox** *(for `rmoutliers`, `knnsearch`, and KDTree search)*
- **Partial Differential Equation (PDE) Toolbox** *(optional: only needed if specifying custom `hmax`/`hmin` 2D remeshing)*

---

## 🚀 Step-by-Step Usage Guide

### **Step 1: Data Cleaning, Alignment & Surface Smoothing**

`OxCM_profilometryData` ingests the two opposing cut surfaces (Plane-A and Plane-B), removes EDM artefacts, performs 3D hierarchical ICP alignment, interpolates both datasets onto a common structured Cartesian grid, averages opposing deformations, applies bivariate B-spline smoothing, and saves the resulting point cloud.

```matlab
OxCM_profilometryData(nCpu, boundarClean, outlierClean, gx, gy, smth, flipSecondData, outputFilename)
```

#### Parameters:

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `nCpu` | Integer | System Cores | Number of parallel worker threads for parallel pool. |
| `boundarClean` | Integer | `1` | Number of perimeter boundary peel passes to purge EDM edge roll-off. |
| `outlierClean` | Integer | `1` | Number of MAD statistical outlier removal passes along the $Z$-axis. |
| `gx`, `gy` | Float | `0.1` | Common Cartesian interpolation grid spacing along $X$ and $Y$ (in mm). |
| `smth` | Float | `0.95` | B-spline smoothing factor $S \in [0, 1]$ ($0$: exact fit, $1$: maximum smoothing). |
| `flipSecondData`| String | `''none''` | Mirror plane-B if required: `''none''`, `''x''`, `''y''`, or `''xy''`. |
| `outputFilename`| String | `''myData.txt''` | Name of output text file containing smoothed $[x, y, z]$ coordinates. |

*When executed, interactive file-selection dialogs will prompt you to select the Plane-A and Plane-B raw `.txt` profilometry files.*

---

### **Step 2: In-Plane Coordinate Rotation (Optional)**

`OxCM_rotateProcessedData` rotates the 2.5D point cloud about the $Z$-axis and re-centers the coordinate system at the origin $(0, 0)$. This ensures the specimen's principal stress axes align precisely with the Cartesian axes of the finite element solver.

```matlab
OxCM_rotateProcessedData(rotationAngle, filename)
```

#### Parameters:

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `rotationAngle` | Float | `0` | In-plane rotation angle about the $Z$-axis in degrees. |
| `filename` | String | `''myData.txt''` | Path to the point cloud text file to rotate and update. |

---

### **Step 3: Arbitrary-Perimeter 3D Volumetric Extrusion**

`OxCM_extrudedPerimeterMesh` recovers the 2D bounding perimeter directly from the point cloud via $\alpha$-shapes, optionally remeshes the 2D planar footprint to target element sizes, extrudes the planar mesh downwards into $N$ structured parallel $Z$-layers over depth $D$, enforces positive Jacobian determinants, and writes a solver-ready FEniCS XML mesh.

```matlab
OxCM_extrudedPerimeterMesh(filename, extrusionDepth, hmax, hmin, scaleFactor, numLayersZ, outputFilename, elementType)
```

#### Parameters:

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `filename` | String | `''myData.txt''` | Input point cloud file generated by Step 1 / Step 2. |
| `extrusionDepth` | Float | `10` | Total extrusion depth along $-Z$ (in mm). |
| `hmax` | Float | `[]` | Target maximum 2D element edge length in XY (pass `[]` to use raw point spacing). |
| `hmin` | Float | `[]` | Target minimum 2D element edge length in XY (pass `[]` to use raw point spacing). |
| `scaleFactor` | Float | `1.0` | Inward scaling factor relative to centroid (typically `1.0`). |
| `numLayersZ` | Integer | `10` | Number of parallel structured $Z$-layers along the extrusion depth. |
| `outputFilename`| String | `''myMesh.xml''` | Name of output FEniCS-compatible XML mesh file. |
| `elementType` | String | `''tetrahedron''` | Mesh topology: `''tetrahedron''` (default simplex for FEniCS) or `''hexahedron''`. |

---

### **Step 4: Mesh Quality Inspection & Visualization**

`OxCM_plotXMLMesh` parses the generated FEniCS XML mesh, detects the cell type (`tetrahedron`, `hexahedron`, or `triangle`), extracts the external boundary faces, and renders an interactive 3D view of the domain for rapid pre-solve validation.

```matlab
OxCM_plotXMLMesh(xmlFilename)
```

#### Parameters:

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `xmlFilename` | String | `''myMesh.xml''` | Path to the FEniCS XML mesh file to inspect and plot. |

---

## ⚡ Master Execution Script

Create a master run script (e.g., `run_pipeline.m`) in your root directory to execute all 4 steps sequentially:

```matlab
%% Master Execution Script for OxCMData Pipeline
clear; close all; clc;

%% Configuration Parameters
nCpu           = 4;              % CPU workers
boundarClean   = 2;              % Perimeter peel passes (EDM roll-off removal)
outlierClean   = 1;              % Z-height outlier passes
gx             = 0.1;            % Grid spacing X (mm)
gy             = 0.1;            % Grid spacing Y (mm)
smth           = 0.95;           % B-spline smoothing parameter S in [0, 1]
flipSecondData = 'none';         % Plane-B coordinate flip ('none', 'x', 'y', 'xy')
pointCloudFile = 'myData.txt';   % Intermediate point cloud file

rotAngle       = 0.0;            % In-plane rotation angle (deg)

extrusionDepth = 50.0;           % Extrusion depth D (mm)
hmax           = 1.0;            % Max element size (mm) (or [] for raw cloud)
hmin           = 0.2;            % Min element size (mm)
numLayersZ     = 12;              % Structured Z-layers (converges around 10-15 layers)
meshFile       = 'myMesh.xml';   % Final FEniCS mesh file
elemType       = 'tetrahedron';  % 'tetrahedron' or 'hexahedron'

%% STEP 1: Process Profilometry Data
disp('--- [STEP 1] Running Profilometry Processing & Alignment ---');
OxCM_profilometryData(nCpu, boundarClean, outlierClean, gx, gy, smth, flipSecondData, pointCloudFile);

%% STEP 2: In-Plane Rotation (Optional)
disp('--- [STEP 2] Rotating and Re-Centering Point Cloud ---');
OxCM_rotateProcessedData(rotAngle, pointCloudFile);

%% STEP 3: Extrude Boundary & Generate 3D Mesh
disp('--- [STEP 3] Generating 3D Volumetric Extrusion Mesh ---');
OxCM_extrudedPerimeterMesh(pointCloudFile, extrusionDepth, hmax, hmin, 1.0, numLayersZ, meshFile, elemType);

%% STEP 4: Visualize & Verify Mesh
disp('--- [STEP 4] Plotting Final FEniCS XML Mesh ---');
OxCM_plotXMLMesh(meshFile);

disp('Pipeline execution complete! Mesh is ready for OxCM / FEniCS.');
```

---

## 📂 Public Datasets Evaluated

The fidelity of this pipeline was validated using three open-access datasets:

1. **pyCM Benchmark Specimen:** Standard baseline test case with minimal cutting artefacts.  
   *Available on GitHub:* [majroy/pyCM](https://github.com/majroy/pyCM)
2. **ENPOWER Austenitic Edge-Welded Beam:** Complex beam geometry with severe EDM wire-bowing and edge roll-off.  
   *Available on Zenodo:* [doi.org/10.5281/zenodo.3373916](https://doi.org/10.5281/zenodo.3373916)
3. **Clad Pressure Vessel Steel:** Non-rectangular dissimilar weldment with steep interfacial stress gradients.  
   *Available on Zenodo:* [doi.org/10.5281/zenodo.4003865](https://doi.org/10.5281/zenodo.4003865)

---

## 🛠️ Active Development & Roadmap

> **Notice:** This repository represents the initial prototype of the **OxCMData** pipeline and is **actively maintained and updated**.

Future releases and enhancements are planned and will be pushed directly to this repository:
- [ ] **Multi-Language Implementations:** Standalone versions in **Python** (NumPy/SciPy) and **C++** for seamless cross-platform deployment without MATLAB dependencies.
- [ ] **FEniCS Support:** Direct export routines to modern HDF5/XDMF formats for native parallel execution in **DOLFIN / FEniCS**.
- [ ] **WEDM Parameter Calibration:** Automated boundary erosion and filtering rules calibrated for different wire electrical discharge machining (WEDM) cut settings.
- [ ] **Direct Gmsh/meshio Interoperability:** Extended export bindings to standard finite element formats (`.inp`, `.msh`, `.vtk`).

---

## 📖 Citation

If you use **OxCMData** or find this pipeline helpful in your research, please cite:

```bibtex
@article{Uzun2024OxCMData,
  title   = {OxCMData: Automated Data Processing and Mesh Generation for Contour Method Residual Stress Evaluation},
  author  = {Uzun, Fatih},
  journal = {International Journal of Pressure Vessels and Piping},
  year    = {2026}
}

@article{Uzun2024OxCM,
  title   = {The OxCM contour method solver for residual stress evaluation},
  author  = {Uzun, Fatih and Korsunsky, Alexander M.},
  journal = {Engineering with Computers},
  volume  = {40},
  pages   = {3059--3072},
  year    = {2024},
  doi     = {10.1007/s00366-024-01959-3}
}
```

---

## 📄 License & Commercial Inquiries

This project is licensed under the **OxCMData Academic & Non-Commercial Research License**:

- **Academic & Research Use:** Completely **free of charge** for universities, researchers, and students, provided that appropriate citations are included in any resulting publications.
- **Commercial & Industrial Use:** Strictly **prohibited without a separate commercial license**. If your company or engineering consultancy wishes to use, deploy, or integrate `OxCMData` for proprietary or commercial projects, please contact:
  - **Fatih Uzun** — 📧 `fatihuzun@me.com` | `fatih.uzun@eng.ox.ac.uk`

---

## ✉️ Contact & Support

For bug reports, feature requests, or questions regarding contour method data processing, please open an issue on GitHub or contact:

- **Fatih Uzun**  
  📧 `fatihuzun@me.com`
