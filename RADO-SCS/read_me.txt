Introduction

RADO-SCS model is the first open source spinal cord model for simulating conventional/High Frequency Spinal Cord Stimulation, Dorsal Root Ganglion (DRG) Stimulation, and transpinal Direct Current Stimulation (tsDCS). The RADO-SCS model emerged over years. In our very first version of the SCS model, RADO-SCS 1.0, cylindrical shapes were used to represent the spinal tissues which later got upgraded with an extensive research and modeling works to RADO-SCS 2.0. This version included detailed spinal tissues anatomy but the spinal roots, rootlets, sympathetic chain, DRG, and vasculatures were not included in that model.

Freely available online, the RADO-SCS will be updated continuously with version control.

Citation

This open-source model is made available based on Creative Commons Attribution 4.0 International (CC BY 4.0) License. Please use the following citation when referring to our open-source model.

Khadka, N., Liu, X., Zander, H., Swami, J., Rogers, E., Lempka, S.F., Bikson, M., 2019. Realistic Anatomically Detailed Open-Source Spinal Cord Stimulation (RADO-SCS) Model. bioRxiv 857946. https://doi.org/10.1101/857946.

Significance of the RADO-SCS model

The RADO-SCS 3.0 can be used to simulate any SCS approach with both unprecedented resolution (precision) and transparency (reproducibility). Compared to prior models, RADO-SCS model meets or exceeds detail for every tissue compartment. The resulting electric fields in white and gray-matter, and axon model activation thresholds are broadly consistent with prior simulations.

Getting Started

Before you download the RADO-SCS model for simulation, we suggest to read “forward modeling workflow for current flow models” which will provide detailed help with a general FEM modeling Pipeline using CAD-derived/MRI-derived models. The current flow modeling pipelines have been published many times in our prior publications (https://doi.org/10.1101/704940). When ready, you can download the latest version of the RADO-SCS model ( V 3.0) and use it for your exploratory clinical or academic research.

Modeling Workflow

- Download all STL files as they represent realistic and anatomically detailed T9-T11 level spinal tissues, vertebraes, IV discs, soft-tissues, vasculatures, and clinical leads (pre-positioned). If you don’t need all the detailed tissues, you can remove them from your final model.

- If you want to position the SCS leads in different location, import all STLs using CAD software like SolidWorks 2016 or later (recommended), Abacus, Mimics, etc., and re-position the leads.

- Save your final model files and import them either in free or commercial meshing software. We recommend Simpleware or other powerful meshing software because the model has sub-millimeter resolution.

- After the model is meshed, import the model into finite element method (FEM) or finite element analysis (FEA) current flow simulating software such as COMSOL, Abaqus, Mimics, etc.

- If interested in fiber activation thresholds simulation, export the extracellular voltages calculated from the FEM model and couple it into the NEURON/ Matlab modeling platform.

- Source axon trajectories and coordinates are also provided/

Notes and Disclaimers

- Version 3.0 is the most updated version of the RADO-SCS model.

- RADO-SCS does not include a pathological human spinal cord, but there are plans to incorporate this in future versions.

- This open-source model only includes CAD derived model files with pre-positioned SCS leads. However, the users can change the lead position using commercial or open-source CAD-software as their need.

- RADO-SCS is not a modeling open-source software. It is a CAD-derived model of realistic human spinal anatomy which can be downloaded and simulated for educational, clinical, or exploratory purposes. See the modeling pipeline on how to use the open-source RADO-SCS for numerical simulation.

- There is a plan to develop an open-source current flow modeling platform using the most-updated version of the RADO-SCS model in future.