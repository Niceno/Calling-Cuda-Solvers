# T-Flows User Guide

## Introduction

T-Flows is a computational fluid dynamics (CFD) program for simulation of turbulent single and multiphase flows.  Numerical method is based on collocated finite volume method on unstructured arbitrary grids and turbulence models include a range of Reynods-Averaged Navier-Stokes (RANS) models, large eddy simulations (LES), as well as hybrid RANS-LES approach.

Multiphse models include an algebraic volume of fluid (VOF) method and Lagrangian particle tracking model.  Three-phase flows sitations (two fluid phases with VOF and one solid phase as particles) are also supported.

In T-Flows, Navier-Stokes equations are discretized in their incompressible form, meaning that pressure and temperatures are not link through an equation of state.  Flows with variable physical propertries, including flows with density variations, are supported.

## Minimum software requirements

- make utility
- Fortran 2008 compiler
- C compiler

T-Flows is mostly written in Fortran 2008 (only one function is written in C) and the compilation is controlled by makefiles.  So, the the requirements listed above are a bare minimum for you to start using the code.

>[!WARNING]
>This is a warning

## Highly desirable software requirements

- GMSH
- any other free or commercial mesh generator exporting ANSYS' .msh format
- visualization software which can read .vtu file format such as Paraview or VisIt, or any commercial tool which can read .vtu file format
- OpenMPI installation (mpif90 for compilation and mpirun for parallel processing)

T-Flows is, in essence, the flow solver without any graphical user interface (GUI).  Although it comes with its own mesh generator, it is very rudimentary and an external software, either free or commercial, would be highly desirable for meshing of complex computational domains.  We regularly use GMSH and would highly recommend it for its inherent scipting ability, but if you have access to any commercial mesh generator which can export meshes in ANSYS' .msh (and .cas, this should be checked) format, that would just fine.  Having no GUI, T-Flows relies on external tools for visualisation of results.  The results are saved in .vtu, Paraview's unstructured data format, and any visualisation software which can read that format is highly desirable for post-processing of results.

From its beginnings, T-Flows was developed for parallel execution with Message Passing Interface (MPI).  If you inted to run it on parallel computational platforms, you will also need an installation of OpenMPI on your system.  

## Optional software packages

- xmgrace
- PETSc

Visualization tools such as ParaView and VisIt are powerful, self-contained and sufficient for all sorts of post-processings, occasionally you might want to extract profiles from your solution fields and compare them agains experiments or direct numerical simulation (DNS) results, so a two-dimensional plotting software might come handy.  We find xmgrace light particularly suitable for that purpose and many test cases which come with T-Flows, already have benchmark cases compared in xmgrace's format.

Although T-Flows comes with its own suite of linear solvers based on Krylov sub-space family of methods (Incomplete Cholesky and Jacobi preonditioned CG, BiCG and  CGS), to have a better scaling with problem size, you may want to have more choice or even use algebraic multigrid preconditioners available through PETSc.  If PETSc is available on your system, T-Flows' makefiles will link with them and you will have all PETSc solvers at your disposal.









## Welcome to GitHub Pages

You can use the [editor on GitHub](https://github.com/Niceno/Calling-Cuda-Solvers/edit/gh-pages/index.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [Basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/Niceno/Calling-Cuda-Solvers/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.
