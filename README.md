<h1 aligh="center">An approach to apply the Jaya optimization algorithm to the nesting of irregular patterns</h1>
<p align="center">An approach to apply the Jaya optimization algorithm to the nesting of irregular patterns</p>

<p align="center">
  <a href="https://doi.org/10.1093/jcde/qwae093">
    <img src="https://img.shields.io/badge/Journal%20of%20Computational%20Design%20and%20Engineering-2024-orange" alt="Feature Requests">
  </a>
</p>

This repository contains a Jaya Nesting algorithm implementation in C++ as seen in *An approach to apply the Jaya optimization algorithm to the nesting of irregular patterns*.

## Contents
This repository contains the code of the implementation used in the research. Inside you can find a Visual Studio project, and a folder with the datasets used in the paper.

Folders:
- root: Code and modules.
- datasets: Folder containing the datasets in dxf format.

## Requirements
- Visual Studio
- CMake (optional, if using CMake for build)

## Installation and use
This project is prepared to be run inside the Visual Studio environment. We do not provide a standalone deployment of the prototype.

1. Clone the repository:
    ```sh
    git clone https://github.com/cloudlab-aia/jaya-nesting.git
    cd jaya-nesting
    ```

2. Open the project in Visual Studio:
    - Open Visual Studio.
    - Select `File > Open > Project/Solution`.
    - Navigate to the `jaya-nesting` directory and select the `.sln` file.

To run the projecto follow this instructions:

1. Ensure you have the project open in Visual Studio.

2. Build the project:
    - Select `Build > Build Solution` or press `Ctrl+Shift+B`.

3. Run the project:
    - Select `Debug > Start Without Debugging` or press `Ctrl+F5`.

4. Run the compiled executable:
    ```sh
    ./jaya-nesting <figures_to_nest>.dxf <area>.dxf <population> <max_iterations>
    ```
    Ensure you provide the correct arguments. If the arguments are incorrect, you will see the following usage message:
    ```sh
    Uso: ./jaya-nesting figuras.dxf area.dxf population iterations
    ```

*To run the nesting you will need some figures in dxf format. You can use one of the datasets provided as is. If you want to use a custom file, the dxf should be structured as follows:*

- *Each figure must be in a separate layer*
- *Each figure must be composed only by simple lines*
- *Any other entity type will be ignored*

## Data

In the datasets folder there are some dxf files with figures from multiple datasets that can be used with the algorithm.

- albano: [Optimal Allocation of Two-Dimensional Irregular Shapes Using Heuristic Search Methods](https://doi.org/10.1109/TSMC.1980.4308483)
- jakobs: [On genetic algorithms for the packing of polygons](https://doi.org/10.1016/0377-2217(94)00166-9)
- swim: [TOPOS – A new constructive algorithm for nesting problems](http://dx.doi.org/10.1007/s002910050105)
- lee: [A heuristic for nesting problems of irregular shapes](https://doi.org/10.1016/j.cad.2008.02.008)

## Acknowledgements
This research has been performed for the research project <a href="https://aia.ua.es/en/proyectos/federated-serverless-architectures-for-heterogeneous-high-performance-computing-in-smart-manufacturing.html" target="_blank">Federated Serverless Architectures for Heterogeneous High Performance Computing in Smart Manufacturing</a>, at the [UniCAD: CAD/CAM/CAE](https://www.researchgate.net/lab/UniCAD-Antonio-Jimeno-Morenilla) Research Group of the University of Alicante (Spain).

Grant <b>Serverless4HPC PID2023-152804OB-I00</b> funded by MICIU/AEI/10.13039/501100011033 and by ERDF/EU.


## Citation
```bibtex
@article{10.1093/jcde/qwae093,
    author = {Duta, Eduard-Andrei and Jimeno-Morenilla, Antonio and Sanchez-Romero, Jose-Luis and Macia-Lillo, Antonio and Mora-Mora, Higinio},
    title = {An approach to apply the Jaya optimization algorithm to the nesting of irregular patterns},
    journal = {Journal of Computational Design and Engineering},
    volume = {11},
    number = {6},
    pages = {112-121},
    year = {2024},
    month = {10},
    abstract = {The problem of nesting frequently arises in the industrial environment, and it has a strong ecological and economic impact in the manufacturing processes. It basically consists of placing a set of pieces (polygons) on a material sheet, making sure that the pieces do not overlap and that they do not exceed the boundaries of the sheet. With regard to irregular 2D polygons, the problem is NP-complete. Therefore, different heuristics have been developed so as to cope with the problem. In this paper, the application of the Jaya metaheuristic algorithm to the nesting problem is proposed. This algorithm has been already applied to several engineering problems and has generally demonstrated better results than most metaheuristic algorithms. In this paper, the Jaya algorithm has been adapted to the specific features of the nesting problem so as to optimize the placement of pieces on a sheet, with the objective of minimizing material waste and computational time. The results of our experimentation demonstrate the algorithm’s effectiveness in reducing the convex hull area across various datasets, showing potential in solving complex, irregular shape nesting problems. This research provides a new application of the Jaya algorithm and opens ways for future work in optimization techniques and parameter-free heuristic algorithms for nesting.},
    issn = {2288-5048},
    doi = {10.1093/jcde/qwae093},
    url = {https://doi.org/10.1093/jcde/qwae093},
    eprint = {https://academic.oup.com/jcde/article-pdf/11/6/112/61212299/qwae093.pdf},
}
```

## License Information
This project is licensed under the <a href="LICENSE.txt">GPL-3 license</a>.
