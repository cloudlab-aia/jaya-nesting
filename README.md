# Jaya Nesting

This project is a Jaya Nesting algorithm implementation in C++.

## Prerequisites

- Visual Studio with C++ support
- CMake (optional, if using CMake for build)

## Installation

1. Clone the repository:
    ```sh
    git clone https://github.com/cloudlab-aia/jaya-nesting.git
    cd jaya-nesting
    ```

2. Open the project in Visual Studio:
    - Open Visual Studio.
    - Select `File > Open > Project/Solution`.
    - Navigate to the `jaya-nesting` directory and select the `.sln` file.

## Running the Algorithm

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

*To run the nesting you will need some figures in dxf format. The dxf should be structured as follows:*

- *Each figure must be in a separate layer*
- *Each figure must be composed only by simple lines*
- *Any other entity type will be ignored*

## Datasets

In the datasets folder there are some dxf files with figures from multiple datasets that can be used with the algorithm.

- albano: [Optimal Allocation of Two-Dimensional Irregular Shapes Using Heuristic Search Methods](https://doi.org/10.1109/TSMC.1980.4308483)
- jakobs: [On genetic algorithms for the packing of polygons](https://doi.org/10.1016/0377-2217(94)00166-9)
- swim: [TOPOS – A new constructive algorithm for nesting problems](http://dx.doi.org/10.1007/s002910050105)
- lee: [A heuristic for nesting problems of irregular shapes](https://doi.org/10.1016/j.cad.2008.02.008)

## Acknowledgements

- This implementation is based on the paper: [An approach to apply the Jaya optimization algorithm to the nesting of irregular patterns](https://doi.org/10.1093/jcde/qwae093)
- Jaya original paper: [A simple and new optimization algorithm for solving constrained and unconstrained optimization problems](https://doi.org/10.5267/j.ijiec.2015.8.004)
