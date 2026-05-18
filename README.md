# Hackathon Squad - Optimization Solver

![C++17](https://img.shields.io/badge/C%2B%2B-17-blue)
![Platform](https://img.shields.io/badge/Platform-Windows-0078D6)
![Compiler](https://img.shields.io/badge/Compiler-MinGW%20g%2B%2B-orange)
![Test Runner](https://img.shields.io/badge/Test%20Runner-Node.js-339933)
![Graph Type](https://img.shields.io/badge/Graph-Sparse-lightgrey)

## Problem Statement

The **Hackathon Squad** problem models the task of assembling a dream team from a very large pool of freshman coders for a university hackathon.

Each coder has a **Skill Rating**, and some pairs of coders are incompatible because of known rivalries or conflicts. A valid team must therefore satisfy two conditions:

- It should include only coders who can work together.
- It should maximize the total sum of selected skill ratings.

Technically, this is the **Maximum Weight Independent Set (MWIS)** problem on a sparse, weighted, undirected graph:

- Each coder is a vertex.
- Each conflict pair is an edge.
- Each skill rating is the weight of a vertex.
- The objective is to find an independent set with maximum total weight.

## Constraints (Crucial)

1 <= N <= 200,000
0 <= M <= (N * (N - 1)) / 2
1 <= Sᵢ <= 1,000,000,000
1 <= u, v <= N
All conflict pairs (u, v) are distinct.
- The graph is therefore **sparse**, which is a key structural constraint exploited by the solver
- Execution time limit: **5 minutes**

## Deliverables

The final submission artifact is a highly optimized, single-file C++ execution engine:

- `merged_solution.cpp`

This merged file is designed for direct compilation and final evaluation. It dynamically routes graph components to exact or heuristic solvers based on component topology, enabling strong performance across forests, bipartite regions, dense pockets, and large irregular sparse subgraphs.

For development and maintenance, the repository also contains modular source files under:

- `src/`

## Algorithm Architecture (Briefly)

The solver first decomposes the full graph into isolated connected components. Each component is then dispatched to the most appropriate solving strategy:

- **Exact Tree DP** for forest components
- **Dinic's Algorithm (Bipartite Matching / Min-Cut formulation)** for bipartite components
- **Branch and Bound** for small dense components
- **Multi-start Simulated Annealing** and **Heavy Insert** heuristics for large, complex components

This hybrid architecture allows the program to remain exact where the topology is favorable and aggressively optimized where the search space is too large for exhaustive methods.

## Prerequisites

Use a Windows machine with the following installed:

- **VS Code**
- **MinGW** with `g++` added to your **System PATH**
- **Node.js** for running the automated test script

Before compiling, open a new **PowerShell** window in the repository root and confirm the tools are visible:

```powershell
g++ --version
node --version
```

If either command is not recognized, install the missing tool and restart PowerShell before continuing.

## Complete Instructions (For a beginner)

Assumption: your PowerShell terminal is already opened in the repository root, and you can see files such as `merged_solution.cpp`, `run_tests.js`, and the `test_suite` folder.

### Step 1: Compiling the Code

Run the following command to compile the final submission file into a Windows executable:

```powershell
g++ -std=c++17 -O2 -o .\hackathon_squad_solver.exe .\merged_solution.cpp
```

What this does:

- `-std=c++17` enables the C++17 language standard
- `-O2` turns on compiler optimizations for faster execution
- `-o .\hackathon_squad_solver.exe` names the output executable

If compilation succeeds, you will see `hackathon_squad_solver.exe` in the repository root.

### Step 2: Running a Single Test

To run the solver on a single input file and print the answer to the terminal:

```powershell
Get-Content .\input.txt | .\hackathon_squad_solver.exe
```

Example using one of the repository test files:

```powershell
Get-Content .\test_suite\input_01_tiny_random.txt | .\hackathon_squad_solver.exe
```

If you want to save the program output to a file while also viewing it, you can use:

```powershell
Get-Content .\test_suite\input_01_tiny_random.txt | .\hackathon_squad_solver.exe | Tee-Object .\actual_output.txt
```

### Step 3: Running the Automated Test Suite

To validate the solver against all bundled test cases and expected outputs, run:

```powershell
node .\run_tests.js .\test_suite
```

This command executes the automated test runner, feeds every input in `test_suite\` to the solver, and compares the produced output with the corresponding expected output files.
