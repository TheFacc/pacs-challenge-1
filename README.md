# PACS Challenge 1 - Gradient Descent

This repo implements a simple C++ optimization framework that aims to solve the optimization problem:
$$\mathbf{x} = \argmin_{\mathbf{y} \in \mathbb{R}^n} f(\mathbf{y})$$

by implementing a simple 2D gradient descent strategy with initial step $\alpha_0$ that gets iteratively updated.

The program allows users to input their objective function and gradients through a text file, and set optional optimization parameters through commandline. The program additionally allows for some options for step size strategy and gradient computation.

## Prerequisites

To run this program, you will need:
- A C++ compiler supporting C++17 or later (e.g., GCC, Clang).
- The muParser library installed on your system for parsing mathematical expressions from text. If muParser is not installed, you can typically install it through your distribution's package manager. For example, on Ubuntu, you might use `sudo apt-get install libmuparser-dev`

Note: For Nix users, a `shell.nix` file is included. Run `nix-shell` in the project directory to enter an environment where all dependencies are available.

## Usage

1. **Edit Makefile**: Note that the Makefile included in this project is set up for NixOS and specifies a specific path to the muParser library:

    ```make
    LDFLAGS = -L/nix/store/...-muparser-2.3.2/lib -lmuparser
    ```

   **You may need to edit this path to match your system's configuration.** A common default path for other distributions might be `/usr/local/lib` or `/usr/lib`:

    ```make
    LDFLAGS = -L/usr/local/lib -lmuparser
    ```

2. **Compile the Program**: Run the `make` command in the directory containing the project files.

3. **Define objective function file** `function.txt`:
    - For gradient strategy `FiniteDifference`, you can include just the objective function:
        ```txt
        x1 * x2 + 4 * x1^4 + x2^2 + 3 * x1
        ```
   - For `Exact` strategy, include both the objective function and its gradients:
        ```txt
        x1 * x2 + 4 * x1^4 + x2^2 + 3 * x1
        x2 + 16 * x1^3 + 3
        x1 + 2 * x2
        ```
4. **Running the Program** (with default parameters):
    ```bash
    ./main function.txt
    ```

## Options

The program supports various optional parameters to control the optimization process. Some can be given via commandline, some are hardcoded as `constexpr` as suggested by the challenge.

#### **Commandline options**:
-   `--x0=<x0>,<y0>`: Initial point for the optimization process. Replace `<x0>` and `<y0>` with the starting values for variables `x1` and `x2`, respectively. (Default: `0,0`)

-   `--alpha0=<alpha0>`: Initial step size for the optimization. Replace `<alpha0>` with your desired initial step size.  (Default: `1.0`)

-   `--beta=<beta>`: Step reduction factor used in all step strategies in various ways.  (Default: `0.5`)

-   `--sigma=<sigma>`: Decrease constant specific for the Armijo rule.  (Default: `1e-4`)

-   `--resTolerance=<resTolerance>`: tolerance $\epsilon_r$ for control on the residual (Default: `1e-6`): $$|f(\mathbf{x}_{k+1})-f(\mathbf{x}_k)|<\epsilon_r$$

-   `--stepTolerance=<stepTolerance>`: tolerance $\epsilon_s$ for control on the step length (Default: `1e-6`): $$|| \mathbf{x}_{k+1} - \mathbf{x}_k||<\epsilon_s $$

-   `--maxIterations=<maxIterations>`: Maximum number of iterations for the optimization. (Default: `1000`)

#### **Hardcoded options**:

Two options are hardcoded  as `constexpr` and must be edited in the `main.cpp` in lines 19-20 before compiling. This is for efficiency reasons.

- `GradientStrategy` can be either:

    - `GradientStrategy::Exact` - use the exact gradient provided by the user in the input file
    - `GradientStrategy::FiniteDifference` - use finite differences to approximate gradient

- `StepSizeStrategy` can be either:

    - `StepSizeStrategy::ExponentialDecay` - implement the Exponential decay rule: $\alpha_k = \alpha_0 e^{-\beta k}$
    - `StepSizeStrategy::InverseDecay` - implement the Inverse decay rule: $\alpha_k = \alpha_0 / (1+\beta k)$
    - `StepSizeStrategy::Armijo` - implement the Armijo rule to solve approximate line search: set $\alpha_k = \alpha_0 \beta^k$ while the following decrease condition is not satisfied: $$f(\mathbf{x}_{k+1})-f(\mathbf{x}_k - \alpha_k \nabla f(\mathbf{x}_k)) \ge \sigma \alpha_k ||\nabla f(\mathbf{x}_K)||$$




## Examples

Run optimization with a custom starting point and maximum 100 iterations:

```bash
./main function.txt --x0=1.,2.0 --maxIterations=100
```
Run optimization with custom settings for step size, reduction factor, and tolerances:

```bash
./main function.txt --alpha0=0.5 --beta=0.8 --sigma=0.001 --resTolerance=1e-3 --stepTolerance=1e-2
```

Example output:
```
$ make clean && make && ./main function.txt
  ** PACS Challenge 1 **  

Reading file: function.txt
 - Objective function: x1 * x2 + 4 * x1^4 + x2^2 + 3 * x1 
 - Gradient component 1: x2 + 16 * x1^3 + 3 
 - Gradient component 2: x1 + 2 * x2 

Optimization parameters:
 - Initial point: x1 = 0, x2 = 0
 - Initial step size (alpha_0): 1
 - Step reduction factor (beta): 0.5
 - Decrease constant (sigma): 0.0001
 - Step size tolerance: 1e-06
 - Function residual tolerance: 1e-06
 - Maximum number of iterations: 1000

Starting optimization...
 - Gradient computation: exact (given from file)
 - Step size strategy: Armijo rule.
 ...completed after 44 iterations.
Argmin: x1 = -0.590438, x2 = 0.293602
```