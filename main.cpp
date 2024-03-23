#include "muParser.h"
#include <vector>
#include <iostream>
#include <functional> // function
#include <numeric> // inner_product
#include <cmath> // sqrt, fabs, pow
#include <fstream> // ifstream, getline
#include <algorithm> // find

enum class StepSizeStrategy {
    Armijo, // Armijo rule: f(x_k) - f(x_k+1) >= sigma * alpha * ||grad f(x_k)||^2
    ExponentialDecay, // Exponential decay: alpha = alpha0 * e^(-beta*k)
    InverseDecay, // Inverse decay: alpha = alpha0 / (1 + beta*k)
};
enum class GradientStrategy {
    Exact, // Exact gradient must be provided by the user in the file (lines 2 and 3)
    FiniteDifference // approximation
};
constexpr StepSizeStrategy stepStrategy = StepSizeStrategy::Armijo;
constexpr GradientStrategy gradStrategy = GradientStrategy::Exact;

struct OptimizationParams {
    std::vector<double> x0; // Starting point for optimization
    double alpha0;          // Initial step size (alpha_0)
    double beta;            // Step reduction factor (default .5)
    double sigma;           // Decrease constant (sigma)
    double stepTolerance;   // eps_s = tolerance for x,y step size
    double resTolerance;    // eps_r = function residual tolerance
    unsigned maxIterations; // Maximum number of iterations allowed

    // Constructor giving all parameters directly
    OptimizationParams(const std::vector<double>& initPoint, 
                       double alpha0, 
                       double beta, 
                       double sigma, 
                       double stepTol, 
                       double relTol, 
                       unsigned int maxIter)
                : x0(initPoint),
                alpha0(alpha0),
                beta(beta),
                sigma(sigma),
                stepTolerance(stepTol),
                resTolerance(relTol),
                maxIterations(maxIter) {}

    // Constructor by optional argc,argv with defaults
    OptimizationParams(int argc, char* argv[]) 
                : x0({0., 0.}),
                alpha0(1.0),
                beta(0.5),
                sigma(1e-4),
                stepTolerance(1e-6),
                resTolerance(1e-6),
                maxIterations(1000) {
        // Parse the command line options
        std::cout << "Optimization parameters:" << std::endl;
        for (int i = 2; i < argc; ++i) {
            std::string option = argv[i];
            if (option.find("--x0=") == 0) {
                std::string x0Str = option.substr(5);
                size_t commaPos = x0Str.find(',');
                if (commaPos != std::string::npos) {
                    std::string x1Str = x0Str.substr(0, commaPos);
                    std::string x2Str = x0Str.substr(commaPos + 1);
                    try {
                        double x1 = std::stod(x1Str);
                        double x2 = std::stod(x2Str);
                        x0 = {x1, x2};
                    } catch (const std::exception& e) {
                        std::cerr << "ERROR: Invalid x0 values: " << x0Str << "   (Expected format: x1,x2      -> using default 0,0)" << std::endl;
                        // return 1;
                    }
                } else {
                    std::cerr << "ERROR: Invalid x0 format: " << x0Str << "   (Expected format: x1,x2      -> using default 0,0)" << std::endl;
                    // return 1;
                }
            } else if (option.find("--alpha0=") == 0) {
                std::string alpha0Str = option.substr(9);
                try {
                    alpha0 = std::stod(alpha0Str);
                } catch (const std::exception& e) {
                    std::cerr << "ERROR: Invalid alpha0 value: " << alpha0Str << "   (Expected a number like 1.0      -> using default)" << std::endl;
                    // return 1;
                }
            } else if (option.find("--beta=") == 0) {
                std::string betaStr = option.substr(7);
                try {
                    beta = std::stod(betaStr);
                } catch (const std::exception& e) {
                    std::cerr << "ERROR: Invalid beta value: " << betaStr << "   (Expected a number like 0.5      -> using default)" << std::endl;
                    // return 1;
                }
            } else if (option.find("--sigma=") == 0) {
                std::string sigmaStr = option.substr(8);
                try {
                    sigma = std::stod(sigmaStr);
                } catch (const std::exception& e) {
                    std::cerr << "ERROR: Invalid sigma value: " << sigmaStr << "   (Expected a number like 1e-4      -> using default)" << std::endl;
                    // return 1;
                }
            } else if (option.find("--resTolerance=") == 0) {
                std::string resTolStr = option.substr(15);
                try {
                    resTolerance = std::stod(resTolStr);
                } catch (const std::exception& e) {
                    std::cerr << "ERROR: Invalid resTolerance value: " << resTolStr << "   (Expected a number like 1e-6      -> using default)" << std::endl;
                    // return 1;
                }
            } else if (option.find("--stepTolerance=") == 0) {
                std::string stepTolStr = option.substr(16);
                try {
                    stepTolerance = std::stod(stepTolStr);
                } catch (const std::exception& e) {
                    std::cerr << "ERROR: Invalid stepTolerance value: " << stepTolStr << "   (Expected a number like 1e-6      -> using default)" << std::endl;
                    // return 1;
                }
            } else if (option.find("--maxIterations=") == 0) {
                std::string maxIterStr = option.substr(16);
                try {
                    maxIterations = std::stoi(maxIterStr);
                } catch (const std::exception& e) {
                    std::cerr << "ERROR: Invalid maxIterations value: " << maxIterStr << "   (Expected an integer like 1000      -> using default)" << std::endl;
                    // return 1;
                }
            } else {
                std::cerr << "ERROR: Invalid option: " << option << "   (Ignoring)" << std::endl;
                // return 1;
            }
        }
        std::cout << " - Initial point: x1 = " << x0[0] << ", x2 = " << x0[1] << std::endl;
        std::cout << " - Initial step size (alpha_0): " << alpha0 << std::endl;
        std::cout << " - Step reduction factor (beta): " << beta << std::endl;
        std::cout << " - Decrease constant (sigma): " << sigma << std::endl;
        std::cout << " - Step size tolerance: " << stepTolerance << std::endl;
        std::cout << " - Function residual tolerance: " << resTolerance << std::endl;
        std::cout << " - Maximum number of iterations: " << maxIterations << std::endl;
        std::cout << std::endl;
    }
};

bool isValidExpression(const std::string& expr) {
    // Check for invalid characters
    const std::string validChars = "0123456789+-*/^() x1x2";
    for (char ch : expr) {
        if (validChars.find(ch) == std::string::npos)
            return false;
    }
    return true;
}
bool checkParserVariables(const mu::Parser& p, const std::vector<std::string>& vars) {
    // check if the parser uses only the correct variables
    const mu::varmap_type& usedVars = p.GetUsedVar();
    for (const auto& var : usedVars)
        if (std::find(vars.begin(), vars.end(), var.first) == vars.end()) {
            std::cerr << "Error: Unexpected variable '" << var.first << "' found." << std::endl;
            return false;
        }
    return true;
}
bool readFunction(const std::string& filename, 
                                mu::Parser& funcParser, 
                                std::vector<mu::Parser>& gradParsers) {
    // Read file with all due checks
    std::ifstream f(filename);
    if(!f.is_open()){
        std::cerr << "Error opening file: " << filename << std::endl;
        return false;
    }
    if(gradParsers.size() != 2){
        std::cerr << "Error: gradient parsers variable not correctly initialized (it must have size 2)" << std::endl;
        return false;
    }

    std::cout << "Reading file: " << filename << std::endl;

    std::string line;
    unsigned lineNumber = 0;
    while (getline(f, line)) {
        if (!isValidExpression(line)) {
            std::cerr << "Error: Invalid expression detected: " << line << std::endl;
            return false;
        }
        if (lineNumber == 0) {
            // Set the objective function
            funcParser.SetExpr(line);
            if(funcParser.GetExpr().empty()){
                std::cerr << "Error reading function from file: " << filename << std::endl;
                return false;
            }
            std::cout << " - Objective function: " << funcParser.GetExpr() << std::endl;
            if constexpr (gradStrategy != GradientStrategy::Exact)
                return true; // We don't need to read the exact gradient components
        } else if(lineNumber < 3){
            // Set the gradient components
            mu::Parser p;
            p.SetExpr(line);
            if(line.empty() || p.GetExpr().empty()){
                std::cerr << " - Error reading gradient component " << lineNumber << std::endl;
                return false;
            }
            std::cout << " - Gradient component " << lineNumber << ": " << p.GetExpr() << std::endl;
            gradParsers[lineNumber - 1] = p;
        } else {
            std::cout << " -- (Ignoring extra lines) --" << std::endl;
            break;
        }
        ++lineNumber;
    }
    f.close();
    // check if the number of gradient components is correct
    if (lineNumber < 3) {
        std::cerr << "Error reading gradient components from file - expected 2 components, got " << lineNumber - 1 << std::endl;
        return false;
    }

    // check if variables in funcParser are the same as in gradParsers
    // get variables from funcParser
    // const mu::varmap_type& usedVars = funcParser.GetUsedVar();
    // std::vector<std::string> vars;
    // for (const auto& var : usedVars)
    //     vars.push_back(var.first);
    // std::cout << "Variables used in objective function: ";
    // for (const auto& var : vars) {
    //     std::cout << var << " ";
    // }
    // std::cout << std::endl;
    
    std::vector<std::string> vars = {"x1", "x2"}; // set manually
    if(!checkParserVariables(funcParser, vars) || !checkParserVariables(gradParsers[0], vars) || !checkParserVariables(gradParsers[1], vars))
        return false;
    std::cout << std::endl;
    return true;
}

// wrappers for muParser
double objectiveWrapper(const std::vector<double>& x, mu::Parser& parser) {
    parser.DefineVar("x1", &const_cast<double&>(x[0]));  // Use const_cast to bypass muParser's non-const interface
    parser.DefineVar("x2", &const_cast<double&>(x[1]));
    return parser.Eval();
}
std::vector<double> gradientWrapper(const std::vector<double>& x, std::vector<mu::Parser>& gradParsers) {
    std::vector<double> grad(gradParsers.size());
    for (size_t i = 0; i < grad.size(); ++i) {
        gradParsers[i].DefineVar("x1", &const_cast<double&>(x[0]));
        gradParsers[i].DefineVar("x2", &const_cast<double&>(x[1]));
        grad[i] = gradParsers[i].Eval();
    }
    return grad;
}

double norm(const std::vector<double>& v) {
    // Vector Euclidean norm
    return sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.0));
}
double norm(const std::vector<double>& x, const std::vector<double>& x_new) {
    // Euclidean norm of the difference
    // Eigen::Vector2d would come in handy here...:)
    if (x.size() != 2 || x_new.size() != 2)
        throw std::invalid_argument("Error computing norm: vectors must be 2D.");
    return std::sqrt(std::pow(x[0] - x_new[0], 2) + std::pow(x[1] - x_new[1], 2));
}
std::vector<double> finiteDifferenceGradient(
                            const std::function<double(const std::vector<double>&)>& f,
                            const std::vector<double>& x,
                            double h = 1e-5) {
    std::vector<double> grad(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        // Create a copy of x and modify only the i-th element for the partial derivative approximation
        std::vector<double> x_plus_h = x;
        x_plus_h[i] += h;
        double f_plus_h = f(x_plus_h);
        double f_x = f(x);
        grad[i] = (f_plus_h - f_x) / h;  // Partial derivative approximation
    }
    return grad;
}



// Armijo optimization
template<StepSizeStrategy stepStrategy = StepSizeStrategy::Armijo, 
         GradientStrategy gradStrategy = GradientStrategy::Exact
         >
std::vector<double> armijoOptimization(
    const std::function<double(const std::vector<double>&)>& f,
    const std::function<std::vector<double>(const std::vector<double>&)>& grad_f,
    const OptimizationParams& params
) {
    // As per challenge requests:
    // We start from the initial point and use gradient descent to update our position.
    // At each step, we calculate the new point and check whether the function value has decreased sufficiently according to the Armijo rule.
    // If the function value hasn't decreased enough, we reduce our step size and try again.
    // We continue this process until either the maximum number of iterations is reached, the gradient norm falls below the set tolerance,
    //  or the change in function value between iterations falls below the relative tolerance.

    std::vector<double> x = params.x0; // Initial point
    double alpha = params.alpha0;
    double sigma = params.sigma;
    unsigned int iter = 0;

    std::cout << "Starting optimization..." << std::endl;
    if constexpr (gradStrategy == GradientStrategy::FiniteDifference)
        std::cout << " - Gradient computation: approx using finite differences." << std::endl;
    else
        std::cout << " - Gradient computation: exact (given from file)" << std::endl;

    if constexpr (stepStrategy == StepSizeStrategy::Armijo)
        std::cout << " - Step size strategy: Armijo rule." << std::endl;
    else if constexpr (stepStrategy == StepSizeStrategy::ExponentialDecay)
        std::cout << " - Step size strategy: exponential decay." << std::endl;
    else if constexpr (stepStrategy == StepSizeStrategy::InverseDecay)
        std::cout << " - Step size strategy: inverse decay." << std::endl;
    else
        std::cout << " - Unknown step size strategy." << std::endl;

    while (iter < params.maxIterations) {
        std::vector<double> grad;
        if constexpr (gradStrategy == GradientStrategy::Exact) {
            grad = grad_f(x);  // Compute exact gradient
        } else if constexpr (gradStrategy == GradientStrategy::FiniteDifference){
            grad = finiteDifferenceGradient(f, x); // approx gradient
        } else {
            std::cerr << "Error: Unknown gradient strategy." << std::endl;
            break;
        }
        
        // Compute step and check
        std::vector<double> x_new(x.size());  // New candidate point
        for (size_t i = 0; i < x.size(); ++i)
            x_new[i] = x[i] - alpha * grad[i];  // gradient descent

        if(norm(x, x_new) < params.stepTolerance)
            break; // check step length
        
        // Evals
        double f_new = f(x_new);  // eval new point
        double f_old = f(x);      // eval old point

        if constexpr (stepStrategy == StepSizeStrategy::Armijo) {
            // Armijo condition
            double grad_norm = norm(grad);
            // std::cout << "Step norm: " << norm(x, x_new) << ", f_old: " << f_old << ", f_new: " << f_new << ", grad_norm: " << grad_norm << std::endl;
            if (f_old - f_new >= params.sigma * alpha * (grad_norm * grad_norm)) {
                x = x_new;
                if (std::fabs(f_old - f_new) < params.resTolerance)
                    break; // Check residual
            } else {
                alpha *= params.beta;
            }
        } else {
            if (f_new < f_old) {
                x = x_new;
                if (std::fabs(f_old - f_new) < params.resTolerance)
                    break;
            } else if constexpr (stepStrategy == StepSizeStrategy::ExponentialDecay){
                // exponential decay: alpha = alpha0 * e^(-beta*k)
                alpha = params.alpha0 * std::exp(-params.beta * iter);
            } else if constexpr (stepStrategy == StepSizeStrategy::InverseDecay){
                // inverse decay: alpha = alpha0 / (1 + beta*k)
                alpha = params.alpha0 / (1 + params.beta * iter);
            } else {
                std::cerr << "Error: Unknown step size strategy." << std::endl;
                break;
            }
        }
        ++iter;
    }
    std::cout << " ...completed after " << iter << " iterations." << std::endl;
    return x;
}
template<StepSizeStrategy stepStrategy = StepSizeStrategy::Armijo, 
         GradientStrategy gradStrategy = GradientStrategy::Exact
         >
std::vector<double> armijoOptimization(
    mu::Parser p,
    std::vector<mu::Parser>& gradParsers,
    const OptimizationParams& params
) {
    // Wrapper function to call the optimization with muParser objects directly
    return armijoOptimization<stepStrategy, gradStrategy>(
        [&](const std::vector<double>& x) { return objectiveWrapper(x, p); },
        [&](const std::vector<double>& x) { return gradientWrapper(x, gradParsers); },
        params
    );
}

// example function
// double objectiveFunction(const std::vector<double>& x) {
//     return x[0] * x[1] + 4 * pow(x[0], 4) + pow(x[1], 2) + 3 * x[0];
// }
// std::vector<double> gradientFunction(const std::vector<double>& x) {
//     std::vector<double> grad(2);
//     grad[0] = x[1] + 16 * pow(x[0], 3) + 3;  // Partial derivative with respect to x1
//     grad[1] = x[0] + 2 * x[1];              // Partial derivative with respect to x2
//     return grad;
// }


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: ./armijo <filename> [options]" << std::endl;
        std::cout << "  <filename> - File containing the function and its derivatives" << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "  --x0=<x1_value>,<x2_value> - Initial point for optimization (default: 0,0)" << std::endl;
        std::cout << "  --alpha0=<alpha_value> - Initial step size (default: 1.0)" << std::endl;
        std::cout << "  --beta=<beta_value> - Step reduction factor (default: 0.5)" << std::endl;
        std::cout << "  --sigma=<sigma_value> - Decrease constant (default: 1e-4)" << std::endl;
        std::cout << "  --resTolerance=<resTol_value> - Function residual tolerance (default: 1e-6)" << std::endl;
        std::cout << "  --stepTolerance=<stepTol_value> - Step size tolerance (default: 1e-6)" << std::endl;
        std::cout << "  --maxIterations=<maxIter_value> - Maximum number of iterations allowed (default: 1000)" << std::endl;
        return 1;
    }

    std::cout << "  ** PACS Challenge 1 **  " << std::endl << std::endl;

    // Set function and derivatives using muParser:
    std::string filename = argv[1];
    mu::Parser p;
    std::vector<mu::Parser> gradParsers(2);  // For two variables, x1 and x2
    // p.SetExpr("x1 * x2 + 4 * x1^4 + x2^2 + 3 * x1");
    // gradParsers[0].SetExpr("x2 + 16 * x1^3 + 3");  // Partial derivative w.r.t. x1
    // gradParsers[1].SetExpr("x1 + 2 * x2");         // Partial derivative w.r.t. x2
    if(!readFunction(filename, p, gradParsers))
        return 1;

   // Set initial conditions and parameters:
   // - directly
    // std::vector<double> x0 = {0., 0.};
    // double alpha0 = 1.0;
    // double beta = 0.5;
    // double sigma = 1e-4;
    // double resTolerance = 1e-6;
    // double stepTolerance = 1e-6;
    // unsigned int maxIterations = 1000;
    // OptimizationParams params(x0, alpha0, beta, sigma, stepTolerance, resTolerance, maxIterations);

   // - using the optional argc,argv with defaults embedded inside
    OptimizationParams params(argc, argv);


   // Run the optimization:

   // - using the hardcoded example functions defined above
    // std::vector<double> res = armijoOptimization(objectiveFunction, gradientFunction, params);

   // - using the muParser functions and wrappers (challenge asked for a function taking 'wrappers' as input)
    // std::vector<double> res = armijoOptimization<StepSizeStrategy::Armijo, GradientStrategy::Exact>(
    //     [&](const std::vector<double>& x) { return objectiveWrapper(x, p); },
    //     [&](const std::vector<double>& x) { return gradientWrapper(x, gradParsers); },
    //     params
    // );

   // - using the overriden function with muParser objects
    std::vector<double> res = armijoOptimization<stepStrategy, gradStrategy>(p, gradParsers, params);

    // Print
    std::cout << "Argmin: x1 = " << res[0] << ", x2 = " << res[1] << std::endl;

    return 0;
}
