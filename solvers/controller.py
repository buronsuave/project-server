from anomalies.classification_anomaly import ClassificationAnomaly
from anomalies.completeness_anomaly import CompletenessAnomaly
from sympy import *
from sympy.abc import x

from parsers.parse_sympy import parseLatex
from solvers.sys_solvers.solve_separable import solveSeparable
from solvers.sys_solvers.solve_linear import solveLinear
from solvers.sys_solvers.solve_homogeneous import solveHomogeneous
from solvers.sys_solvers.solve_exact import solveExact
from solvers.sys_solvers.solve_reducible_linear import solveReducibleToLinear
from solvers.sys_solvers.solve_n_order import solveNLinear
from classficators.classificator import classify

def solve(inputString, user_type):
    # Parse Latex expression 
    try:
        equation = parseLatex(inputString)
    except Exception as e:
        raise e

    # Trying to catch a completeness anomaly 
    try:
        # Classify ODE. Could raise a Classification Anomaly
        odeType = classify(str(equation) + "= 0")
        print("ODE TYPE: " + odeType)
        print()

        # Global difficulty 
        global global_difficulty
        global_difficulty = 0

        if odeType == "separable":
            solveArray = solveSeparable(str(equation) + "= 0", 'y', user_type)
            print("Global Difficulty: " + str(global_difficulty))
            return solveArray[1]

        elif odeType == "linear":
            solveArray = solveLinear(str(equation) + "= 0", 'y', user_type)
            print("Global Difficulty: " + str(global_difficulty))
            return solveArray[1]

        elif odeType == "reducible":
            solveArray = solveReducibleToLinear(str(equation) + "= 0", user_type)
            print("Global Difficulty: " + str(global_difficulty))
            return solveArray[1]

        elif odeType == "homogeneous":
            solveArray = solveHomogeneous(str(equation) + "= 0", 'y', user_type)
            print("Global Difficulty: " + str(global_difficulty))
            return solveArray[1]

        elif odeType == "exact":
            solveArray = solveExact(str(equation) + "= 0", user_type)
            print("Global Difficulty: " + str(global_difficulty))
            return solveArray[1]

        elif odeType == "superior":
            solveArray = solveNLinear(str(equation) + "= 0", "y", user_type)
            print("Global Difficulty: " + str(global_difficulty))
            return solveArray[1]

    # Classification error on classify call
    except ClassificationAnomaly as clsa:
        print("undefined classification")
        try:
            # Launch DSolve intervention for solving 
            # an undefined ODE type on server
            solveSingle = dsolve(Eq(equation, 0), Function('y')(x))

            # Create single step in case of found a solution
            solveArray = []
            step = []
            step.append("- Solve with DSolve (backup system): " +  
            "\\\\ \\\\")
            subSteps = []
            h0 = "The server was not able to build the steps for " + \
            "the solution." + "\\\\ \\\\" + \
            "However, the solution found was the following:" + \
            "\\\\ \\\\"
            eq0 = "$" + latex(solveSingle) + "$" + "\\\\ \\\\"
            subSteps.append(h0)
            subSteps.append(eq0)
            step.append(subSteps)
            solveArray.append(step)
            clsa.set_final_solve(solveArray)
        except Exception as e:
            print(e.args[0])
        finally:
            raise clsa

    # Completeness error on some solve call
    except CompletenessAnomaly as ca:
        print("unsolvable by the system")
        print("Global Difficulty: " + str(global_difficulty))
        try:
            # Launch DSolve intervention for solving an uncompleted ODE on server
            solveSingle = dsolve(Eq(equation, 0), Function('y')(x))
            step = []
            step.append(latex("- Solve with DSolve (backup system): ") + "\\\\ \\\\")
            subSteps = []
            h0 = latex("The server was not able to complete the steps for the solution.") + "\\\\ \\\\" + latex("However, the solution found was the following:") + "\\\\ \\\\"
            eq0 = "$" + latex(solveSingle) + "$" + "\\\\ \\\\"
            subSteps.append(h0)
            subSteps.append(eq0)
            step.append(subSteps)
            ca.append_final_solve(step)
        except Exception as e:
            print(e.args[0])
        finally:
            raise ca
    
    except Exception as ga:
        print("Non controlled Exception")
        print(ga.args[0])
        # Launch DSolve intervention for solving 
        # an undefined ODE type on server
        solveSingle = dsolve(Eq(equation, 0), Function('y')(x))
        clsa = CompletenessAnomaly(None)

        # Create single step in case of found a solution
        solveArray = []
        step = []
        step.append("- Solve with DSolve (backup system): " +  
        "\\\\ \\\\")
        subSteps = []
        h0 = "The server was not able to build the steps for " + \
        "the solution." + "\\\\ \\\\" + \
        "However, the solution found was the following:" + \
        "\\\\ \\\\"
        eq0 = "$" + latex(solveSingle) + "$" + "\\\\ \\\\"
        subSteps.append(h0)
        subSteps.append(eq0)
        step.append(subSteps)
        solveArray.append(step)
        clsa.set_partial_solution(solveArray)

        raise clsa        