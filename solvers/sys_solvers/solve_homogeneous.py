from sympy import * 
from sympy.abc import x 
from sympy.parsing import parse_expr 
from solvers.sys_solvers.solve_separable import *

def solveHomogeneous(odeString, functionName, user_type):
  
  '''
    ------------------------------------------------------
    # Init solve
    ------------------------------------------------------
  '''
  # Init solve array
  solveArray = []

  try:
  
    '''
    ------------------------------------------------------
    # Initial algebraic analysis
    ------------------------------------------------------
    '''


    odeLeftString = odeString.split("=")[0]
    odeRightString = odeString.split("=")[1]

    odeLeftSym = parse_expr(odeLeftString)
    odeRightSym = parse_expr(odeRightString)

    y = Function(functionName)
    equation = Eq(odeLeftSym - odeRightSym, 0)

    # Step 1
    left = equation.args[0]
    exp = solve(left, Derivative(y(x), x))
    aux = expand(exp[0])

    left = Derivative(y(x), x)
    # Define the change of variable
    functionF = aux


    '''
    ------------------------------------------------------
    # Step 01: Propose the appropriate variable change to reduce to separable
    ------------------------------------------------------
    '''
    solveArray.append([])
    step = solveArray[0]
    step.append("- Propose the appropriate variable change to reduce to separable" + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]

    #Define function u(x)
    u = Function('u')

    h0 = "Since it is homogeneous, the derivative can be expressed as a function that is also homogeneous, that is:" + "\\\\ \\\\"
    subSteps.append(h0)
    
    eq0 = "$" + latex(Eq(Derivative(y(x), x), functionF)) + "$" + "\\\\ \\\\"
    subSteps.append(eq0)

    #Perform the substitution
    functionF = functionF.subs(y(x), Mul(u(x), x))
    left = Add(Mul(Derivative(u(x), x), x), u(x))

    h1 = "Using the change of variable: " + "\\\\ \\\\"
    subSteps.append(h1)
    
    eq1 = "$" + latex(Eq(y(x), Mul(u(x), x))) + "$" + "\\\\ \\\\"
    subSteps.append(eq1)  

    h2 = "Whose derivative is: " + "\\\\ \\\\"
    subSteps.append(h2)
    
    eq2 = "$" + latex(Eq(Derivative(y(x), x), Add(u(x), Mul(Derivative(u(x), x))))) + "$" + "\\\\ \\\\"
    subSteps.append(eq2)
    
    h3 = "Carrying out the changes for the function and its derivative in the original equation: " + "\\\\ \\\\"
    subSteps.append(h3)

    eq3 = "$" + latex(Eq(left, functionF)) + "$" + "\\\\ \\\\"
    subSteps.append(eq3)

    left = Add(left, Mul(functionF, Integer(-1)))
    left = expand(left)
    separableODE = Eq(left, Integer(0))

    h4 = "Simplifying: " + "\\\\ \\\\"
    subSteps.append(h4)
    
    eq4 = "$" + latex(separableODE) + "$" + "\\\\ \\\\"
    subSteps.append(eq4)
    
    h5 = "Wich is first order separable" + "\\\\ \\\\"
    subSteps.append(h5)
    
    solutionSeparable = solveSeparable(str(separableODE.args[0]) + "= 0", 'u', user_type)
    solveArray += solutionSeparable[1]

    #Final Solution
    solveForU = solutionSeparable[2]
    print(solveForU)

    '''
    ------------------------------------------------------
    # Step 02: Get Explicit Solve
    ------------------------------------------------------
    '''
    solveArray.append([])
    step = solveArray[len(solveArray) - 1]
    step.append("- Undo the variable change" + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]
    
    global finalSolve
    finalSolve = []

    if (len(solveForU)) > 0:
        try:
          for singleSolveForU in solveForU:
            finalSolve.append(alg_mul(singleSolveForU, x))          

          for singleSolve in finalSolve:
            eq1s6 = Eq(y(x), singleSolve)
            subSteps.append("$" + latex(eq1s6) + "$" + "\\\\ \\\\") 

            # Analytic intervention for all the single solves if is teacher
            if (user_type == 'teacher'):
              print("Teacher")
              try:
                roots = []
                roots_process = PropagatingThread(target = get_roots, args = [singleSolve, roots])
                roots_process.start()
                roots_process.join(timeout = 3)

                h0 = "Whose roots are: " + "\\\\ \\\\"
                subSteps.append(h0) 
                subIndex = 1
                for root in roots:
                  eq0 = "$" + "x_{" + str(subIndex) + "} = " + latex(root) + "$" + "\\\\ \\\\"
                  subIndex = subIndex + 1
                  subSteps.append(eq0)

              except Exception as e:
                print("Error with roots")
                print(e)

              try:
                critics = []
                critics_process = PropagatingThread(target = max_min, args = [singleSolve, critics])
                critics_process.start()
                critics_process.join(timeout = 3)

                h0 = "Whose critics are: " + "\\\\ \\\\"
                subSteps.append(h0)
                subIndex = 1
                for critic in critics:
                  eq0 = "$" + "x_{" + str(subIndex) + "} = " + latex(critic) + "$" + "\\\\ \\\\"
                  subIndex = subIndex + 1
                  subSteps.append(eq0)

              except Exception as e:
                print("Error with critics")
                print(e)

              try:
                inflexions = []
                inflexions_process = PropagatingThread(target = inflexion_points, args = [singleSolve, inflexions])
                inflexions_process.start()
                inflexions_process.join(timeout = 3)

                h0 = "Whose inflexions are: " + "\\\\ \\\\"
                subSteps.append(h0)
                subIndex = 1
                for inflexion in inflexions:
                  eq0 = "$" + "x_{" + str(subIndex) + "} = " + latex(inflexion) + "$" + "\\\\ \\\\"
                  subIndex = subIndex + 1
                  subSteps.append(eq0)

              except Exception as e:
                print("Error with inflexions")
                print(e)

          if (user_type == "teacher"):
            '''
            ------------------------------------------------------
            # Step 07: Generate Plot
            ------------------------------------------------------
            '''
            solveArray.append([])
            step = solveArray[len(solveArray) - 1]
            step.append("- Graphs" + "\\\\ \\\\")
            step.append([])
            subSteps = step[1]

            for singleSolve in finalSolve:
              # Add plot step to solution
              print("Creating plot")

              try:
                plot_string = create_plot(singleSolve)[1:]
                plot_string = plot_string.replace("\\n", "")
              except Exception as e:
                print(e)

              subSteps.append(plot_string)
              print("Plot appended")      
          
        except:
          subSteps.append("Can not get the explicit solution solving for " + functionName + "\\\\ \\\\")

    def display_step(step):
      stepStr = ""
      for subStep in step:
        stepStr += str(subStep)
      return stepStr

    def display_solve(solveArray):
      solveStr = ""
      for stepAux in solveArray:
        if len(stepAux) != 0:
          solveStr += stepAux[0]
          solveStr += display_step(stepAux[1])
        else:
          solveArray.remove(stepAux)
      return solveStr    

    return [ display_solve(solveArray), solveArray ]  
  
  except CompletenessAnomaly as ca:
    
    if ca.partial_solution[0][0] == "partial integral":
      step = solveArray[len(solveArray) - 1]
      subSteps = step[1]
      subSteps.append("-------------------------------" + "\\\\ \\\\")

      for int_substep in ca.partial_solution[0][1]:
        subSteps.append(int_substep["text"] + "\\\\ \\\\")
        subSteps.append(int_substep["symbol"] + "\\\\ \\\\")
        subSteps.append("-------------------------------" + "\\\\ \\\\")

    ca.set_partial_solution(solveArray)

    raise ca