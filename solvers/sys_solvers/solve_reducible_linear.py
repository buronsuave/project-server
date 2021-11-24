from sympy import * 
from sympy.abc import x 
from sympy.parsing import parse_expr 
from sympy.parsing.latex import parse_latex

from solvers.sys_solvers.solve_linear import *

def solveReducibleToLinear(odeString, user_type):  
  try:
    odeLeftString = odeString.split("=")[0]
    odeRightString = odeString.split("=")[1]

    odeLeftSym = parse_expr(odeLeftString)
    odeRightSym = parse_expr(odeRightString)

    y = Function('y')
    equation = Eq(odeLeftSym - odeRightSym, 0)

    solveArray = []

    left = equation.args[0]
    exp = solve(left, Derivative(y(x), x))
    aux = expand(exp[0])

    left = Derivative(y(x), x)

    functionF = parse_expr("0")
    functionG = parse_expr("0")

    n = Integer(0)

    aux = Mul(aux, Pow(y(x), Integer(-1)))
    aux = simplify(aux)

    for term in aux.args:
      if 'y' in str(term):
        for subTerm in term.args:
          if 'y' in str(subTerm):
            if type(subTerm) is Pow:
              n = Add(subTerm.args[1], Integer(1))
              subG = Mul(term, Pow(subTerm, Integer(-1)))
              functionG = Add(functionG, subG)
            else:
              n = 2
              subG = Mul(term, Pow(subTerm, Integer(-1)))
              functionG = Add(functionG, subG)
      else:
        functionF = Add(functionF, term)
    
    print(functionF)
    print(functionG)

    step = []
    step.append("- Identify the reducible to linear equation, its parts and degree" + "\\\\ \\\\")
    subSteps = []
    h0 = "From the equation its degree is given by: " + str(n)  + "\\\\ \\\\"
    subSteps.append(h0)

    h0 = "Dividing the original equation by y raised to " + str(n)  + ": \\\\ \\\\"
    subSteps.append(h0)

    newEquation = simplify(Mul(odeLeftSym, Pow(Function('y')(x), -1*n)))
    e0 = "$" + latex(newEquation) +" = 0$" + "\\\\ \\\\"
    subSteps.append(e0)

    functionF2 = Mul(functionF, Integer(Add(Integer(1), Mul(Integer(-1), n))))
    functionG2 = Mul(functionG, Integer(Add(Integer(1), Mul(Integer(-1), n))))
    print(functionF2)
    print(functionG2)

    u = Function('u')

    h1 = "Finding the right substitution in the parameter u(x)" + "\\\\ \\\\ "
    subSteps.append(h1)

    e1 = "$" + "u(x) = y^{1-n}" +"$" + "\\\\ \\\\"
    subSteps.append(e1)

    e2 = "$" + "\\frac{1}{1-n}\\frac{du}{dx} = \\frac{dy}{dx}y^{-n}" + "$" + "\\\\ \\\\"
    subSteps.append(e2)

    haux = "Applying with n = " + str(n)
    subSteps.append(haux)

    e3 = "$" + "u(x) = y^{" + str(1-n) + "}$" + "\\\\ \\\\"
    subSteps.append(e3)

    e4 = "$" + "\\frac{1}{" + str(1-n) + "}\\frac{du}{dx} = \\frac{dy}{dx}y^{" + str(-1*n) + "}" + "$" + "\\\\ \\\\"
    subSteps.append(e4)

    h2 = "Substituting into the equation, yields " + "\\\\ \\\\ "
    subSteps.append(h2)

    equation = Eq(Add(Derivative(u(x), x), Mul(functionF2, u(x))), -1*functionG2)

    h3 = "$" + latex(equation) + "$" + "\\\\ \\\\ "
    subSteps.append(h3)

    h4 = "Which is linear" + "\\\\ \\\\ "
    subSteps.append(h4)

    step.append(subSteps)
    solveArray.append(step)

    odeStringEqLeft = equation.args[0]
    odeStringEqRigth = equation.args[1]
    odeStringLinear = str(odeStringEqLeft) + "=" + str(odeStringEqRigth)

    solveFromLinear = solveLinear(odeStringLinear, 'u', user_type)
    solveArray += solveFromLinear[1]
    solveForU = solveFromLinear[2]
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
            finalSolve.append(Pow(singleSolveForU, 1/(1-n)))          

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
          subSteps.append("Can not get the explicit solution solving for y" + "\\\\ \\\\")

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
