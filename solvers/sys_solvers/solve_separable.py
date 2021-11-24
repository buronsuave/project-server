from timers.custom_threads import PropagatingThread
from anomalies.completeness_anomaly import CompletenessAnomaly
from sympy import * 
from sympy.abc import x 
from sympy.parsing import parse_expr 

from algebraics.operations import *
from integrals.integrator import *
from analytics.investigator import *

def solveSeparable(odeString, functionName, user_type): 
  '''
    ------------------------------------------------------
    # Init solve
    ------------------------------------------------------
  '''
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

    equation = Eq(alg_subs(odeLeftSym, odeRightSym), 0)
    left = equation.args[0]
    express = alg_solve(Eq(left, Integer(0)), Derivative(y(x), x))
    aux = alg_simplify(express[0])
    aux = alg_expand(aux)
    aux = alg_factor(aux)

    '''
    ------------------------------------------------------
    # Step 01: Detect separable structure
    ------------------------------------------------------
    '''
    solveArray.append([])
    step = solveArray[0]
    step.append("- Identify the separable equation and its parts" + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]

    functionF = Integer(1)
    functionG = Integer(1)

    if type(aux) is Add:
      if functionName in str(aux):
        functionG = aux
      else:
        functionF = aux
    else:
      if not (functionName in str(aux)):
        functionF = aux
      else:
        for term in aux.args:
          if functionName in str(term):
            functionG = alg_mul(functionG, term)
          else:
            functionF = alg_mul(functionF, term)

    functionG = alg_mul_inv(functionG)
    expr = alg_div(functionF, functionG)

    h0 = "With algebra, transform the expression: " + "\\\\ \\\\"
    subSteps.append(h0)

    eq0 = "$" + latex(Eq(odeLeftSym, odeRightSym)) + "$" + "\\\\ \\\\"
    subSteps.append(eq0)

    h1 = "into the equation: " + "\\\\ \\\\"
    subSteps.append(h1)

    eq1 = "$" + latex(Derivative(y(x), x)) + " = " + latex(expr) + "$" + "\\\\ \\\\"
    subSteps.append(eq1)
    
    h2 = "wich has the form: " + "\\\\ \\\\"
    subSteps.append(h2)

    eq2 = "$" + latex(Derivative(y(x), x)) + " = " + latex(alg_div(Function('f')(x), \
      Function('g')(Symbol(functionName)))) + "$" + "\\\\ \\\\"
    subSteps.append(eq2)
    
    h3 = "where: " + "\\\\ \\\\"
    subSteps.append(h3)

    eq3 = "$g{\\left("+functionName+"\\right)} = " + latex(functionG) + "$ \\\\ \\\\"
    subSteps.append(eq3)

    eq4 = "$f{\\left(x \\right)} = " + latex(functionF) + "$ \\\\ \\\\"
    subSteps.append(eq4)

    h4 = "So, it is 1st order separable" + "\\\\ \\\\"
    subSteps.append(h4)    

    '''
    ------------------------------------------------------
    # Step 02: Separate functions
    ------------------------------------------------------
    '''
    solveArray.append([])
    step = solveArray[1]
    step.append("- Separate functions" + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]

    functionG = alg_substitution(functionG, y(x), Symbol(functionName))
    left = alg_mul(functionG, Symbol('(d'+functionName+')'))
    right = alg_mul(functionF, Symbol('(dx)'))
    
    h0 = "Multiply by the differential of x and multiply by " + \
      f"g({functionName})" + ", so the result is " + \
      f"g({functionName})" + " and " + "f(x)" + \
      " with their respective differentials" + "\\\\ \\\\"
      
    subSteps.append(h0)

    eq0 = "$" + latex(left) + " = " + latex(right)+ "$" + "\\\\ \\\\"
    subSteps.append(eq0)    
    
    '''
    ------------------------------------------------------
    # Step 03: Integrate Left Side
    ------------------------------------------------------
    '''
    solveArray.append([])
    step = solveArray[2]
    step.append("- Solve left" + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]

    left = alg_div(left, Symbol('(d'+functionName+')'))

    h0 = "Integrate left side with respect to " + functionName + "\\\\ \\\\"
    subSteps.append(h0)

    eq0 = "$\int{" +  latex(left) +"} d"+functionName+"$" 
    subSteps.append(eq0 + "\\\\ \\\\")

    left = alg_expand(left)
    left_int_solve = int_solve(left, Symbol(functionName))
    left = left_int_solve["solution"]

    subSteps.append("-------------------------------" + "\\\\ \\\\")
    for int_substep in left_int_solve["steps"]:
      subSteps.append(int_substep["text"] + "\\\\ \\\\")
      subSteps.append(int_substep["symbol"] + "\\\\ \\\\")
      subSteps.append("-------------------------------" + "\\\\ \\\\")

    '''
    ------------------------------------------------------
    # Step 04: Integrate Right Side
    ------------------------------------------------------
    '''
    solveArray.append([])
    step = solveArray[3]
    step.append("- Solve right" + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]

    h0 = "Integrate right side with respect to x" + "\\\\ \\\\"
    subSteps.append(h0)

    right = alg_div(right, Symbol('(dx)'))
    eq0 = "$\int{" +  latex(right) +"} dx$" 
    subSteps.append(eq0 + "\\\\ \\\\")

    right = alg_expand(right)
    right_int_solve = int_solve(right, x)
    right = right_int_solve["solution"]

    subSteps.append("-------------------------------" + "\\\\ \\\\")
    for int_substep in right_int_solve["steps"]:
      subSteps.append(int_substep["text"] + "\\\\ \\\\")
      subSteps.append(int_substep["symbol"] + "\\\\ \\\\")
      subSteps.append("-------------------------------" + "\\\\ \\\\")

    '''
    ------------------------------------------------------
    # Step 05: Equate Both Sides
    ------------------------------------------------------
    '''
    solveArray.append([])
    step = solveArray[4]
    step.append("- Get implicit solution" + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]

    h1s5 = "Equate both sides" + "\\\\ \\\\"
    subSteps.append(h1s5)

    
    exp1s5 ="$" + latex(left) + " = " + latex(right) + "$" + "\\\\ \\\\"    
    subSteps.append(exp1s5)

    h2s5 = "Substract right side from both sides and add the arbitrary constant C. " + \
      "The implicit answer is: " + "\\\\ \\\\"
    subSteps.append(h2s5)

    express = Add(left, Mul(right, Integer(-1)), Symbol('C'))
    eq1s5 = "$" + latex(express) + "$ = 0"+ "\\\\ \\\\"
    subSteps.append(eq1s5)    

    '''
    ------------------------------------------------------
    # Step 06: Get Explicit Solve
    ------------------------------------------------------
    '''
    solveArray.append([])
    step = solveArray[5]
    step.append("- Get the explicit solution solving for " + functionName + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]

    global finalSolve
    finalSolve = []
    
    def final_solve_timeout(expression, symbol):
      global finalSolve
      finalSolve = solve(expression, symbol)

    try:
      process = PropagatingThread(target = final_solve_timeout, args=(express, Symbol(functionName)))
      process.start()
      process.join(timeout=5)

      for singleSolve in finalSolve:
        eq1s6 = Eq(y(x), singleSolve)
        subSteps.append("$" + latex(eq1s6) + "$" + "\\\\ \\\\") 

        # Analytic intervention for all the single solves if is teacher
        if (user_type == 'teacher' and functionName == 'y'):
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

      if (user_type == "teacher" and functionName == 'y'):
        '''
        ------------------------------------------------------
        # Step 07: Generate Plot
        ------------------------------------------------------
        '''
        solveArray.append([])
        step = solveArray[6]
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

    def display_solve():
        solveStr = ""
        for stepAux in solveArray:
          solveStr += stepAux[0]
          solveStr += display_step(stepAux[1])
        return solveStr    
    return [display_solve(), solveArray, finalSolve]

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




