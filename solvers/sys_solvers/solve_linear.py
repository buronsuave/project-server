from timers.custom_threads import PropagatingThread
from anomalies.completeness_anomaly import CompletenessAnomaly
from sympy import * 
from sympy.abc import x 
from sympy.parsing import parse_expr 

from algebraics.operations import *
from integrals.integrator import *
from analytics.investigator import *

def solveLinear(odeString, functionName, user_type):
  
  '''
    ------------------------------------------------------
    # Init solve
    ------------------------------------------------------
  '''
  solveArray = []
  
  try: 
    odeLeftString = odeString.split("=")[0]
    odeRightString = odeString.split("=")[1]
    odeLeftSym = parse_expr(odeLeftString)
    odeRightSym = parse_expr(odeRightString)
    y = Function(functionName)
    equation = Eq(odeLeftSym - odeRightSym, 0)
    left = equation.args[0]
    exp = alg_solve(left, Derivative(y(x), x))
    aux = alg_expand(exp[0])

    left = Derivative(y(x), x)

    functionF = parse_expr("0")
    functionG = parse_expr("0")

    for term in aux.args:
      if functionName in str(term):
        functionF = Add(functionF, Mul(term, Pow(y(x), Integer(-1))))
      else:
        functionG = Add(functionG, term)

    functionF = Mul(functionF, Integer(-1))
    functionF = simplify(functionF)
    functionG = simplify(functionG)

    right = alg_add(functionG, Mul(Integer(-1), functionF, y(x)))

    '''
    ------------------------------------------------------
    # Step 01: Identify the linear equation
    ------------------------------------------------------
    '''
    solveArray.append([])
    step = solveArray[0]
    step.append("- Identify the linear equation and its parts" + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]

    h0 = "With algebra, transform the expression: " + "\\\\ \\\\"
    subSteps.append(h0)
    
    eq0 = "$" + latex(Eq(odeLeftSym, odeRightSym)) + "$" + "\\\\ \\\\"
    subSteps.append(eq0)

    h1 = "into the equation: " + "\\\\ \\\\"
    subSteps.append(h1)

    eq1 = "$" + latex(Derivative(y(x), x)) + " = " + latex(exp) + "$" + "\\\\ \\\\"
    subSteps.append(eq1)

    h2 = "wich has the form: " + "\\\\ \\\\"
    subSteps.append(h2)

    eq2 = "$" + latex(Derivative(y(x), x)) + " = " + latex(Add(Function('g')(x), Mul(Function('f')(x), Integer(-1), y(x)))) + "$" + "\\\\ \\\\"
    subSteps.append(eq2)

    h3 = "where: " + "\\\\ \\\\"
    subSteps.append(h3)

    eq3 = "$g{\\left(x \\right)} = " + latex(functionG) + "$ \\\\ \\\\"
    subSteps.append(eq3)

    eq4 = "$f{\\left(x \\right)} = " + latex(functionF) + "$ \\\\ \\\\"
    subSteps.append(eq4)

    h4 = "So, it is 1st order linear" + "\\\\ \\\\"
    subSteps.append(h4)

    '''
    ------------------------------------------------------
    # Step 02: Calculate integral factor
    ------------------------------------------------------
    '''
    solveArray.append([])
    step = solveArray[1]
    step.append("- Calculate integral factor" + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]

    h1s2 = "Lets propose a function M(x) such that: " + "\\\\ \\\\"
    subSteps.append(h1s2)

    eqAux = Eq(Mul(Function('M')(x), Function('f')(x)), Derivative(Function('M')(x), x))
    eq1s2 = "$" +  latex(eqAux) + "$" + "\\\\ \\\\"
    subSteps.append(eq1s2)

    h2s2 = "Substituting: " + "\\\\ \\\\"
    subSteps.append(h2s2)
    
    eq2s2 = "$" + latex(Eq(Mul(Function('M')(x), functionF), Derivative(Function('M')(x), x))) + "$" + "\\\\ \\\\"
    subSteps.append(eq2s2)
    
    h3s2 = "Which is a 1st order separable diferential equation. Hence solving for M(x)" + "\\\\ \\\\"
    subSteps.append(h3s2)
    
    functionM = Pow(E, Integral(functionF, x))
    eq3s2 = "$" + latex(Eq(Mul(Pow(Function('M')(x), Integer(-1)), Symbol('dM(x)')), Mul(functionF, Symbol('(dx)')))) + "$" + "\\\\ \\\\"
    subSteps.append(eq3s2)
    
    eq4s2 = "$" + latex(Eq(log(Function('M')(x)), Integral(functionF, x))) + "$" + "\\\\ \\\\"
    subSteps.append(eq4s2)    
    
    functionF = expand(functionF)
    
    integralM = int_solve(expand(functionF), x)
    exponentM = integralM["solution"]
    subSteps.append("-------------------------------" + "\\\\ \\\\")
    for int_substep in integralM["steps"]:
      subSteps.append(int_substep["text"] + "\\\\ \\\\")
      subSteps.append(int_substep["symbol"] + "\\\\ \\\\")
      subSteps.append("-------------------------------" + "\\\\ \\\\")

    eq5s2 = "$" + latex(Eq(log(Function('M')(x)), exponentM)) + "$" + "\\\\ \\\\"
    subSteps.append(eq5s2)    
    
    eq6s2 = "$" + latex(Function('M')(x)) + " = " + latex(Pow(E, exponentM))+ "$" + "\\\\ \\\\"
    subSteps.append(eq6s2)    

    '''
    ------------------------------------------------------
    # Step 03: Reduce to 1st order separable ODE
    ------------------------------------------------------
    '''
    solveArray.append([])
    step = solveArray[2]
    step.append("- Reduce to 1st Order Separable ODE" + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]

    functionM = functionM.replace(Integral(functionF, x), exponentM)
    functionM = Pow(E, exponentM)
    functionM = alg_simplify(functionM)

    h1s3 = "Multiplying the original equation by M(x):" + "\\\\ \\\\"
    subSteps.append(h1s3)
    
    equation = Eq(left, right)
    left = alg_mul(left, functionM)
    right = alg_add(Mul(Integer(-1), functionF, y(x), functionM), Mul(functionG, functionM))
    equation = Eq(left, right)

    left = alg_add(left, Mul(functionF, y(x), functionM))
    right = alg_add(right, Mul(functionF, y(x), functionM))
    equation = Eq(left, right)
    equationaux = Eq(expand(Mul( Function('M')(x), left, pow(functionM, Integer(-1)))), Mul( Function('M')(x), right, pow(functionM, Integer(-1))))

    eq1s3 = "$" + latex(equationaux) + "$" + "\\\\ \\\\"
    subSteps.append(eq1s3)
    
    h2s3 = "By the definition of M(x) this is equivalent to: " +  "\\\\ \\\\"
    subSteps.append(h2s3)
    
    eq2s3 = "$" + latex(Add(Mul(Derivative(y(x),x), Function('M')(x)),Mul(y(x), Derivative(Function('M')(x), x)))) + " = " + latex(Mul(functionG, Function('M')(x))) +  "$" + "\\\\ \\\\"
    subSteps.append(eq2s3)
    
    h3s3 = "Notice that the left hand side can be reduce by the chain rule to: " + "\\\\ \\\\"
    subSteps.append(h3s3)
    
    eq3s3 = "$" + latex(Add(Mul(Derivative(y(x),x), Function('M')(x)),Mul(y(x), Derivative(Function('M')(x), x)))) + " = " + latex(Derivative(Mul(y(x), Function('M')(x)),x)) +  "$" + "\\\\ \\\\"
    subSteps.append(eq3s3)  
    
    h4s3 = "Therefore: " + "\\\\ \\\\"
    subSteps.append(h4s3)
    
    eq4s3 = "$" + latex(Derivative(Mul(y(x), Function('M')(x)),x)) + " = " + latex(Mul(functionG, Function('M')(x))) +  "$" + "\\\\ \\\\"
    subSteps.append(eq4s3)
    
    h5s3 = "Which again is a 1st order separable diferential equation" + "\\\\ \\\\"
    subSteps.append(h5s3)   
    
    equationaux = Eq(Symbol('dM(x)'+functionName+'(x)'), factor(Mul(Symbol('(dx)'),functionG, Function('M')(x))))
    eq5s3 = "$" + latex(equationaux)+  "$" + "\\\\ \\\\"
    subSteps.append(eq5s3)

  
    '''
    ------------------------------------------------------
    # Step 04: Get implicit solve
    ------------------------------------------------------
    '''
    solveArray.append([])
    step = solveArray[3]
    step.append("- Get implicit solution" + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]
    
    left = Derivative(Mul(functionM, y(x)), x)

    equation = Eq(left, right)
    left = alg_mul(left, Pow(Derivative(Mul(functionM, y(x)), x), Integer(-1)), Symbol('d'), Mul(y(x), functionM))
    right = alg_mul(right, Symbol('dx'))
    equation = Eq(left, right)
    h6s3 = "Integrating the left hand side, and indicating the integral at right hand side: " + "\\\\ \\\\"
    subSteps.append(h6s3)   
 
    left = alg_mul(y(x), functionM)
    right = alg_mul(right, Pow(Symbol('dx'), Integer(-1)))
    eq6s3 = "$" + latex(Symbol('M(x)'+functionName+'(x)'))+ " = " + latex(Integral(Mul(Mul(right, Pow(functionM, Integer(-1))),Function('M')(x)),x)) + "$" + "\\\\ \\\\"
    subSteps.append(eq6s3)
 
    h7s3 = "Substituting M(x) = " + "\\\\ \\\\"
    subSteps.append(h7s3)
    
    eqAux1 = "$" + latex(functionM) + "$" +  "\\\\ \\\\"
    subSteps.append(eqAux1) 
    
    equationaux = Eq(left, Integral(right,x))
    eq7s3 = "$" + latex(equationaux)+ "$" + "\\\\ \\\\"
    subSteps.append(eq7s3)
    
    right = expand(right)

    right_integral = int_solve(right, x)
    right = right_integral["solution"]
    subSteps.append("-------------------------------" + "\\\\ \\\\")
    for int_substep in right_integral["steps"]:
      subSteps.append(int_substep["text"] + "\\\\ \\\\")
      subSteps.append(int_substep["symbol"] + "\\\\ \\\\")
      subSteps.append("-------------------------------" + "\\\\ \\\\")

    right = alg_add(right, Symbol('C'))
    equation = Eq(left, right)
    
    h8s3 = "Integrating the right hand side: " + "\\\\ \\\\"
    subSteps.append(h8s3) 
    
    eq8s3 = "$" + latex(equation)+ "$" + "\\\\ \\\\"
    subSteps.append(eq8s3)

    '''
    ------------------------------------------------------
    # Step 05: Obtain solution
    ------------------------------------------------------
    '''
    solveArray.append([])
    step = solveArray[4]
    step.append("- Solve for " + functionName + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]

    left = y(x)
    right = alg_mul(right, Pow(functionM, Integer(-1)))
    right = alg_simplify(right)
    equation = Eq(left, right)
    h9s3 = "Solve for " + "\\\\ \\\\"
    subSteps.append(h9s3) 

    eqAux = latex(Symbol(functionName+'(x)')) + + "\\\\ \\\\"
    eq9s3 = "$" + latex(equation)+ "$" + "\\\\ \\\\"
    subSteps.append(eq9s3)
    
    global finalSolve
    finalSolve = []
    
    def final_solve_timeout(expression, symbol):
      global finalSolve
      finalSolve = solve(expression, symbol)

    try:
      process = PropagatingThread(target = final_solve_timeout, args=(equation, Symbol(functionName)))
      process.start()
      process.join(timeout=5)

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

    def display_solve(solveArray):
      solveStr = ""
      for stepAux in solveArray:
        solveStr += stepAux[0]
        solveStr += display_step(stepAux[1])
      return solveStr    
    return [display_solve(solveArray), solveArray]

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
