from sympy import * 
from sympy.abc import x 
from sympy.parsing import parse_expr 
from sympy.parsing.latex import parse_latex

from solvers.sys_solvers.solve_separable import *

def solveHomogeneous(odeString, functionName, user_type):
  
  '''
    ------------------------------------------------------
    # Init solve
    ------------------------------------------------------
  '''
  # Init solve array
  solveArray = []

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

  h0 = "Since it's homogeneous, the derivative can be expressed as a function that is also homogeneous, that is:" + "\\\\ \\\\"
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

  solveForU = solutionSeparable[2]

  '''
  ------------------------------------------------------
  # Step 02: Get Explicit Solve
  ------------------------------------------------------
  '''
  solveArray.append([])
  step = solveArray[1]
  step.append("- Undo the variable change" + "\\\\ \\\\")
  step.append([])
  subSteps = step[1]
  
  global finalSolve
  finalSolve = []

  def final_solve_timeout(expression, symbol):
    global finalSolve
    finalSolve = solve(expression, symbol)


  if (len(solveForU)) > 0:
    h6 = "Multiplying both sides by x and using the change of variable for each solution" + "\\\\ \\\\"
    subSteps.append(h6)
    
    for particularSolveForU in solveForU:
      particularSolve = Eq(y(x), Mul(particularSolveForU, x))
      eq5 = "$" + latex(particularSolve) + "$ \\\\ \\\\"
      subSteps.append(eq5)
    
    step.append(subSteps)
    solveArray.append(step)

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
  return [ display_solve(solveArray), solveArray ]  