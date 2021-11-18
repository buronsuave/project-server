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
  step.append("- Identify the separable equation and its parts" + "\\\\ \\\\")
  step.append([])
  subSteps = step[1]


  h0 = "Since it's homogeneous, the derivative can be expressed as a function that is also homogeneous, that is:" + "\\\\ \\\\"
  eq0 = "$" + latex(Eq(Derivative(y(x), x), functionF)) + "$" + "\\\\ \\\\"

  u = Function('u')

  functionF = functionF.subs(y(x), Mul(u(x), x))
  left = Add(Mul(Derivative(u(x), x), x), u(x))
  
  h3 = "Carrying out the changes for the function and its derivative in the original equation: " + "\\\\ \\\\"
  eq3 = "$" + latex(Eq(left, functionF)) + "$" + "\\\\ \\\\"

  left = Add(left, Mul(functionF, Integer(-1)))
  left = expand(left)
  separableODE = Eq(left, Integer(0))

  h4 = "Simplifying: " + "\\\\ \\\\"
  eq4 = "$" + latex(separableODE) + "$" + "\\\\ \\\\"
  h5 = "Wich is first order separable" + "\\\\ \\\\"

  h1 = "Using the change of variable: " + "\\\\ \\\\"
  eq1 = "$" + latex(Eq(y(x), Mul(u(x), x))) + "$" + "\\\\ \\\\"

  h2 = "Whose derivative is: " + "\\\\ \\\\"
  eq2 = "$" + latex(Eq(Derivative(y(x), x), Add(u(x), Mul(Derivative(u(x), x))))) + "$" + "\\\\ \\\\"

  step = []
  subSteps = []
  subSteps.append(h0)
  subSteps.append(eq0)
  subSteps.append(h1)
  subSteps.append(eq1)
  subSteps.append(h2)
  subSteps.append(eq2)
  subSteps.append(h3)
  subSteps.append(eq3)
  subSteps.append(h4)
  subSteps.append(eq4)
  subSteps.append(h5)
  step.append(subSteps)
  solveArray.append(step)
    
  solutionSeparable = solveSeparable(str(separableODE.args[0]) + "= 0", 'u', user_type)
  solveArray += solutionSeparable[1]

  solveForU = solutionSeparable[2]
  if (len(solveForU)) > 0:
    h6 = "Multiplying both sides by x and using the change of variable for each solution" + "\\\\ \\\\"
    step = []
    step.append("- Undo the variable change" + "\\\\ \\\\")
    subSteps = []
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