from sympy import * 
from sympy.abc import x 
from sympy.parsing import parse_expr 
from sympy.parsing.latex import parse_latex

from solvers.sys_solvers.solve_linear import *

def solveReducibleToLinear(odeString, user_type):  
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
          n = Add(subTerm.args[1], Integer(1))
          subG = Mul(term, Pow(subTerm, Integer(-1)))
          functionG = Add(functionG, subG)
    else:
      functionF = Add(functionF, term)
  
  step = []
  step.append("- Identify the reducible to linear equation, its parts and degree" + "\\\\ \\\\")
  subSteps = []
  h0 = "From the equation its degree is given by: " + n  + "\\\\ \\\\"
  subSteps.append(h0)

  functionF2 = Mul(functionF, Integer(Add(Integer(1), Mul(Integer(-1), n))))
  functionG2 = Mul(functionG, Integer(Add(Integer(1), Mul(Integer(-1), n))))

  u = Function('u')
  h1 = "Finding the right substitution in the parameter u(x)" + "\\\\ \\\\ "
  h2 = "Substituting into the equation, yields " + "\\\\ \\\\ "
  equation = Eq(Add(Derivative(u(x), x), Mul(functionF2, u(x))), functionG2)
  h3 = "$" + latex(equation) + "$" + "\\\\ \\\\ "
  h4 = "Which is linear" + "\\\\ \\\\ "
  
  subSteps.append(h1)
  subSteps.append(h2)
  subSteps.append(h3)
  subSteps.append(h4)
  step.append(subSteps)
  solveArray.append(step)

  odeStringEqLeft = equation.args[0]
  odeStringEqRigth = equation.args[1]
  odeStringLinear = str(odeStringEqLeft) + "=" + str(odeStringEqRigth)

  solveFromLinear = solveLinear(odeStringLinear, 'u', user_type)

  solveArray += solveFromLinear[1]

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
