from sympy import * 
from sympy.abc import x 
from sympy.parsing import parse_expr 
from sympy.parsing.latex import parse_latex
from anomalies.completeness_anomaly import CompletenessAnomaly

from algebraics.operations import *
from integrals.integrator import *
from analytics.investigator import *

def solveExact(odeString, functionName, user_type):  
  
  '''
  ------------------------------------------------------
  # Init solve
  ------------------------------------------------------
  '''
  try:
    
    '''
    ------------------------------------------------------
    # Initial algebraic analysis
    ------------------------------------------------------
    '''
    # init solve array
    solveArray = [] 

    odeLeftString = odeString.split("=")[0]
    odeRightString = odeString.split("=")[1]

    odeLeftSym = parse_expr(odeLeftString)
    odeRightSym = parse_expr(odeRightString)


    y = Function('y')
    equation = Eq(odeLeftSym - odeRightSym, 0)
    equation = equation.subs(y(x), Symbol('y'))


    functionP = Integer(0)
    functionQ = Integer(0)
    functionF = Integer(0)
    leftPartial = Integer(0)
    functionG = Function('g')
  
    '''
    ------------------------------------------------------
    # Step 01: Detect exact structure
    ------------------------------------------------------
    '''
    solveArray.append([])
    step = solveArray[0]
    step.append("- Identify the exact equation and its parts" + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]

    for term in equation.args[0].args:
      if 'Derivative' in str(term):
        functionQ = Add(functionQ, Mul(term, Pow(Derivative(Symbol('y'), x), Integer(-1))))
      else:
        functionP = Add(functionP, term)

    h0 = "With algebra, transform the expression: " + "\\\\ \\\\"
    subSteps.append(h0)
    eq0 = "$" + latex(Eq(odeLeftSym, odeRightSym)) + "$" + "\\\\ \\\\"
    subSteps.append(eq0)
    h1 = "into the equation: " + "\\\\ \\\\"
    subSteps.append(h1)
    eq1 = "$" + latex(equation) + "$" + "\\\\ \\\\"
    subSteps.append(eq1)
    h2 = "wich has the form: " + "\\\\ \\\\"
    subSteps.append(h2)
    eq2 = "$" + latex(Function('P')(Symbol('x,y'))) + " + " + latex( Mul(Function('Q')(Symbol('x,y')), Derivative(y(x),x))) + " = 0" +  "$" + "\\\\ \\\\"
    subSteps.append(eq2)
    h3 = "where: " + "\\\\ \\\\"
    subSteps.append(h3)
    eq3 = "$" + latex(Function('Q')(Symbol('x,y'))) + " = " + latex(functionQ) + "$ \\\\ \\\\"
    subSteps.append(eq3)
    eq4 = "$" + latex(Function('P')(Symbol('x,y'))) + " = " + latex(functionP) + "$ \\\\ \\\\"
    subSteps.append(eq4)
    h4 = "So, it is an exact differential equation" + "\\\\ \\\\"
    subSteps.append(h4)    
  

    #Step 2
    '''
    ------------------------------------------------------
    # Step 02: Obtain F(x,y) as a result of P(x,y)
    ------------------------------------------------------
    '''
    solveArray.append([])
    step = solveArray[1]
    step.append("- Obtain F(x,y) as a result of P(x,y) " + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]

    functionF = integrate(functionP, x)
    functionF = Add(functionF, functionG(Symbol('y')))

    h1s2 = "Lets integrate the function " + "\\\\ \\\\"
    subSteps.append(h1s2)
    
    eqAux1 = "$" + latex(Function('P')(Symbol('x'), Symbol('y'))) + "$" + "\\\\ \\\\" 
    subSteps.append(eqAux1)
    
    eqAux = Eq(Function('F')(Symbol('x'), Symbol('y')), Integral(Function('P')(Symbol('x'), Symbol('y')),x))
    eq1s2 = "$" +  latex(eqAux) + "$" + "\\\\ \\\\"
    subSteps.append(eq1s2)
    
    h2s2 = "Substituting " + "\\\\ \\\\"
    subSteps.append(h2s2)
    
    eqAux2 = "$" + latex(Function('P')(Symbol('x'), Symbol('y'))) + " = " + latex(functionP) + "$"  + "\\\\ \\\\"  
    subSteps.append(eqAux2)
    
    eqAux = Eq(Function('F')(Symbol('x'), Symbol('y')), Integral(functionP,x))
    eq2s2 = "$" + latex(eqAux) + "$" + "\\\\ \\\\"
    subSteps.append(eq2s2)
    
    h3s2 = latex("Integrating: ")  + "\\\\ \\\\"
    subSteps.append(h3s2)

    eqAux = Eq(Function('F')(Symbol('x'), Symbol('y')), functionF)
    eq3s2 = "$" + latex(eqAux) + "$" + "\\\\ \\\\"
    subSteps.append(eq3s2)


    #Step 3
    '''
    ------------------------------------------------------
    # Step 03: Use the properties of the exact equation to find (y)
    ------------------------------------------------------
    '''
    solveArray.append([])
    step = solveArray[2]
    step.append(latex("- Use the properties of the exact equation to find (y)") + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]
    
    partialF = diff(functionF, Symbol('y'))
    leftPartial = Add(partialF, Mul(functionQ, Integer(-1)))
    rightGSolveSide = solve(leftPartial, Derivative(functionG(Symbol('y')), Symbol('y')))

    h1s3 = "To find g(y) differentiate respect y the function: " +  "\\\\ \\\\"
    subSteps.append(h1s3)

    eq1s3 = "$" +  latex(Derivative(Function('F')(Symbol('x'), Symbol('y')), Symbol('y')))  +  " = " + latex(partialF) + "$" + "\\\\ \\\\"
    subSteps.append(eq1s3)

    eqAux = "$" + latex(Function('F')(Symbol('y'), Symbol('x'))) +  "$" + "\\\\ \\\\"
    subSteps.append(eqAux)

    h2s3 = "This by definition of a exact diferential equation must be equal to " + "\\\\ \\\\"
    subSteps.append(h2s3)

    eqAux1 = "$" + latex(Function('Q')(Symbol('x'), Symbol('y'))) + "$" + "\\\\ \\\\"
    subSteps.append(eqAux1)
    
    h3s3 = "Hence equating: "  + "\\\\ \\\\"
    subSteps.append(h3s3)

    eq2s3 = "$" +  latex(partialF)  +  " = " + latex(functionQ) + "$" + "\\\\ \\\\"
    subSteps.append(eq2s3)

    h4s3 = "Solving for: " + "\\\\ \\\\"
    subSteps.append(h4s3) 
    
    eqAux2 =  "$" + latex(Derivative(Function('g')(Symbol('y')), Symbol('y'))) +  "$"  + "\\\\ \\\\"
    subSteps.append(eqAux2)
    
    eq3s3 = "$" +  latex(Derivative(Function('g')(Symbol('y')), Symbol('y')))  +  " = " + latex(rightGSolveSide[0]) + "$" + "\\\\ \\\\"
    subSteps.append(eq3s3)

    #Step 4
    '''
    ------------------------------------------------------
    # Step 04: Get g(y) and particular F(x,y)
    ------------------------------------------------------
    '''
    solveArray.append([])
    step = solveArray[3]
    step.append("- Get g(y) and particular F(x,y) " + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]


    gIntValue = integrate(rightGSolveSide[0], Symbol('y'))
    functionF = Add(functionF, Mul(functionG(Symbol('y')), Integer(-1)), gIntValue)

    h1s4 = "Integrating both sides to get " + "$" + latex(Function('g')(Symbol('y'))) +  "$" + "\\\\ \\\\"
    subSteps.append(h1s4)

    eq1s4 = "$" + latex(Function('g')(Symbol('y'))) +  " = " + latex(Integral(rightGSolveSide[0], Symbol('y'))) + "$" + "\\\\ \\\\"
    subSteps.append(eq1s4)

    eq2s4 = "$" + latex(Function('g')(Symbol('y'))) +  " = " + latex(gIntValue) + "$" + "\\\\ \\\\"
    subSteps.append(eq2s4)

    h2s4 = "Substituting this result into F(x,y)" +  "\\\\ \\\\"
    subSteps.append(h2s4)

    eq3s4 = "$" + latex(Function('F')(Symbol('x,y'))) +  " = " + latex(functionF) + "$" + "\\\\ \\\\"
    subSteps.append(eq3s4)

    functionF = Add(functionF, Symbol('C'))
    h3s4 = "Simplifying and using F(x,y) as a constant C, we get: " + "\\\\ \\\\"
    subSteps.append(h3s4)

    functionF = simplify(functionF)
    
    eq4s4 = "$" + latex(Integer(0)) +  " = " + latex(functionF) + "$" + "\\\\ \\\\"
    subSteps.append(eq4s4)


    '''
    ------------------------------------------------------
    # Step 05: Get the explicit solution solving for y
    ------------------------------------------------------
    '''
    solveArray.append([])
    step = solveArray[4]
    step.append("- Get the explicit solution solving for y" + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]

    global finalSolve
    finalSolve = []
    
    def final_solve_timeout(expression, symbol):
      global finalSolve
      finalSolve = solve(expression, symbol)


    try:
      process = PropagatingThread(target = final_solve_timeout, args=(functionF, Symbol(functionName)))
      process.start()
      process.join(timeout=5)

      solveY = solve(functionF, Symbol('y'))
      for singleSolve in solveY:
        eq1s5 = Eq(y(x), singleSolve)
        subSteps.append("$" + latex(eq1s5) + "$" + "\\\\ \\\\") 
      step.append(subSteps)
      solveArray.append(step)
    except:
      step = []
      step.append("- Can not get the explicit solution solving for y" + "\\\\ \\\\")
      subSteps = []
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
