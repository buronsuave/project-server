from sympy import * 
from sympy.abc import x 
from sympy.parsing import parse_expr 
from anomalies.completeness_anomaly import CompletenessAnomaly

def solveNLinear(odeString, functionName, user_type):
  
  print("In the N linear solver")

  try:
    '''
      ------------------------------------------------------
      # Init solve
      ------------------------------------------------------
    '''
    # init solve array
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
    equation = equation.subs(y(x), Symbol('y'))
    equationsolve = parse_expr("E**(r*x)")

    '''
    ------------------------------------------------------
    # Step 01: Subtitute the potential solution 
    ------------------------------------------------------
    '''
    print("Step 1")
    solveArray.append([])
    step = solveArray[0]
    step.append("- Arrange the equation and subsitute the potential solution"+ "\\\\ \\\\")
    step.append([])
    subSteps = step[1]

    h0 = "Initial equation is given by: " + "\\\\ \\\\"
    subSteps.append(h0)

    eq0 = "$" + latex(Eq(odeLeftSym, odeRightSym)) + "$" + "\\\\ \\\\"
    subSteps.append(eq0)   
    
    eq1 = "$" + latex(equation) + "$" + "\\\\ \\\\"
    subSteps.append(eq1)

    h1 = "Lets propose the solution: " + "\\\\ \\\\"
    subSteps.append(h1) 

    eq2 = "$" + latex(Eq(Symbol('y(x)'),equationsolve)) + "$" + "\\\\ \\\\"
    subSteps.append(eq2)
    
    equation = equation.subs(Symbol('y'), equationsolve)
    h2 = "Substituting" + "\\\\ \\\\"
    subSteps.append(h2) 

    eq3 = "$" + latex(equation) + "$" + "\\\\ \\\\"
    subSteps.append(eq3)  
    
    '''
    ------------------------------------------------------
    # Step 02: Evaluate the solution on the SDE
    ------------------------------------------------------
    '''
    print("Step 2")
    solveArray.append([])
    step = solveArray[1]
    step.append("- Apply derivatives and simplify" + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]

    #Declare temp functions
    functionP = Integer(0)
    functionQ = Integer(0)
    functionF = Integer(0)
    leftPartial = Integer(0)
    functionG = Integer(0)
    functionS = Integer(0)
    functionT = Integer(0)
    #Mul(pow(x, -2), equationsolve.subs(Symbol('r'), 1))
    for term in equation.args[0].args:
      if 'Derivative' in str(term):
        if type(term) is Mul:
          for subterm in term.args:
            if 'Derivative' in str(subterm):
              try:
                expression = term.subs(Derivative(equationsolve, x, subterm.args[1].args[1]), diff( equationsolve, x, subterm.args[1].args[1]))
                functionF = Add( functionF , expression)
              except:
                expression = term.subs(Derivative(equationsolve, x), diff( equationsolve, x))
                functionF = Add( functionF , expression)
        else:
          expression = term.subs(Derivative(equationsolve, x, term.args[1].args[1]), diff( equationsolve, x , term.args[1].args[1]))
          functionF = Add( functionF , expression)
      else:
        if 'r' in str(term):
          functionF = Add( functionF , term)
        else:
          functionT = Mul(term, -1)

    #Derivatives are subsituted
    equation = Eq(functionF, 0)

    h1 = "Applying derivatives: " + "\\\\ \\\\"
    subSteps.append(h1)

    eq1 = "$" + latex(equation) + "$" + "\\\\ \\\\"
    subSteps.append(eq1)

    #Factorize
    equation = factor(equation)
    h2 = "Simplifying: " + "\\\\ \\\\"
    subSteps.append(h2)

    eq2 = "$" + latex(equation) + "$" + "\\\\ \\\\"
    subSteps.append(eq2)     

    '''
    ------------------------------------------------------
    # Step 03: Obtaining the roots 
    ------------------------------------------------------
    '''
    print("Step 3")
    solveArray.append([])
    step = solveArray[2]
    step.append("- Obtain roots" + "\\\\ \\\\")
    step.append([])
    subSteps = step[1]

    #Since e¨x cannot be 0, then function Q must be zero 
    functionQ = expand(Mul(functionF, pow(equationsolve, -1)))
    equationaux = Eq(functionQ, 0)

    #Simplify if posible

    h0 =  "Since the exponential function can not be zero, then: " + "\\\\ \\\\"
    subSteps.append(h0)
    
    equationaux = factor(equationaux)
    eq1 = "$" + latex(equationaux) + "$" + "\\\\ \\\\"
    subSteps.append(eq1)

    #Solve equation for r in Q
    solutions = []
    solutionlist = roots(functionQ)
    for solution in solutionlist:
      for i in range(0, solutionlist[solution]):
        solutions.append(solution)
    solutionlist = solutions 

    if not solutionlist: 
      solutionlist = Poly(functionQ.as_numer_denom()[0]).nroots(10, 80)
      solutionlist = [ round(number, 3) for number in solutionlist ]

    functionF = Integer(0)
    i = 0
    real = 0
    imag = 0
    lastsolution = None
    coef = Integer(1)
    FuncArray = []
    Func = []
    Rows = []
    for solution in solutionlist:
      if solution.compare(lastsolution):
        coef = Integer(1)
      else:
        coef = Mul(coef, x)

      # if isinstance(solution, float):
          #solution = round(solution, 3)
      real = re(solution)
      imag = im(solution)

      subSteps.append("$" + latex(Eq(Indexed(Symbol('r'), i), solution )) + "$" + "\\\\ \\\\")  

      equationsolvere = equationsolve.subs(Symbol('r'), real)
      equationsolveim = equationsolve.subs(Symbol('r'), (imag *I)).rewrite(cos)
    
      functionP = Mul(equationsolve.subs(Symbol('r'), solution), coef)
      Rows.append(functionP)

      functionP = Mul(functionP, Indexed(Symbol("C", real = True), i))
      functionQ = Mul(equationsolvere, equationsolveim ,  Indexed(Symbol("C", real = True), i), coef)
      functionG = Add(functionG, functionQ)
      functionF = Add(functionF, functionP)

      imagpart = Mul(im(equationsolveim.subs(x, re(x))).subs(re(x), x),Indexed(Symbol("K"), i), equationsolvere, coef)
      realpart = Mul(re(equationsolveim.subs(x, re(x))).subs(re(x), x),Indexed(Symbol("C", real = True), i), equationsolvere, coef)
    
      if imag == 0:
        equationsolveim = realpart
      else: 
        equationsolveim = Add(realpart, imagpart)

      functionQ = equationsolveim
      functionS = Add(functionS, functionQ)
      i+=1
      lastsolution = solution
    

    '''
    ------------------------------------------------------
    # Step 04: Substitute the roots and Final homogeneous solution
    ------------------------------------------------------
    '''
    print("Step 4")
    solveArray.append([])
    step = solveArray[3]
    step.append("- Substitute the roots and Final homogeneous solution"+ "\\\\ \\\\")
    step.append([])
    subSteps = step[1]

    #Get Matrices
    expression = 0 
    Rowsaux2 = Rows.copy()
    for i in range(0, len(Rowsaux2)+1):
      Func = []
      Rows = Rowsaux2.copy() 
      if i == len(Rowsaux2):
        Func.append(Rows)
        for j in range(0, len(Rowsaux2) - 1):
          Rowsaux = []
          for funcs in Rows:  
            Rowsaux.append(diff(funcs, x))
          Rows = Rowsaux    
          Func.append(Rowsaux)  
      else:
        Rows[i] = 0
        Func.append(Rows)
        for j in range(0, len(Rowsaux2) - 1):
          Rowsaux = []
          for funcs in Rows:
            expression = Add(expression, diff(funcs, x))
            Rowsaux.append(diff(funcs, x))
          Rows = Rowsaux.copy()
          if j ==  len(Rowsaux2) - 2:
            Rowsaux[i] = 1
          else:
            Rowsaux[i] = 0
          Func.append(Rowsaux)
      FuncArray.append(Func)
    
    h0 = "Substituting the obtained solutions in the proposed solution and adding up: " + "\\\\ \\\\"
    subSteps.append(h0)



    equation = Eq(Symbol('y(x)'), functionF)
    subSteps.append("$" + latex(equation) + "$" + "\\\\ \\\\")
    
    equationHomAux = Eq(Symbol('y(x)'), factor(expand(functionG)))
    if (str(equation) != str(equationHomAux)):
      subSteps.append("$" + latex(equationHomAux) + "$" + "\\\\ \\\\")
    equation = equationHomAux

    equationHomAux = Eq(Symbol('y(x)'), expand(functionS))
    if (str(equation) != str(equationHomAux)):
      subSteps.append("$" + latex(equationHomAux) + "$" + "\\\\ \\\\")
    equation = equationHomAux

    if functionT != 0:
      #Create and DIsplay Matrices
      '''
      ------------------------------------------------------
      # Step 05: Equate Both Sides
      ------------------------------------------------------
      '''
      print("Step 5")
      solveArray.append([])
      step = solveArray[4]
      step.append("- Obtain non constant coeficcients of complementary function"+ "\\\\ \\\\")
      step.append([])
      subSteps = step[1]

      h0 = "Searching up for the system of matrices: " + "\\\\ \\\\"
      subSteps.append(h0)

      Matrices = []
      for matrix in FuncArray:
        Matrices.append(Matrix(matrix))
        subSteps.append("$" + latex(Matrix(matrix)) + "$" +  "\\\\ \\\\ \\\\")

    #Calculate the determinant of each Matrix
      solveArray.append([])
      step = solveArray[5]
      step.append("- Obtaining determinants of each matrix:"+ "\\\\ \\\\")
      step.append([])
      subSteps = step[1]

      h0 = "Obtaining determinants of each function: " + "\\\\ \\\\"
      subSteps.append(h0)
      Dets = []
      j=0
      for matrix in Matrices:
        deti = matrix.det()
        subSteps.append("$" + latex(Eq(Indexed(Symbol('P'), j), deti)) + "$" +  "\\\\ \\\\ \\\\")
        Dets.append(deti)
        j+=1

      '''
      ------------------------------------------------------
      # Step 06: Get Explicit Solve
      ------------------------------------------------------
      '''
      print("Step 6")
      solveArray.append([])
      step = solveArray[5]
      step.append("- Obtaining Integral factors"+ "\\\\ \\\\")
      step.append([])
      subSteps = step[1]

      #Calculate Integral Factors
      Factors = [] 
      subSteps.append("Taking into account that: " + "\\\\ \\\\")
      subSteps.append("$" + latex(Eq(Function('T')(x), functionT)) + "$" + "\\\\ \\\\")
      for i in range(0, len(Dets)-1):

        subSteps.append("$" + latex(Eq(Indexed(Symbol('U'), i), Integral(nsimplify(Mul(Function('T')(x), Indexed(Symbol('P'), i), pow(Indexed(Symbol('P'), len(Dets)-1),-1))), x))) + "$" + "\\\\ \\\\")
        subSteps.append("$" + latex(Eq(Indexed(Symbol('U'), i), Integral(nsimplify(Mul(Function('T')(x), Dets[i], pow(Dets[len(Dets)-1], -1))), x))) + "$" + "\\\\ \\\\")
        subSteps.append("$" + latex(Eq(Indexed(Symbol('U'), i), Integral(nsimplify(Mul(functionT, Dets[i], pow(Dets[len(Dets)-1], -1))), x))) + "$" + "\\\\ \\\\")

        auxExpr = expand(Mul(functionT.simplify(), Dets[i], pow(Dets[len(Dets)-1], -1)))
        auxExpr = auxExpr.rewrite(cos)
        auxExpr = auxExpr.subs(x, re(x))
        auxExpr = im(auxExpr)
        auxExpr = simplify(Mul(auxExpr, I))

        # if functionT.is_real:
        expandAux = expand(Mul(functionT.simplify(), Dets[i], pow(Dets[len(Dets)-1], -1)))
        u = integrate(expandAux, x)
        # else:
        #   u = integrate(expand(Mul(simplify(functionT.rewrite(cos)), Dets[i], pow(Dets[len(Dets)-1], -1)), x))

        Factors.append(u)
        subSteps.append("$" + latex(Eq(Indexed(Symbol('U'), i),u)) + "$" + "\\\\ \\\\")

    
      '''
      ------------------------------------------------------
      # Step 07: Get Explicit Solve
      ------------------------------------------------------
      '''
      print("Step 7")
      solveArray.append([])
      step = solveArray[6]
      step.append("- Obtain final complement by substituting the integral factors in the constants of the homogeneous solution" + "\\\\ \\\\")
      step.append([])
      subSteps = step[1]
      
      functionT = functionF
      for i in range(0, len(Dets)-1):
        functionT  = functionT.subs(Indexed(Symbol('C', real = True), i), Factors[i])    

      h0 = "Substitution: " + "\\\\ \\\\"
      subSteps.append(h0)
      subSteps.append("$" + latex(expand(Eq(Symbol('Y(x)'), functionT))) + "$" + "\\\\ \\\\" )
      subSteps.append("$" + latex(expand(Eq(Symbol('Y(x)'), simplify(functionT)))) + "$" + "\\\\ \\\\" )

      if not (functionT.is_real): 
        h1 = "Rewriting complex functions in terms of cos and sin: " + "\\\\ \\\\"
        subSteps.append(h1)
        subSteps.append("$" + latex(Eq(Symbol('Y(x)'), simplify(functionT.rewrite(cos)).rewrite(cos))) + "$" + "\\\\ \\\\")

      
      '''
      ------------------------------------------------------
      # Step 08: Subtitute integral factors
      ------------------------------------------------------
      '''
      print("Step 8")
      solveArray.append([])
      step = solveArray[7]
      step.append("- Obtain final complement by substituting the integral factors in the constants of the homogeneous solution" + "\\\\ \\\\")
      step.append([])
      subSteps = step[1]
  
      
      subSteps.append("Finally the general solution is the adition of the homogeneous solution and the particular solution:"+ "\\\\ \\\\")
      
      functionF = Add(simplify(functionF), functionT)
      h0 = "Adding both the final solution is: " + "\\\\ \\\\"
      subSteps.append(h0)

      eq1 = "$" + latex(simplify(expand(Eq(Symbol('y(x)'), functionF)))) + "$" + "\\\\ \\\\"
      subSteps.append(eq1)

      print("Step 8.1")
      if not functionF.is_real:  
        functionF = trigsimp(logcombine(simplify(functionF.rewrite(cos)), force=True))

        h1 = "Expressing complex terms as sin and cos: " + "\\\\ \\\\"
        eq2 = "$" + latex(expand(Eq(Symbol('y(x)'), functionF))) + "$" + "\\\\ \\\\"
    
        functionFAux = expand(functionF)
        
        realpart = re(functionFAux.subs(x, re(x))).subs(re(x), x)
        imagpart = im(functionFAux.subs(x, re(x))).subs(re(x), x)

        print("Step 8.2")
        for i in range(0, len(Dets)-1):
          imagpart = imagpart.subs(Indexed(Symbol("C", real=True), i),Indexed(Symbol("K", real=True), i))

        print("Step 8.3")
        functionFAux = imagpart + realpart
        eq3 = "$" + latex(Eq(Symbol('y(x)'), simplify(simplify(functionFAux)))) + "$" + "\\\\ \\\\"
        
        imagpart = abs(imagpart)
        if (len(imagpart.args) > 1):
          print("Step 8.4")
          subSteps.append(h1)
          subSteps.append(eq2)
          subSteps.append(eq3)  

    def display_step(step):
      stepStr = ""
      for subStep in step:
        stepStr += str(subStep)
      return stepStr

    def display_solve(solveArray):
      solveStr = ""
      for stepAux in solveArray:
        if (len(stepAux) > 0):
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
  
  except Exception as e:
    print(e)