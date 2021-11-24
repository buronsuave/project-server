from anomalies.classification_anomaly import ClassificationAnomaly
from sympy import *
from sympy.abc import x,z
from sympy.parsing import parse_expr

#FIXME: Check with some examples our methods for classification (sume bugs)

def checkSeparable(odeString, functionName):
    # Testing completeness anomaly
    # return False

    odeLeftString = odeString.split("=")[0]
    odeRightString = odeString.split("=")[1]
    odeLeftSym = parse_expr(odeLeftString)
    odeRightSym = parse_expr(odeRightString)
    y = Function(functionName)
    equation = Eq(odeLeftSym - odeRightSym, 0)
    left = equation.args[0]
    try:
        express = solve(Eq(left, Integer(0)), Derivative(y(x), x))
    except:
        return False

    if (len(express) == 0):
        return False

    aux = simplify(express[0]) 
    aux = aux.expand(force=True)
    aux = factor(aux)

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
                    functionG = Mul(functionG, term)
                else:
                    functionF = Mul(functionF, term)

    if functionName in str(functionF):
        return False

    functionG = functionG.subs(y(x), Symbol(functionName))
    # functionG = functionG.subs(Symbol('x'), Symbol(functionName))

    if 'x' in str(functionG):
        return False

    # functionG = functionG.subs(Symbol(functionName), Symbol('x'))
    functionG = functionG.subs(Symbol(functionName), y(x))

    if (Mul(functionF, functionG)) == aux:
        return True

    return False

def checkLinear(odeString):
    # Testing completeness anomaly
    # return False

    odeLeftString = odeString.split("=")[0]
    odeRightString = odeString.split("=")[1]

    odeLeftSym = parse_expr(odeLeftString)
    odeRightSym = parse_expr(odeRightString)

    y = Function('y')
    equation = Eq(odeLeftSym - odeRightSym, 0)
    left = equation.args[0]
    try:
        express = solve(Eq(left, Integer(0)), Derivative(y(x), x))
    except:
        return False
    
    if (len(express) == 0):
        return False
    
    aux = expand(express[0])

    left = Derivative(y(x), x)

    functionF = parse_expr("0")
    functionG = parse_expr("0")

    for term in aux.args:
        if "y" in str(term):
            functionF = Add(functionF, Mul(term, Pow(y(x), Integer(-1))))
        else:
            functionG = Add(functionG, term)

    functionF = Mul(functionF, Integer(-1))
    functionF = simplify(functionF)
    functionG = simplify(functionG)

    if 'y' in str(functionF):
        return False
    
    if 'y' in str(functionG):
        return False

    right = Add(functionG, Mul(Integer(-1), functionF, y(x)))
    right = right.expand(force=True)

    if Add(aux, Mul(Integer(-1), right)).simplify() == Integer(0):
        return True

    return False

def checkReducibleLinear(odeString):
    # Testing completeness anomaly
    # return False

    odeLeftString = odeString.split("=")[0]
    odeRightString = odeString.split("=")[1]

    odeLeftSym = parse_expr(odeLeftString)
    odeRightSym = parse_expr(odeRightString)

    y = Function('y')
    equation = Eq(odeLeftSym - odeRightSym, 0)
    
    left = equation.args[0]
    exp = solve(left, Derivative(y(x), x))
    aux = expand(exp[0])

    functionF = parse_expr("0")
    functionG = parse_expr("0")
    
    aux = Mul(aux, Pow(y(x), Integer(-1)))
    aux = simplify(aux)

    for term in aux.args:
        if 'y' in str(term):
            for subTerm in term.args:
                if 'y' in str(subTerm):
                    try:
                        n = Add(subTerm.args[1], Integer(1))
                        subG = Mul(term, Pow(subTerm, Integer(-1)))
                        functionG = Add(functionG, subG)
                    except:   
                        n = None             
        else:
            functionF = Add(functionF, term)
    
    print(n)
    print(functionF)
    print(functionG)

    return False

    if functionG != 0 and functionF != 0:
        return True
    else:
        return False

def checkExact(odeString):
    # Testing completeness anomaly
    # return False

    odeLeftString = odeString.split("=")[0]
    odeRightString = odeString.split("=")[1]

    odeLeftSym = parse_expr(odeLeftString)
    odeRightSym = parse_expr(odeRightString)

    y = Function('y')
    equation = Eq(odeLeftSym - odeRightSym, 0)
    equation = equation.subs(y(x), Symbol('y'))

    functionP = Integer(0)
    functionQ = Integer(0)

    for term in equation.args[0].args:
        if 'Derivative' in str(term):
            if len(term.args[1]) > 1:
                return False
                
            functionQ = Add(functionQ, Mul(term, Pow(Derivative(Symbol('y'), x), Integer(-1))))
        else:
            functionP = Add(functionP, term)

    if (functionP == 0 or functionQ == 0):
        return False

    partialP = diff(functionP, Symbol('y'))
    partialQ = diff(functionQ, Symbol('x'))

    if Add(partialP, Mul(Integer(-1), partialQ)).simplify() == Integer(0):
        return True                                                                                                                           

    return False

def checkHomogeneous(odeString):
    # Testing completeness anomaly
    # return False

    odeLeftString = odeString.split("=")[0]
    odeRightString = odeString.split("=")[1]

    odeLeftSym = parse_expr(odeLeftString)
    odeRightSym = parse_expr(odeRightString)

    y = Function('y')
    equation = Eq(odeLeftSym - odeRightSym, 0)
    
    left = equation.args[0]
    exp = solve(left, Derivative(y(x), x))
    aux = expand(exp[0])
    
    left = Derivative(y(x), x)
    functionF = aux
  
    u = Function('u')

    functionF = functionF.subs(y(x), Mul(u(x), x))
    left = Add(Mul(Derivative(u(x), x), x), u(x))
    
    left = Add(left, Mul(functionF, Integer(-1)))
    left = expand(left)
    separableODE = Eq(left, Integer(0))
    try:
        print(str(separableODE.args[0]) + "= 0")
        return checkSeparable(str(separableODE.args[0]) + "= 0", "u")
    except: 
        return False

def checkSuperiorOrder(odeString):
    # Testing completeness anomaly
    # return False

    odeLeftString = odeString.split("=")[0]
    odeRightString = odeString.split("=")[1]

    odeLeftSym = parse_expr(odeLeftString)
    odeRightSym = parse_expr(odeRightString)

    y = Function('y')
    equation = Eq(odeLeftSym - odeRightSym, 0)
    equation = equation.subs(y(x), Symbol('y'))

    for term in equation.args[0].args:
        if 'Derivative' in str(term):
            print(term)
            if type(term) is Mul:
                for subterm in term.args:
                    if not ('Derivative' in str(subterm)):
                        if 'x' in str(subterm) or 'y' in str(subterm):
                            return False            
                    else:
                        continue
            else:
                continue  
        else:
            if 'x' in str(term) and 'y' in str(term):
                return False

    return True

def classify(odeString):
    
    # Check if the equation is separable
    try:
        if checkSeparable(odeString, "y"):
            return "separable"
    except:
        print("Non Separable")

    # Check if the equation is 1st order linear
    try:
        if checkLinear(odeString):
            return "linear"
    except Exception as e:
        print(e)
        print("Non Linear")
    
    # Check if the equation is reducible to linear
    try:
        if checkReducibleLinear(odeString):
            return "reducible"
    except:
        print("Non reducible")

    # Check if the equation is 1st order exact
    try:
        if checkExact(odeString):
            return "exact"
    except:
        print("Non Exact")
    
    # Check if the equation is 1st order homogeneous
    try:
        if checkHomogeneous(odeString):
            return "homogeneous"
    except:
        print("Non Homogeneous")
    
    # Check if the equation is linear of superior order with constant coefficients
    try:
        if checkSuperiorOrder(odeString):
            return "superior"
    except Exception as excep:
        print(excep.args[0])
        print("Non Superior Order")

    # Use dsolve to classify the equation
    try:
        odeSym = parse_expr(odeString.split('=')[0])
        odeEq = Eq(odeSym, 0)

        # Gets all the possible ways to see the equation
        odeClass = classify_ode(odeSym, Function('y')(x))
        
        # Iterate in odeClass list and check if match
        for odeType in odeClass:
            print(odeType)
            if odeType is 'separable':
                return 'separable'
            elif odeType is '1st_exact':
                return 'exact'
            elif odeType is '1st_linear':
                return 'linear'
            elif odeType is 'Bernoulli_Integral':
                return 'reducible'
    
    except Exception as e:
        print(e.args[0])
    
    raise ClassificationAnomaly()