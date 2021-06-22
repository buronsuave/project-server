from sympy import *
from sympy.abc import x
from sympy.parsing import parse_expr

def checkSeparable(odeString):
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

    aux = simplify(express[0])
    aux = aux.expand(force=True)
    aux = factor(aux)

    functionF = Integer(1)
    functionG = Integer(1)

    if type(aux) is Add:
        if 'y' in str(aux):
            functionG = aux
        else:
            functionF = aux
    else:
        if not ('y' in str(aux)):
            functionF = aux
        else:
            for term in aux.args:
                if 'y' in str(term):
                    functionG = Mul(functionG, term)
                else:
                    functionF = Mul(functionF, term)

    if 'y' in str(functionF):
        return False

    functionG = functionG.subs(y(x), Symbol('y'))
    functionG = functionG.subs(Symbol('x'), Symbol('u'))

    if 'u' in str(functionG):
        return False

    functionG = functionG.subs(Symbol('u'), Symbol('x'))
    functionG = functionG.subs(Symbol('y'), y(x))

    if (Mul(functionF, functionG)) == aux:
        return True

    return False

def checkLinear(odeString):
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

def checkExact(odeString):
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
            functionQ = Add(functionQ, Mul(term, Pow(Derivative(Symbol('y'), x), Integer(-1))))
        else:
            functionP = Add(functionP, term)
    
    partialP = diff(functionP, Symbol('y'))
    partialQ = diff(functionQ, Symbol('x'))
    if Add(partialP, Mul(Integer(-1), partialQ)).simplify() == Integer(0):
        return True                                                                                                                           

    return False

def checkHomogeneous(odeString):
    return true;

def classify(odeString):
    try:
        if checkSeparable(odeString):
            return "separable"
    except:
        print("Non Separable")

    try:
        if checkLinear(odeString):
            return "linear"
    except:
        print("Non Linear")

    try:
        if checkExact(odeString):
            return "exact"
    except:
        print("Non Exact")
    
    try:
        if checkHomogeneous(odeString):
            return "homogeneous"
    except:
        print("Non Homogeneous")

    return "undefined"
