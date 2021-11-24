from sympy import *
from sympy.abc import x
from sympy.parsing.latex import parse_latex

import re

def parseSympy(inputString):
    inputString = inputString.replace("D", "Derivative")
    inputString = inputString.replace("^", "**")
    inputString = inputString.replace("y", "y(x)")
    inputString = inputString.replace(
        "Derivative(y(x),x)", "(Derivative(y(x),x))")

    print(inputString)

    therms = inputString.split("=")
    if (len(therms) != 2):
        raise Exception("Invalid equation. Must be exactly one equal sign (=)")

    leftString = inputString.split("=")[0]
    rightString = inputString.split("=")[1]
    leftSymbol = parse_expr(leftString)
    rightSymbol = parse_expr(rightString)

    leftSide = Add(leftSymbol, Mul(rightSymbol, Integer(-1)))

    return leftSide

def regexCorrect(variable):
    regexps = [r'(\\sin[^)]*\)\})', r'(\\cos[^)]*\)\})', r'(\\tan[^)]*\)\})', 
        r'(\\csc[^)]*\)\})', r'(\\sec[^)]*\)\})', r'(\\cot[^)]*\)\})']

    for regexp in regexps:
        foo = re.compile(regexp)
        searchAll = foo.findall(variable)
        for oldValue in searchAll:
            variable = variable.replace(
                oldValue, "\\left(" + oldValue + "\\right)")
    return variable

def newRegexCorrect(variable):
    regexp = r'y\^\{[^\}]*\}\{\\left\(x \\right\)\}'
    foo = re.compile(regexp)
    searchAll2 = foo.findall(variable)
    for oldValue in searchAll2:
        startPowIndex = oldValue.index("{")
        endPowIndex = oldValue.index("}")
        powString = oldValue[(startPowIndex+1): endPowIndex]
        variable = variable.replace(
            oldValue, "M^{" + powString + "}"
        )
    return variable

def fixMFunction(expression):
    for item in expression.args:
        if type(item) is Function('M'):
            arguments = item.args[0]
            expression = expression.subs(Function('M')(
                arguments), Mul(Symbol('M'), arguments))
        elif "Function('M')" in str(srepr(item)):
            newItem = fixMFunction(item)
            expression = expression.subs(item, newItem)
    return expression

def fixLogarithm(expression):
    for item in expression.args:
        if type(item) is log:
            arguments = item.args[0]
            expression = expression.subs(item, ln(arguments))
        elif "log" in str(srepr(item)):
            newItem = fixLogarithm(item)
            expression = expression.subs(item, newItem)
    return expression

def fixXFunction(expression):
    for item in expression.args:
        if type(item) is Function('x'):
            arguments = item.args[0]
            expression = expression.subs(Function('x')(
                arguments), Mul(Symbol('x'), arguments))
        elif "Function('x')" in str(srepr(item)):
            newItem = fixXFunction(item)
            expression = expression.subs(item, newItem)
    return expression

def parseLatex(inputString):
    inputString = inputString.replace("\\\\", "\\")
    inputString = inputString.replace(
        "y{\\left(x \\right)}", "M"
    )
    inputString.replace(
        "log", "ln"
    )
    inputString = regexCorrect(inputString)
    inputString = newRegexCorrect(inputString)

    therms = inputString.split("=")
    if (len(therms) != 2):
        raise Exception("Invalid equation. Must be exactly one equal sign (=)")

    leftString = inputString.split("=")[0]
    leftSymbol = parse_latex(leftString)

    leftSymbol = fixMFunction(leftSymbol)
    leftSymbol = fixLogarithm(leftSymbol)
    leftSymbol = fixXFunction(leftSymbol)

    leftSymbolStr = str(leftSymbol)
    leftSymbolStr = leftSymbolStr.replace("M", "y(x)")
    leftSymbol = parse_expr(leftSymbolStr)

    leftSymbol = leftSymbol.subs(Symbol('e'), E)

    regexp =  r'd\*\*[2-9]\*y\(x\)\/dx\*\*[2-9]'
    foo= re.compile(regexp)
    leftSymbolString = str(leftSymbol)
    Searchall= foo.findall(leftSymbolString)

    if len(Searchall) > 0:
        for oldvalue in Searchall:
            leftSymbolString = leftSymbolString.replace(oldvalue, "Derivative(y(x),x," + oldvalue[-1] + ")")
        leftSymbol = parse_expr(leftSymbolString)

    return leftSymbol
