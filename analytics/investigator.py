from sympy import *
import base64
import uuid
import os
x = Symbol('x')

# Researcher Properties and parameters 

def get_roots(expression, solution_list):
  print(expression)
  solutions = []

  # expression = expression.subs(Symbol('C'), 0)

  try:
    solutions = solve(Eq(expression, 0), x)

    if len(solutions) == 0:
      raise Exception()

    for solution in solutions:
      solution_list.append(solution)

  except:
    try:
      solutions = Poly(expression.as_numer_denom()[0]).nroots(10, 80)
      print(solutions)
      solutions = [ round(number, 4) for number in solutions ]
    except:
      try:
        solutions = nsolve(expression, 0, verify = false)
        for solution in solutions:
          solution_list.append(solution)
      # solutionlist = [ round(number, 3) for number in solutionlist ]
      except Exception as e:
        raise e
  # if not solutionlist: 
  #   solutionlist = Poly(expression.as_numer_denom()[0]).nroots(10, 80)
  #   solutionlist = [ round(number, 3) for number in solutionlist ]


def validate(expression, solution_list):
  tolerance = 0.5
  if solution_list is None:
    return False

  for solution in solution_list:
    eval = expression.evalf(subs = {x:solution})
    if re(eval) > tolerance or re(eval) < -1 * tolerance:
      return False
  return True

def max_min(expression, solution_list):
  # expression = parse_expr(expression)

  try:
    dexpression = diff(expression, x)
    solutions = []

    get_roots(dexpression, solutions)
    if validate(dexpression, solutions):
      for solution in solutions:
        solution_list.append(solution)
  except Exception as e:
    raise e

def inflexion_points(expression, solution_list):
  # expression = parse_expr(expression)

  try:
    dexpression = diff(expression, x, 2)
    solutions = []

    get_roots(dexpression, solutions)
    if validate(dexpression, solutions):
      for solution in solutions:
        solution_list.append(solution)
  except Exception as e:
    raise e

def create_plot(expression):
  p0 = plot(expression.subs(Symbol('C'), -10), show = False)
  min = -9
  max = 10
  for i in range(min, max, 1):
    pi = plot(expression.subs(Symbol('C'), i), show = False, ylim = (-300, 300))
    p0.append(pi[0])
  
  try:
    id = str(uuid.uuid1())
    file = f'analytics/images/{id}.png'
    p0.save(file)
    
    base64_string = None
    with open(file, 'rb') as image_file:
      base64_string = str(base64.encodebytes(image_file.read()))

    os.remove(file)
    return base64_string

  except Exception as e:
    print(e)