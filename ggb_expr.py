#!/usr/bin/env python
# coding: utf-8

# List of token names.   This is always required
tokens = (
    'NUMBER',
    'BSYMBOL0',
    'BSYMBOL1',
    'BSYMBOL2',
    'BSYMBOL3',
    'DEGREE',
    'MINUS',
    'LPAREN',
    'RPAREN',
    'LBRACKET',
    'RBRACKET',
    'COMMA',
    'WORD',
)

t_MINUS    = r'-'
t_BSYMBOL0 = r'[≟⊥∈∥]'
t_BSYMBOL1 = r'\+'
t_BSYMBOL2 = r'[*/]'
t_BSYMBOL3 = r'\^'
t_DEGREE   = r'°'

symbol_to_name = {
    '≟': "Equality",
    '⊥': "ArePerpendicular",
    '∈': "ContainedBy",
    '∥': "AreParallel",
    '+': "Sum",
    '-': "Minus",
    '*': "Product",
    '/': "Ratio",
    '^': "Power",
}

t_LPAREN   = r'\('
t_RPAREN   = r'\)'
t_LBRACKET = r'\['
t_RBRACKET = r'\]'
t_COMMA    = r','

def t_NUMBER(t):
    r'[\d]+(\.[\d]*)?'
    if '.' in t.value: t.value = float(t.value)
    else: t.value = int(t.value)    
    return t

def t_WORD(t):
    r"[A-Za-zα-ω][A-Za-zα-ω_'0-9]*"
    return t

t_ignore  = ' \t'

def t_error(t):
    print("Illegal character '%s'" % t.value[0])
    t.lexer.skip(1)

# Build the lexer
import ply.lex as lex
lex.lex()

import ply.yacc as yacc

class Expression:
    def __init__(self, name, subexpr_list):
        self.name = name
        self.subexpr_list = subexpr_list

    def indent_print(self, indent = 0):
        print("  "*indent + self.name)
        for subexpr in self.subexpr_list:
            subexpr.indent_print(indent+1)

class ConstExpr(Expression):
    def __init__(self, value, in_degrees = False):
        self.in_degrees = in_degrees
        self.value = value
        self.name = str(value)
        if in_degrees: self.name += '°'
        self.subexpr_list = []

class VariableExpr(Expression):
    def __init__(self, name):
        self.name = name
        self.subexpr_list = []

precedence = (
    ('nonassoc', 'BSYMBOL0'),       # equality, relations
    ('left', 'BSYMBOL1', 'MINUS'),  # addition, subtraction
    ('left', 'BSYMBOL2'),           # multiplication, division
    ('right', 'BSYMBOL3'),          # powers
)

def p_expression_b0(p):
    '''expression : expression BSYMBOL0 expression
                  | expression MINUS expression
                  | expression BSYMBOL1 expression
                  | expression BSYMBOL2 expression
                  | expression BSYMBOL3 expression
    '''
    p[0] = Expression(symbol_to_name[p[2]], [p[1], p[3]])

def p_expr_uminus(p):
    'expression : MINUS expression'
    p[0] = Expression(symbol_to_name[p[1]], [p[2]])

def p_nakedlist_expr(p):
    'nakedlist : expression'
    p[0] = [p[1]]
def p_nakedlist_comma(p):
    'nakedlist : nakedlist COMMA expression'
    p[0] = p[1] + [p[3]]

def p_expression_function(p):
    'expression : WORD LBRACKET nakedlist RBRACKET'
    p[0] = Expression(p[1], p[3])
def p_expression_paren(p):
    'expression : LPAREN expression RPAREN'
    p[0] = p[2]

def p_expression_const(p):
    'expression : NUMBER'
    p[0] = ConstExpr(p[1])
def p_expression_const_angle(p):
    'expression : NUMBER DEGREE'
    p[0] = ConstExpr(p[1], in_degrees = True)
def p_expression_one_degree(p):
    'expression : DEGREE'
    p[0] = ConstExpr(1, in_degrees = True)
def p_expression_var(p):
    'expression : WORD'
    p[0] = VariableExpr(p[1])

def p_error(p):
    print("Syntax error in input!")

def make_parser():
    parser = yacc.yacc(tabmodule='ggb_parsetab')
    return parser
 
if __name__ == '__main__':
    parser = make_parser()
    s = "AreEqual[Polygon[Translate[A, Vector[Point[yAxis], (3, 1)]], Translate[B, Vector[Point[yAxis], (3, 1)]], Translate[C, Vector[Point[yAxis], (3, 1)]], Translate[D, Vector[Point[yAxis], (3, 1)]], Translate[E, Vector[Point[yAxis], (3, 1)]], Translate[F, Vector[Point[yAxis], (3, 1)]]], poly1]"
    tree = parser.parse(s, debug = False)
    #tree.indent_print(0)
