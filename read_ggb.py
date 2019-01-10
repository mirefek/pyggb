from zipfile import ZipFile
from xml.etree.ElementTree import parse
import os

import ply.yacc as yacc
import ggb_expr
parser = ggb_expr.make_parser()

def get_values(node):
    d = node.attrib
    i = 0
    key = "a{}".format(i)
    result = []
    while key in d:
        result.append(d[key])
        i+=1
        key = "a{}".format(i)

    return result

class Element:
    def __init__(self, label, state, data):
        self.label = label
        self.data = data
        self.state = state # "to_update", "to_keep", or "set"
    def __repr__(self):
        return "'{}':{}".format(self.label, self.data.__repr__())

class Command:
    def __init__(self, operation, inputs, outputs):
        self.operation = operation
        self.inputs = inputs
        self.outputs = outputs
    def __str__(self):
        inputs_str = ' '.join(self.inputs)
        outputs_str = ' '.join([s if s is not None else "_" for s in self.outputs])
        return "{} : {} -> {}".format(
            self.operation, inputs_str, outputs_str
        )

class ConstCommand:
    def __init__(self, datatype, value, label):
        self.datatype = datatype
        self.value = value
        self.label = label
    def __str__(self):
        if self.datatype == int: datatype_str = "int"
        elif self.datatype == AngleSize: datatype_str = "AngleSize"
        elif self.datatype == Measure: datatype_str = "Measure"
        elif self.datatype == Boolean: datatype_str = "Boolean"
        else: raise Exception("Unknown const datatype '{}'".format(self.datatype))
        return "const {} {} -> {}".format(datatype_str, self.value, self.label)

def tweak_command_name(name): # "AreCollinear" -> "are_collinear"
    altered_name = [name[0].lower()]
    for x in name[1:]:
        if x.isupper(): altered_name+=['_', x.lower()]
        else: altered_name.append(x)
    return ''.join(altered_name)

from geo_types import *

type_to_shortcut = {
    int       : 'i',
    Boolean   : 'b',
    Measure   : 'm',
    Point     : 'p',
    Polygon   : 'P',
    Circle    : 'c',
    Arc       : 'C',
    Line      : 'l',
    Ray       : 'r',
    Segment   : 's',
    Angle     : 'a',
    AngleSize : 'A',
    Vector    : 'v',
}

def command_types_name(name, params):
    return "{}_{}".format(name, ''.join(type_to_shortcut[type(x)] for x in params))

import commands as commands_module
from inspect import getmembers, isfunction
command_dict = dict(o for o in getmembers(commands_module) if isfunction(o[1]))

def apply_command(name, params):
    name_with_types = command_types_name(name, params)
    if name_with_types in command_dict: name = name_with_types
    elif name not in command_dict:
        raise Exception("Command '{}' not found".format(name_with_types))

    f = command_dict[name]

    result = f(*params)
    if not isinstance(result, (tuple, list)): result = (result,)
    try:
        for res in result:
            assert(type(res) in type_to_shortcut)
    except:
        raise Exception("Command '{}' returned an unsupported type {}".format(name, type(res)))
    return result

def xml_element(node):
    t = node.attrib["type"]
    label = node.attrib["label"]
    if t in ("conic", "conicpart"):
        matrix = node.find('matrix').attrib
        A = [float(matrix['A'+str(i)]) for i in range(6)]
        assert(A[3] == 0)
        if A[0] == 0:
            assert(A[1] == 0)
            return Element(label,
                "to_keep" if t == "conic" else "to_update",
                Line(A[4:6], -A[2]/2)
            )
        else:
            assert(A[0] == 1 and A[1] == 1)
            return Element(label,
                "to_keep" if t == "conic" else "to_update",
                Circle([-A[4], -A[5]], np.sqrt(A[4]**2 + A[5]**2 - A[2])),
            )
    elif t in ("line", "segment", "ray", "point", "vector"):
        coords = node.find('coords').attrib
        x, y, z = (float(coords[key]) for key in ('x', 'y', 'z'))
        if t == "point":
            assert(z != 0)
            return Element(label, "to_keep", Point([x/z,y/z]))
        elif t == "vector":
            assert(z == 0)
            return Element(label, "to_update", Vector(((0,0),(x,y))))
        return Element(label,
            "to_keep" if t == "line" else "to_update",
            Line([x,y], -z),
        )
    elif t in ("numeric", "angle"):
        value = float(node.find('value').attrib['val'])
        if t == "numeric": return Element(label, "to_update", Measure(value, 0))
        else: return Element(label, "to_update", AngleSize(value))
    elif t in ("polygon", "boolean"): return Element(label, "empty", None)
    elif t in ("function", "locus"): return None
    else: raise Exception("element of unknown type '{}'".format(t))

command_set = set()
conic_files = set()
type_set = set()

def process_ggb(filename):
    print(filename)
    ggb = ZipFile(filename, 'r').open("geogebra.xml")
    construction = parse(ggb).getroot().find('construction')

    element_dict = dict()

    for element in construction.iter('element'):
        el = xml_element(element)
        if el is not None: element_dict[el.label] = el

    command_list = []
    def get_element(label):
        if label not in element_dict:
            raise Exception("label '{}' not declared".format(label))
        el = element_dict[label]
        if el.state == "set": return el.data
        elif type(el.data) == Point:
            el.state = "set"
            command_list.append(Command("point", [], [label]))
        elif type(el.data) == Measure:
            el.state = "set"
            command_list.append(ConstCommand(Measure, el.x, label))
        return el.data

    expr_index = 0
    def get_expr_node(node, label = None):
        if type(node) == ggb_expr.VariableExpr:
            assert(label is None)
            return get_element(node.name), node.name
        else:
            params = [get_expr_node(subexpr) for subexpr in node.subexpr_list]
            if label is None:
                nonlocal expr_index
                label = "expr{}".format(expr_index)
                expr_index += 1
            if type(node) == ggb_expr.ConstExpr:
                value = node.value
                if node.in_degrees:
                    datatype = AngleSize
                    value *= np.pi
                    value /= 180
                else: datatype = type(node.value)

                if datatype == float: datatype = Measure
                command_list.append(ConstCommand(datatype, value, label))
                return datatype(value), label
            else:
                command_name = tweak_command_name(node.name)
                param_data, param_labels = zip(*params)
                command_list.append(Command(command_name, param_labels, [label]))
                (result,) = apply_command(command_name, param_data)
                return result, label

    def get_expr(expr_str, root_label = None):
        root = parser.parse(expr_str)
        return get_expr_node(root, label = root_label)

    for node in construction:
        if node.tag == "command":
            command_name = tweak_command_name(node.attrib["name"])
            if command_name == "locus": continue
            input_expr = get_values(node.find("input"))
            output_labels = get_values(node.find("output"))
            for x in output_labels:
                if x != "": break
            else: continue

            input_data, input_labels = zip(*tuple(get_expr(expr) for expr in input_expr))
            if type(input_data[-1]) == int and \
               command_types_name(command_name, input_data) not in command_dict and \
               command_types_name(command_name, input_data[:-1]) in command_dict and \
               len(output_labels) == 1:
                input_data = input_data[:-1]
                input_labels = input_labels[:-1]
                command_list.pop()

            #print(command_name, input_labels)
            result = apply_command(command_name, input_data)
            output_elements = [element_dict[x] for x in output_labels if x != ""]

            # match result with output elements

            has_empty = False
            for el in output_elements:
                if el.state == "empty": has_empty = True

            if has_empty:
                assert(len(result) == len(output_labels))
                result_labels = [(label if label != "" else None) for label in output_labels]
                for result_el, label in zip(result, output_labels):
                    if label == "": continue
                    ori_el = element_dict[label]
                    assert(ori_el.state != "set")
                    if ori_el.state in ("to_keep", "to_update"):
                        if not result_el.equivalent(ori_el.data):
                            #print(output_labels)
                            #print(result)
                            print(ori_el)
                            print(result_el)
                            raise Exception("Incompatible expected output and output")
                    if ori_el.state in ("empty", "to_update"): ori_el.data = result_el
                    ori_el.state = "set"
            elif command_name == "point" and len(input_data) == 1:
                assert(len(output_elements) == 1)
                assert(len(output_labels) == 1)
                el = output_elements[0]
                assert(el.state == "to_keep")
                assert(type(el.data) == Point)
                assert(input_data[0].contains(el.data.a))
                result_labels = output_labels
                el.state = "set"
            else:
                assert(len(result) >= len(output_elements))
                result_labels = [None for _ in result]
                for el in output_elements:
                    try:
                        i, obj = next(filter(lambda i_obj: el.data.equivalent(i_obj[1]), enumerate(result)))
                    except StopIteration:
                        print(output_elements)
                        print(result)
                        raise Exception("No equivalent result found")
                    assert(result_labels[i] == None)
                    result_labels[i] = el.label
                    assert(el.state != "set")
                    if el.state == "to_update": el.data = obj
                    el.state = "set"

            command_list.append(Command(command_name, input_labels, result_labels))

        if node.tag == "expression":
            label = node.attrib["label"]
            expr = node.attrib["exp"]
            result, _ = get_expr(expr, label)
            ori_el = element_dict[label]
            if ori_el.state in ("to_keep", "to_update"):
                if not result.equivalent(ori_el.data):
                    print(ori_el)
                    print(result)
                    raise Exception("Incompatible expected output and output")
            if ori_el.state in ("empty", "to_update"): ori_el.data = result
            ori_el.state = "set"

    return command_list

def process_dir(datadir):
    for filename in os.listdir(datadir):
        if not filename.endswith(".ggb"): continue
        command_list = process_ggb(os.path.join(datadir, filename))
        f = open(os.path.join(datadir, "decoded", filename[:-4]+".txt"), 'w')
        for command in command_list: f.write("{}\n".format(command))

if __name__ == "__main__":
    process_dir("../tests/true/")
    process_dir("../tests/false/")
