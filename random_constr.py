import os
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

class Element:
    def __init__(self, label, element_dict):
        self.data = None
        self.label = label
        assert(label not in element_dict)
        element_dict[label] = self

    def drawable(self):
        return isinstance(self.data, (Point, Line, Angle, Polygon, Circle, Vector))
    def draw(self, cr, corners):
        if self.drawable(): self.data.draw(cr, corners)
    def important_points(self):
        if self.drawable(): return self.data.important_points()
        else: return []

    def has_value(self):
        return isinstance(self.data, (Measure, Boolean, AngleSize, Angle, Segment, Polygon))
    def value(self):
        if isinstance(self.data, (Measure, AngleSize)): return self.data.x
        elif isinstance(self.data, Boolean): return float(self.data.b)
        elif isinstance(self.data, Angle): return self.data.angle
        elif isinstance(self.data, Segment): return self.data.length
        elif isinstance(self.data, Polygon): return command_module.area_P(self.data).x
        else: return None

class Command:
    def __init__(self, command_name, input_elements, output_elements):
        self.name = command_name
        self.input_elements = input_elements
        self.output_elements = output_elements

    def apply(self):
        input_data = [x.data for x in self.input_elements]
        name = command_types_name(self.name, input_data)
        if name not in command_dict: name = self.name
        f = command_dict[name]
        output_data = f(*input_data)
        if not isinstance(output_data, (tuple, list)): output_data = (output_data,)
        assert(len(output_data) == len(self.output_elements))
        for x,o in zip(output_data, self.output_elements):
            if o is not None: o.data = x

    def __str__(self):
        inputs_str = ' '.join([x.label for x in self.inputs])
        outputs_str = ' '.join([x.label if x is not None else "_" for x in self.outputs])
        return "{} : {} -> {}".format(
            self.name, inputs_str, outputs_str
        )

const_type_to_str = {
    int : "int",
    AngleSize : "AngleSize",
    Measure : "Measure",
}
str_to_const_type = dict((s,t) for (t,s) in const_type_to_str.items())

class ConstCommand:
    def __init__(self, datatype, value, element):
        self.datatype = datatype
        self.value = value
        self.element = element

    def apply(self):
        self.element.data = self.datatype(self.value)

    def __str__(self):
        datatype_str = const_type_to_str[self.datatype]
        return "const {} {} -> {}".format(datatype_str, self.value, self.label)

def parse_command(line, element_dict):
    tokens = line.split()
    if tokens[0] == "const":
        assert(len(tokens) == 5)
        datatype = str_to_const_type[tokens[1]]
        value = float(tokens[2])
        assert(tokens[3] == "->")
        label = tokens[4]
        element = Element(label, element_dict)
        command = ConstCommand(datatype, value, element)
        element.command = command
        return command
    else:
        command_name = tokens[0]
        assert(tokens[1] == ":")
        labels = [None if token == "_" else token for token in tokens[2:]]
        arrow_index = labels.index("->")
        input_labels = labels[:arrow_index]
        output_labels = labels[arrow_index+1:]
        input_elements = [element_dict[label] for label in input_labels]
        def element_or_none(label):
            if label is None: return None
            else: return Element(label, element_dict)
        output_elements = list(map(element_or_none, output_labels))
        command = Command(command_name, input_elements, output_elements)
        for el in output_elements:
            if el is not None: el.command = command
        return command

class Construction:
    def __init__(self, display_size = (100,100), min_border = 0.1, max_border = 0.25):
        self.corners = np.array(((0,0), display_size))
        self.min_border = min_border
        self.max_border = max_border
        self.nc_commands = []
        self.to_prove = None
        self.element_dict = dict()
        self.elements = []

    def render(self, cr, elements = None): # default: render all elements
        if elements is None: elements = self.elements

        for el in elements:
            el.draw(cr, self.corners)

    def render_to_numpy(self, elements = None):
        surface = cairo.ImageSurface(cairo.FORMAT_A8, self.width, self.height)
        cr = cairo.Context(surface)
        self.render(cr, elements)

        data = surface.get_data()
        data = np.array(data, dtype = float)/255
        data = data.reshape([self.height, surface.get_stride()])
        data = data[:,:self.width]
        return data

    def load(self, filename):
        self.nc_commands = []
        self.to_prove = None
        self.element_dict = dict()
        with open(filename, 'r') as f:
            for line in f:
                command = parse_command(line, self.element_dict)
                if isinstance(command, ConstCommand): command.apply()
                elif isinstance(command, Command):
                    if command.name == "prove":
                        [inp] = command.input_elements
                        [out] = command.output_elements
                        if out is not None: del self.element_dict[out.label]
                        assert(self.to_prove is None)
                        self.to_prove = inp

                    else: self.nc_commands.append(command)

        assert(self.to_prove is not None)
        self.elements = list(self.element_dict.values())

    def run_commands(self):
        for command in self.nc_commands: command.apply()

    def generate(self, require_theorem = True, max_attempts = 100): # max_attempts = 0 -> inf
        while True:
            try:
                self.run_commands()
            except:
                max_attempts -= 1
                if max_attempts == 0: raise
                continue
            if require_theorem and not self.to_prove.data.b: continue
            break

        self.fit_to_window()

    def fit_to_window(self):
        important_points = []
        for el in self.elements: important_points += el.important_points()
        if len(important_points) == 0: return
        src_corners = np.stack([
            np.min(important_points, axis = 0),
            np.max(important_points, axis = 0),
        ])
        src_size = np.maximum(0.01, src_corners[1] - src_corners[0])

        dest_size = self.corners[1] - self.corners[0]
        dest_corners_shift = np.random.random(size = [2,2])
        dest_corners_shift *= self.max_border - self.min_border
        dest_corners_shift += self.min_border
        dest_corners_shift *= np.array((1,-1)).reshape((2,1)) * dest_size
        dest_corners = self.corners + dest_corners_shift
        dest_size = dest_corners[1] - dest_corners[0]

        scale = np.min(dest_size / src_size)
        src_corners *= scale
        shift = np.average(dest_corners, axis = 0) - np.average(src_corners, axis = 0)
        for el in self.elements:
            if isinstance(el.data, (int, float)): continue
            el.data.scale(scale)
            el.data.translate(shift)

        important_points = []
        for el in self.elements: important_points += el.important_points()
        corners = np.stack([
            np.min(important_points, axis = 0),
            np.max(important_points, axis = 0),
        ])

    def test(self, num_tests = 100):
        constr_fail, check_fail, success = 0, 0, 0
        for _ in range(num_tests):
            try:
                self.run_commands()
                if self.to_prove.data.b: success += 1
                else: check_fail += 1
            except:
                constr_fail += 1

        constr_fail, check_fail, success = [
            100*x / num_tests
            for x in (constr_fail, check_fail, success)
        ]
        print("{:.2f}% failed constructions, {:.2f}% false, {:.2f}% true".format(
            constr_fail, check_fail, success
        ))

if __name__ == "__main__":
    #datadir = "ggb-benchmark/true"
    datadir = "patrik"
    construction = Construction()
    for filename in os.listdir(datadir):
        if not filename.endswith(".txt"): continue
        construction.load(os.path.join(datadir, filename))
        construction.test()
