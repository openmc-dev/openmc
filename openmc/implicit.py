from __future__ import annotations
from abc import ABC, abstractmethod

import lxml.etree as ET
import numpy as np

# ------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------

def _to_function(f: int | float | ImplicitFunction) -> ImplicitFunction:
    return Constant(float(f)) if isinstance(f, (int, float)) else f

# ------------------------------------------------------------------
# Main ImplicitFunction Abstract class
# ------------------------------------------------------------------

class ImplicitFunction(ABC):

    @abstractmethod
    def __repr__(self) -> str: ...

    @abstractmethod
    def evaluate(self, point) -> float: ...

    @abstractmethod
    def to_xml_element(self, _cached: list[int] | None = None) -> ET.Element: ...

    @classmethod
    def from_xml_element(cls, element: ET.Element, _cached: dict | None = None) -> ImplicitFunction:
        if _cached is None: _cached = {}
        tag    = element.tag
        attrib = element.attrib

        # Handle cache reference before parsing children (no children to parse)
        if tag == "from_cache":
            return _cached[int(attrib["id"])]

        # Recursively parse children for all other tags
        children = [cls.from_xml_element(child, _cached) for child in element]

        # Register a new cached node and return it
        if tag == "to_cache":
            node = Cached(children[0])
            _cached[int(attrib["id"])] = node
            return node

        dispatch = {
            "x":        lambda: X(),
            "y":        lambda: Y(),
            "z":        lambda: Z(),
            "constant": lambda: Constant(float(attrib["value"])),
            "add":      lambda: Add(*children),
            "neg":      lambda: Neg(*children),
            "sub":      lambda: Sub(*children),
            "scale":    lambda: Scale(children[0], float(attrib["value"])),
            "mul":      lambda: Mul(*children),
            "div":      lambda: Div(*children),
            "pow":      lambda: Pow(children[0], int(attrib["value"])),
            "sin":      lambda: Sin(*children),
            "cos":      lambda: Cos(*children),
            "sqrt":     lambda: Sqrt(*children),
            "exp":      lambda: Exp(*children),
            "log":      lambda: Log(*children),
            "abs":      lambda: Abs(*children),
        }

        if tag not in dispatch:
            raise ValueError(f"Unknown tag '{tag}'")

        return dispatch[tag]()

    # ------------------------------------------------------------------
    # Operator overloading — enables natural Python expression syntax
    # e.g.  Sin(X()) * Cos(Z()) + Constant(1)
    # ------------------------------------------------------------------

    def __add__(self, other: float | ImplicitFunction) -> ImplicitFunction:
        return Add(_to_function(self), _to_function(other))
    
    def __radd__(self, other: float | ImplicitFunction) -> ImplicitFunction:
        return Add(_to_function(other), _to_function(self))
    
    def __neg__(self) -> ImplicitFunction:
        return Neg(_to_function(self))
    
    def __sub__(self, other: float | ImplicitFunction) -> ImplicitFunction:
        return Sub(_to_function(self), _to_function(other))
    
    def __rsub__(self, other: float | ImplicitFunction) -> ImplicitFunction:
        return Sub(_to_function(other), _to_function(self))
    
    def __truediv__(self, other: float | ImplicitFunction) -> ImplicitFunction:
        return Div(_to_function(self), _to_function(other))
    
    def __rtruediv__(self, other: float | ImplicitFunction) -> ImplicitFunction:
        return Div(_to_function(other), _to_function(self))
    
    def __pow__(self, exp: int) -> ImplicitFunction:
        if not isinstance(exp, (int)):
            raise TypeError(f"Pow exponent must be an int, got {type(exp)}")
        if exp <= 0.:
            raise TypeError(f"Pow exponent must be strictly positive, got {exp}")
        return Pow(_to_function(self), exp)

    def __mul__(self, other: float | ImplicitFunction) -> ImplicitFunction:
        if isinstance(other, (int, float)):
            return Scale(self, float(other))
        return Mul(self, other)
    
    def __rmul__(self, other: float | ImplicitFunction) -> ImplicitFunction:
        if isinstance(other, (int, float)):
            return Scale(self, float(other))

# ---------------------------------------------------------------------------
# Terminal nodes
# ---------------------------------------------------------------------------

class X(ImplicitFunction):
    def __repr__(self): return "X"
    def evaluate(self, point): return point[0]
    def to_xml_element(self, _cached=None): return ET.Element("x")

class Y(ImplicitFunction):
    def __repr__(self): return "Y"
    def evaluate(self, point): return point[1]
    def to_xml_element(self, _cached=None): return ET.Element("y")

class Z(ImplicitFunction):
    def __repr__(self): return "Z"
    def evaluate(self, point): return point[2]
    def to_xml_element(self, _cached=None): return ET.Element("z")

class Constant(ImplicitFunction):
    def __repr__(self): return f"{self.value}"
    def __init__(self, value: float) -> None:
        self.value = float(value)
    def evaluate(self, point): return self.value
    def to_xml_element(self, _cached=None):
        element = ET.Element("constant")
        element.set("value", str(self.value))
        return element

# ---------------------------------------------------------------------------
# Operations
# ---------------------------------------------------------------------------

class Add(ImplicitFunction):

    def __repr__(self): return f"{self.f} + {self.g}"
    def __init__(self, f:ImplicitFunction, g:ImplicitFunction):
        self.f = f
        self.g = g
    def evaluate(self, point): return self.f.evaluate(point) + self.g.evaluate(point)
    def to_xml_element(self, _cached=None):
        if _cached is None: _cached = []
        element = ET.Element("add")
        element.append(self.f.to_xml_element(_cached))
        element.append(self.g.to_xml_element(_cached))
        return element

class Neg(ImplicitFunction):
    def __repr__(self): return f"-{self.f}"
    def __init__(self, f:ImplicitFunction):
        self.f = f
    def evaluate(self, point): return -self.f.evaluate(point)
    def to_xml_element(self, _cached=None):
        if _cached is None: _cached = []
        element = ET.Element("neg")
        element.append(self.f.to_xml_element(_cached))
        return element

class Sub(ImplicitFunction):
    def __repr__(self): return f"{self.f} - {self.g}"
    def __init__(self, f:ImplicitFunction, g:ImplicitFunction):
        self.f = f
        self.g = g
    def evaluate(self, point): return self.f.evaluate(point) - self.g.evaluate(point)
    def to_xml_element(self, _cached=None):
        if _cached is None: _cached = []
        element = ET.Element("sub")
        element.append(self.f.to_xml_element(_cached))
        element.append(self.g.to_xml_element(_cached))
        return element

class Scale(ImplicitFunction):
    def __repr__(self): return f"{self.scalar} * {self.f}"
    def __init__(self, f:ImplicitFunction, scalar: float):
        self.f = f
        self.scalar = float(scalar)
    def evaluate(self, point): return self.scalar * self.f.evaluate(point)
    def to_xml_element(self, _cached=None):
        if _cached is None: _cached = []
        element = ET.Element("scale")
        element.set("value", str(self.scalar))
        element.append(self.f.to_xml_element(_cached))
        return element
    
class Mul(ImplicitFunction):
    def __repr__(self): return f"{self.f} * {self.g}"
    def __init__(self, f:ImplicitFunction, g:ImplicitFunction):
        self.f = f
        self.g = g
    def evaluate(self, point): return self.f.evaluate(point) * self.g.evaluate(point)
    def to_xml_element(self, _cached=None):
        if _cached is None: _cached = []
        element = ET.Element("mul")
        element.append(self.f.to_xml_element(_cached))
        element.append(self.g.to_xml_element(_cached))
        return element
    
class Div(ImplicitFunction):
    def __repr__(self): return f"{self.f} / {self.g}"
    def __init__(self, f:ImplicitFunction, g:ImplicitFunction):
        self.f = f
        self.g = g
    def evaluate(self, point): 
        v = self.g.evaluate(point)
        if v == 0.: raise ValueError(f"Denominator value cannot be zero.")
        return self.f.evaluate(point) / v
    def to_xml_element(self, _cached=None):
        if _cached is None: _cached = []
        element = ET.Element("div")
        element.append(self.f.to_xml_element(_cached))
        element.append(self.g.to_xml_element(_cached))
        return element
    

class Pow(ImplicitFunction):
    def __repr__(self): return f"{self.f} ** {self.exp}"
    def __init__(self, f:ImplicitFunction, exp: int):
        if not isinstance(exp, int) or exp <= 0:
            raise TypeError(
                f"Pow exponent must be a strictly positive integer, got {exp!r}. "
                f"For exp=-1 use Div, for exp=0.5 use Sqrt.")
        self.f   = f
        self.exp = exp
    def evaluate(self, point): return self.f.evaluate(point) ** self.exp
    def to_xml_element(self, _cached=None):
        if _cached is None: _cached = []
        element = ET.Element("pow")
        element.set("value", str(self.exp))
        element.append(self.f.to_xml_element(_cached))
        return element
    
# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------

class Sin(ImplicitFunction):
    def __repr__(self): return f"Sin({self.arg})"
    def __init__(self, arg:ImplicitFunction):
        self.arg = arg
    def evaluate(self, point): return np.sin(self.arg.evaluate(point))
    def to_xml_element(self, _cached=None):
        if _cached is None: _cached = []
        element = ET.Element("sin")
        element.append(self.arg.to_xml_element(_cached))
        return element

class Cos(ImplicitFunction):
    def __repr__(self): return f"Cos({self.arg})"
    def __init__(self, arg:ImplicitFunction):
        self.arg = arg
    def evaluate(self, point): return np.cos(self.arg.evaluate(point))
    def to_xml_element(self, _cached=None):
        if _cached is None: _cached = []
        element = ET.Element("cos")
        element.append(self.arg.to_xml_element(_cached))
        return element
    
class Sqrt(ImplicitFunction):
    def __repr__(self): return f"Sqrt({self.arg})"
    def __init__(self, arg:ImplicitFunction):
        self.arg = arg
    def evaluate(self, point):
        v = self.arg.evaluate(point)
        if v < 0.: raise ValueError(f"Sqrt input value cannot be negative but is: '{v}'.")
        return np.sqrt(v)
    def to_xml_element(self, _cached=None):
        if _cached is None: _cached = []
        element = ET.Element("sqrt")
        element.append(self.arg.to_xml_element(_cached))
        return element

class Exp(ImplicitFunction):
    def __repr__(self): return f"Exp({self.arg})"
    def __init__(self, arg:ImplicitFunction):
        self.arg = arg
    def evaluate(self, point): return np.exp(self.arg.evaluate(point))
    def to_xml_element(self, _cached=None):
        if _cached is None: _cached = []
        element = ET.Element("exp")
        element.append(self.arg.to_xml_element(_cached))
        return element

class Log(ImplicitFunction):
    def __repr__(self): return f"Log({self.arg})"
    def __init__(self, arg:ImplicitFunction):
        self.arg = arg
    def evaluate(self, point):
        v = self.arg.evaluate(point)
        if v < 0.: raise ValueError(f"Log input value cannot be negative but is: '{v}'.")
        return np.log(v)
    def to_xml_element(self, _cached=None):
        if _cached is None: _cached = []
        element = ET.Element("log")
        element.append(self.arg.to_xml_element(_cached))
        return element

class Abs(ImplicitFunction):
    def __repr__(self): return f"|{self.arg}|"
    def __init__(self, arg:ImplicitFunction):
        self.arg = arg
    def evaluate(self, point): return np.abs(self.arg.evaluate(point))
    def to_xml_element(self, _cached=None):
        if _cached is None: _cached = []
        element = ET.Element("abs")
        element.append(self.arg.to_xml_element(_cached))
        return element

# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------

class Cached(ImplicitFunction):
    """
    Marks a subexpression for memoisation in C++.

    The Python object identity (id()) is used to detect shared nodes during
    XML serialisation. This means the SAME Python object must be reused
    wherever you want sharing to occur — do NOT construct a new Cached()
    at each use site.

    Correct — one object, two references:
        cx = Cached(2 * np.pi * X())
        f  = Sin(cx) * Cos(cx)          # cx serialised once as <to_cache>

    Wrong — two objects, same expression:
        f  = Sin(Cached(2 * np.pi * X())) * Cos(Cached(2 * np.pi * X()))
    """
    def __repr__(self): return f"@[{self.f}]"
    def __init__(self, f:ImplicitFunction):
        self.f = f
    def evaluate(self, point): return self.f.evaluate(point)
    def to_xml_element(self, _cached=None):
        if _cached is None: _cached = []
        key = id(self)
        if key in _cached:
            node = ET.Element("from_cache")
            node.set("id", str(_cached.index(key)))
            return node
        else:
            _cached.append(key)
            xml_id = len(_cached) - 1
            node = ET.Element("to_cache")
            node.set("id", str(xml_id))
            node.append(self.f.to_xml_element(_cached))
            return node