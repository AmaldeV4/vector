# importing essentials for vector computation
from math import pi, sqrt, acos, isclose

# Constants
COMPONENT_LIMIT = 3
UNIT_VECTORS = ['i', 'j', 'k']
DIRECTIONS = ['x', 'y', 'z']

class vec:
    def __init__(self, *args):
        '''
        initializes vector with components passed as arguments
        '''
        self.components = [float(n) for n in args]
        if (d:=(COMPONENT_LIMIT)-len(self.components)) > 0:
            for _ in range(d):
                self.components.append(0.0)
        else:
            self.components = self.components[0:COMPONENT_LIMIT]
        cartesian = [f"{magnitude}{direction}" for magnitude, direction in zip(self.components, UNIT_VECTORS) if magnitude!=0]
        self.cartesian = " + ".join(cartesian).replace("+ -", "- ")

    @classmethod
    def from_input(cls, vector_name: str = 'vector'):
        '''
        creates a vector with coorinate x, y, z separated by spaces.
        '''
        comps = [float(eval(c)) for c in input(f"Enter components of {vector_name}: ").split()]
        return cls(*comps)

    @property
    def magnitude(self):
        return round(sqrt(sum([c*c for c in self.components])), 2)

    def dot(self, B):
        ''' Returns the dot product of two vectors '''
        product = [p*q for p, q in zip(self.components, B.components)]
        return sum(product)

    def cross(self, B):
        ''' Returns the cross product of two vectors '''
        a1, a2, a3 = self.components
        b1, b2, b3 = B.components
        return vec(
            a2*b3 - a3*b2,
            a3*b1 - a1*b3,
            a1*b2 - a2*b1
        )

    def theta(self, B):
        '''Returns the angle between two vectors'''
        modA = self.magnitude
        modB = B.magnitude
        AdotB = self.dot(B)
        angle = acos(AdotB/(modA*modB))*180/pi
        precision = 3 # decimal digits in angle
        return round(angle, precision)

    # returns the cartesian coordinates when vector object is passed into print()
    def __repr__(self):
        '''Returns the cartesian of a vector when vector instance is passed into print()'''
        return self.cartesian

    # addition and subtraction of vectors via + and - operators
    def __add__(self, v): 
        return vec(*[a+b for a,b in zip(self.components, v.components)])
    
    def __sub__(self, v):
        return vec(*[a-b for a,b in zip(self.components, v.components)])
    
    # dot product via * operator
    def __mul__(self, v):
        if isinstance(v, vec): # if v is a scalar, then scale the vector
            return self.dot(v)
        else:
            return vec(*[c*v for c in self.components]) # if v is vector, return the dot product
    def __rmul__(self, v):
        return self.__mul__(v)
    
    # cross product via ^ operator
    def __xor__(self, v):
        if not isinstance(v, vec):
            return NotImplemented
        return self.cross(v)
    
    # division by scalar
    def __truediv__(self, n):
        if isinstance(n, vec):
            return NotImplemented
        dividedComponents = [round(c/n, 2) for c in self.components]
        return vec(*dividedComponents)
    def __floordiv__(self, n):
        if isinstance(n, vec):
            return NotImplemented
        dividedComponents = [c//n for c in self.components]
        return vec(*dividedComponents)
    # defines equality of two vectors. Two vectors are equal if all of their components are equal
    def __eq__(self, v):
        if not isinstance(v, vec):
            return NotImplemented
        return all(isclose(a, b, rel_tol=1e-9, abs_tol=1e-9) for a, b in zip(self.components, v.components))

    # returns the component of vector in position i. v[0] => first component
    def __getitem__(self, i):
        return self.components[i]

    # returns x,y,z components of vector v => v.x, v.y, v.z
    def __getattr__(self, name):
        if name in DIRECTIONS:
            idx = DIRECTIONS.index(name)
            return self.components[idx]
        raise AttributeError(f"'Vector' has no attribute '{name}'")

def resultant(vecA: vec, vecB: vec):
    R_components = [a+b for a, b in zip(vecA.components, vecB.components)]
    R = sqrt(sum([c*c for c in R_components]))
    return R
