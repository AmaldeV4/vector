# importing essentials for vector computation
from math import sqrt, acos, isclose

# Constants for vector operations
COMPONENT_LIMIT = 3           # Maximum number of components (x, y, z)
UNIT_VECTORS = ['i', 'j', 'k'] # Unit vector symbols for display
DIRECTIONS = ['x', 'y', 'z']   # Direction names for attribute access

class vec:
    def __init__(self, *args):
        '''
        Initializes a vector with given components.
        Pads with zeros if fewer than 3 components are provided.
        Stores a formatted cartesian string for display.
        '''
        self.components = [float(n) for n in args]
        # Pad with zeros if not enough components
        if (d:=(COMPONENT_LIMIT)-len(self.components)) > 0:
            for _ in range(d):
                self.components.append(0.0)
        else:
            self.components = self.components[0:COMPONENT_LIMIT]
        # Build cartesian string (e.g., "3i + 2j - 4k")
        cartesian = [f"{magnitude}{direction}" for magnitude, direction in zip(self.components, UNIT_VECTORS) if magnitude!=0]
        self.cartesian = " + ".join(cartesian).replace("+ -", "- ")

    @classmethod
    def from_input(cls, vector_name: str = 'vector'):
        '''
        Creates a vector from user input.
        Accepts space-separated values or math expressions for components.
        Example: "1 sqrt(3) 2"
        '''
        comps = [float(eval(c)) for c in input(f"Enter components of {vector_name}: ").split()]
        return cls(*comps)

    @property
    def magnitude(self):
        '''
        Returns the magnitude (length) of the vector, rounded to 2 decimals.
        '''
        return round(sqrt(sum([c*c for c in self.components])), 2)

    @property
    def unit(self):
        '''
        Returns the unit vector (direction only) of the vector.
        Raises ValueError for zero vector.
        '''
        mag = self.magnitude
        if mag == 0:
            raise ValueError("Zero vector has no unit vector")
        return self.__class__(*[round(c/mag, 2) for c in self.components])

    def dot(self, B):
        '''
        Returns the dot product of self and another vector B.
        '''
        product = [p*q for p, q in zip(self.components, B.components)]
        return sum(product)

    def cross(self, B):
        '''
        Returns the cross product of self and another vector B as a new vec.
        '''
        a1, a2, a3 = self.components
        b1, b2, b3 = B.components
        return vec(
            a2*b3 - a3*b2,
            a3*b1 - a1*b3,
            a1*b2 - a2*b1
        )

    def theta(self, B):
        '''
        Returns the angle (in radians) between self and another vector B.
        '''
        modA = self.magnitude
        modB = B.magnitude
        AdotB = self.dot(B)
        angle = acos(AdotB/(modA*modB))
        return angle
    
    def __repr__(self):
        '''
        Returns the cartesian string representation of the vector.
        Example: "3i + 2j - 4k"
        '''
        return self.cartesian

    def __add__(self, v): 
        '''
        Adds two vectors component-wise.
        '''
        return vec(*[a+b for a,b in zip(self.components, v.components)])
    
    def __sub__(self, v):
        '''
        Subtracts vector v from self component-wise.
        '''
        return vec(*[a-b for a,b in zip(self.components, v.components)])
    
    def __mul__(self, v):
        '''
        If v is a vector, returns dot product.
        If v is a scalar, scales the vector.
        '''
        if isinstance(v, vec):
            return self.dot(v)
        else:
            return vec(*[c*v for c in self.components])
    def __rmul__(self, v):
        '''
        Supports scalar multiplication from the left.
        '''
        return self.__mul__(v)
    
    def __xor__(self, v):
        '''
        Returns the cross product using ^ operator.
        '''
        if not isinstance(v, vec):
            return NotImplemented
        return self.cross(v)
        
    def __truediv__(self, n):
        '''
        Divides the vector by a scalar n (element-wise, rounded to 2 decimals).
        '''
        if isinstance(n, vec):
            return NotImplemented
        dividedComponents = [round(c/n, 2) for c in self.components]
        return vec(*dividedComponents)
    def __floordiv__(self, n):
        '''
        Divides the vector by a scalar n using floor division.
        '''
        if isinstance(n, vec):
            return NotImplemented
        dividedComponents = [c//n for c in self.components]
        return vec(*dividedComponents)
 
    def __matmul__(self, v):
        '''
        Returns the angle between two vectors using @ operator.
        '''
        return self.theta(v)
 
    def __eq__(self, v):
        '''
        Checks if two vectors are equal (all components are close).
        '''
        if not isinstance(v, vec):
            return NotImplemented
        return all(isclose(a, b, rel_tol=1e-9, abs_tol=1e-9) for a, b in zip(self.components, v.components))

    def __getitem__(self, i):
        '''
        Allows indexing to access components (e.g., v[0] for x component).
        '''
        return self.components[i]

    def __getattr__(self, name):
        '''
        Allows access to components by direction (e.g., v.x, v.y, v.z).
        '''
        if name in DIRECTIONS:
            idx = DIRECTIONS.index(name)
            return self.components[idx]
        raise AttributeError(f"'Vector' has no attribute '{name}'")

    def __len__(self):
        '''
        Returns the number of components in the vector.
        '''
        return len(self.components)