#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Pure Python implementations of some OpenMM objects.

For now, we provide pure Python implementations of System, Force, and some derived classes.
In addition, these objects can be constructed from their equivalent Swig proxies, and an
asSwig() method is provided to construct and return a Swig proxy for each class.

COPYRIGHT

@author John D. Chodera <jchodera@gmail.com>

EXAMPLES

First, some useful imports:

>>> import simtk.openmm as openmm # Swig wrappers for OpenMM objects
>>> import pyopenmm # pure Python OpenMM objects (containers, really)
>>> import testsystems # some test systems to play with
>>> import simtk.unit as units # unit system

We can easily create pure Python versions of objects created with simtk.openmm by passing
them to the System constructor:

>>> [system, coordinates] = testsystems.LennardJonesFluid(mm=openmm) # create a Swig wrapped object
>>> pyopenmm_system = pyopenmm.System(system) # creates pure Python forms of System and all Force objects

The Python implementation provides a complete implementation of the Swig OpenMM API, which
follows from the C++ API:

>>> (system.getNumParticles(), pyopenmm_system.getNumParticles())
(216, 216)
>>> (system.getForce(0).getParticleParameters(0)[0], pyopenmm_system.getForce(0).getParticleParameters(0)[0])
(Quantity(value=0.0, unit=elementary charge), Quantity(value=0.0, unit=elementary charge))
    
But the pure Python pyopenmm implementation is a 'thick' wrapper.  In addition to the Swig
OpenMM API, it provides a more 'pythonic' way to access information within objects:

>>> pyopenmm_system.nparticles
216
>>> pyopenmm_system.forces[0].particles[0].charge
Quantity(value=0.0, unit=elementary charge)

Each class also provides an asSwig() method for creating a Swig-wrapped object from the pure Python object:

>>> swig_force = pyopenmm_system.forces[0].asSwig()
>>> swig_system = pyopenmm_system.asSwig() # all Force objects are also converted to Swig

We don't need to do this explicitly, however, because we also provide a new Context implementation
which automatically creates Swig-wrapped objects only when necessary:

>>> integrator = pyopenmm.VerletIntegrator(1.0 * units.femtosecond) # all simtk.openmm symbols are imported first; integrator is actually a Swig-wrapped object from simtk.openmm
>>> context = pyopenmm.Context(pyopenmm_system, integrator) # create a context from a pure Python System object, automatically converting to Swig-wrapped object as needed

Pure Python objects can be deep-copied:

>>> import copy
>>> system_copy = copy.deepcopy(pyopenmm_system)

and pickled:

>>> import pickle
>>> import tempfile # for creating a temporary file for the purpose of this example
>>> file = tempfile.NamedTemporaryFile()
>>> pickle.dump(pyopenmm_system, file)

The operator '+' has been overloaded to allow System objects to be concatenated, provided it is
sensible to do so (i.e. they have the same Force objects in the same order, with the same settings):

>>> import simtk.pyopenmm.amber.amber_file_parser as amber # for reading AMBER files
>>> import os.path
>>> path = os.path.join(os.getenv('PYOPENMM_SOURCE_DIR'), 'test', 'additional-tests', 'systems', 'T4-lysozyme-L99A-implicit')
>>> receptor_system = amber.readAmberSystem(os.path.join(path, 'receptor.prmtop'), mm=pyopenmm) # create System from AMBER prmtop
>>> receptor_coordinates = amber.readAmberCoordinates(os.path.join(path, 'receptor.crd')) # read AMBER coordinates
>>> ligand_system = amber.readAmberSystem(os.path.join(path, 'ligand.prmtop'), mm=pyopenmm) # create System from AMBER prmtop
>>> ligand_coordinates = amber.readAmberCoordinates(os.path.join(path, 'ligand.crd')) # read AMBER coordinates
>>> complex_system = receptor_system + ligand_system # concatenate atoms in systems
>>> complex_coordinates = numpy.concatenate([receptor_coordinates, ligand_coordinates])

TODO

* Rework Pythonic extensions to use OpenMM API, so they can be used to extend OpenMM Swig implementation?
* Debug small discrepancies in explicit solvent energies.
* Fix System _copyDataUsingInterface and addForce methods.
* Debug tiny differences in energies when system is copied repeatedly.
* Use decorators @accepts, etc. everywhere.
* Update docstrings.
* Add more doctests for Force classes.
* Add support for CustomBondForce
* Interoperability with OpenMM App topology classes?
* Change *Info structures to encapsulated classes since we don't have to worry about network transport of these objects anymore.
* return values: Many of the OpenMM interface methods return an integer index of the particle, force, or
parameter added, which doesn't make much sense in the Pythonic way of doing things, and causes
emission of the index to the terminal output unless the argument is eaten by assignment. Can
we modify this behavior without breaking anything?
* pickling: Pickling only works if all classes are defined at top-level for the module.
We can move classes like NonbondedForce.ParticleInfo out of NonbondedForce, calling it something
like NonbondedForceParticleInfo, to allow pickling to work.  Alternatively, we can turn these
ParticleInfo classes into tuples.
* Add support for enumerating and manipulating molecules in System based on constraints and bonds.
* Add support for identifying and manipulating waters or other solvents?
* Add 'validate' method to System to make sure all masses, constraints, and Force objects have 
  consistent particle numbers and entries within bounds.

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import re
import copy
import numpy 

import simtk.unit as units
import simtk.openmm as openmm 

from exceptions import *
from decorators import *

# Import everything that we don't override directly from pyopenmm.
from simtk.openmm import *

#=============================================================================================
# System
#=============================================================================================

class System(object):
    """
    This class represents a molecular system.  The definition of a System involves
    four elements:
    
    <ol>
    <li>The set of particles in the system</li>
    <li>The forces acting on them</li>
    <li>Pairs of particles whose separation should be constrained to a fixed value</li>
    <li>For periodic systems, the dimensions of the periodic box</li>
    </ol>
    
    The particles and constraints are defined directly by the System object, while
    forces are defined by objects that extend the Force class.  After creating a
    System, call addParticle() once for each particle, addConstraint() for each constraint,
    and addForce() for each Force.

    Create a System object.

    >>> system = System()

    Add a particle.

    >>> mass = 12.0 * units.amu
    >>> system.addParticle(mass)
    0

    Add a NonbondedForce.

    >>> nonbondedForce = NonbondedForce()
    >>> system.addForce(nonbondedForce)
    0

    Create a system from a Swig proxy.

    >>> import testsystems
    >>> import pyopenmm
    >>> [system, coordinates] = testsystems.LennardJonesFluid(mm=pyopenmm)
    >>> proxy_system = system.asSwig()
    >>> system_from_proxy = System(proxy_system)

    Create a deep copy.
    
    >>> import copy
    >>> system_copy = copy.deepcopy(system_from_proxy)

    Create a Swig proxy.

    >>> proxy_system_from_proxy = system_from_proxy.asSwig()

    Construct from a pyre Python class.

    >>> system_from_system = System(system_from_proxy)
        
    """

    PYOPENMM_API_EXTENSIONS = True # signal that members of this class have pyopenmm API extensions

    def __init__(self, system=None):
        """
        Create a new System.

        If an openmm.System object is specified, it will be queried to construct the class.

        """

        # Set defaults.
        self.masses      = list() # masses[i] is a Quantity denoting the mass of particle i
        self.constraints = list() # constraints[i] is the ith ConstraintInfo entry
        self.forces      = list() # forces[i] is the ith force term
        self.periodicBoxVectors = [ units.Quantity((2.,0.,0.), units.nanometer), units.Quantity((0.,2.,0.), units.nanometer), units.Quantity((0.,0.,2.), units.nanometer) ] # periodic box vectors (only for periodic systems)
        # TODO: Store periodicBoxVectors as units.Quantity(numpy.array([3,3], numpy.float64), units.nanometer)?

        # Populate the system from a provided Swig proxy or Python system, if given.
        if system is not None:
            self._copyDataUsingInterface(self, system)
    
        return

    @classmethod
    def _copyDataUsingInterface(cls, dest, src):
        """
        Use the public interface to populate 'dest' from 'src'.
        
        """
        dest.__init__() # TODO: Do we need this?
        for index in range(src.getNumParticles()):
            mass = src.getParticleMass(index)            
            dest.addParticle(mass)
        for index in range(src.getNumConstraints()):
            args = src.getConstraintParameters(index)
            dest.addConstraint(*args)
        for index in range(src.getNumForces()):
            force = src.getForce(index)
            if not hasattr(dest, 'PYOPENMM_API_EXTENSIONS'):
                # We're copying a Python to a Swig System object.  Make a Swig copy.
                force = force.asSwig()
                dest.addForce(force)
            else:
                import pyopenmm
                force_name = type(force).__name__ # get name
                #force = globals()[force_name](force) # create pure Python force object
                ForceSubclass = getattr(pyopenmm, force_name) # get corresponding pyopenmm function
                force = ForceSubclass(force=force) # construct pure Python version
                dest.addForce(force)

        box_vectors = src.getDefaultPeriodicBoxVectors()
        # TODO: Regularize how box vectors are returned to math OpenMM interface (using Vec3?).
        #print box_vectors
        #box_vectors = [box_vectors[0].in_units_of(units.nanometers), box_vectors[1].in_units_of(units.nanometers), box_vectors[2].in_units_of(units.nanometers)] # DEBUG
        dest.setDefaultPeriodicBoxVectors(*box_vectors)
        
        return
        
    def asSwig(self):
        """
        Create a Swig proxy for this system object.

        """
        # Create an empty swig proxy object.
        system = openmm.System()
        # Copy data using interface.
        self._copyDataUsingInterface(system, self)

        return system

    def getNumParticles(self):
        """
        Get the number of particles in this System.
        
        """
        return len(self.masses)

    @accepts_compatible_units(units.amu)
    def addParticle(self, mass):
        """
        Add a particle to the System.

        @param mass   the mass of the particle (in atomic mass units)
        @return the index of the particle that was added

        """

        self.masses.append(mass);
        return len(self.masses)-1;

    def getParticleMass(self, index):
        """
        Get the mass (in atomic mass units) of a particle.
    
        @param index the index of the particle for which to get the mass

        """        
        return self.masses[index]

    @accepts_compatible_units(None, units.amu)
    def setParticleMass(self, index, mass):
        """
        Set the mass (in atomic mass units) of a particle.

        @param index  the index of the particle for which to set the mass
        @param mass   the mass of the particle

        """
        masses[index] = mass
        return        

    def getNumConstraints(self):
        """
        Get the number of distance constraints in this System.

        """
        return len(self.constraints)

    @accepts_compatible_units(None, None, units.nanometer)
    def addConstraint(self, particle1, particle2, distance):
        """
        Add a constraint to the System.
        
        @param particle1 the index of the first particle involved in the constraint
        @param particle2 the index of the second particle involved in the constraint
        @param distance  the required distance between the two particles, measured in nm
        @return the index of the constraint that was added

        """
        if particle1 not in range(self.getNumParticles()):
            raise ValueError("particle1 must be in range(0, getNumParticles())")
        if particle2 not in range(self.getNumParticles()):
            raise ValueError("particle1 must be in range(0, getNumParticles())")        
        constraint = self.ConstraintInfo(particle1, particle2, distance)
        self.constraints.append(constraint)

        return

    def getConstraintParameters(self, index):
        """
        Get the parameters defining a distance constraint.
        
        @param index     the index of the constraint for which to get parameters
        @return a tuple of (particle1, particle2, distance) for the given constraint index

        """
        constraint = self.constraints[index]
        return (constraint.particle1, constraint.particle2, constraint.distance)

    @accepts_compatible_units(None, None, None, units.nanometer)
    def setConstraintParameters(self, index, particle1, particle2, distance):
        """
        Set the parameters defining a distance constraint.
        
        @param index     the index of the constraint for which to set parameters
        @param particle1 the index of the first particle involved in the constraint
        @param particle2 the index of the second particle involved in the constraint
        @param distance  the required distance between the two particles, measured in nm

        """
        if particle1 not in range(self.getNumParticles()):
            raise ValueError("particle1 must be in range(0, getNumParticles())")
        if particle2 not in range(self.getNumParticles()):
            raise ValueError("particle1 must be in range(0, getNumParticles())")
        constraint = self.ConstraintInfo(particle1,particle2,distance)
        self.constraints[index] = constraint

        return

    def addForce(self, force):
        """
        Add a Force to the System.

        @param force   the Force object to be added
        @return        the index within the System of the Force that was added

        NOTES

        If a Swig object is specified, a pure Python deep copy will be constructed.
        If a Python object is specified, the actual object will be added, not a deep copy.

        """

        # If not pure Python, make a Python copy.
        #if not hasattr(force, 'PYOPENMM_API_EXTENSIONS'):
        #    import pyopenmm
        #    # TODO: Fix this magic to make sure we're looking up the pyopenmm force.
        #    force_name = type(force).__name__ # get name
        #    #force = globals()[force_name](force) # create pure Python force object
        #    ForceSubclass = getattr(pyopenmm, force_name) # get corresponding pyopenmm function
        #    force = ForceSubclass(force) # construct pure Python version

        # Append the force.
        self.forces.append(force)
        return len(self.forces)-1

    def getNumForces(self):
        """
        Get the number of Force objects that have been added to the System.

        """
        return len(self.forces)

    def getForce(self, index):
        """
        Get a const reference to one of the Forces in this System.

        @param index  the index of the Force to get

        """
        return self.forces[index]

    def getDefaultPeriodicBoxVectors(self):
        """
        Get the default values of the vectors defining the axes of the periodic box (measured in nm).  Any newly
        created Context will have its box vectors set to these.  They will affect
        any Force added to the System that uses periodic boundary conditions.

        Currently, only rectangular boxes are supported.  This means that a, b, and c must be aligned with the
        x, y, and z axes respectively.  Future releases may support arbitrary triclinic boxes.
        
        @returns a      the vector defining the first edge of the periodic box
        @returns b      the vector defining the second edge of the periodic box
        @returns c      the vector defining the third edge of the periodic box
     
        """        
        return self.periodicBoxVectors

    def setDefaultPeriodicBoxVectors(self, a, b, c):
        """
        Set the default values of the vectors defining the axes of the periodic box (measured in nm).  Any newly
        created Context will have its box vectors set to these.  They will affect
        any Force added to the System that uses periodic boundary conditions.

        Currently, only rectangular boxes are supported.  This means that a, b, and c must be aligned with the
        x, y, and z axes respectively.  Future releases may support arbitrary triclinic boxes.
        
        @param a      the vector defining the first edge of the periodic box
        @param b      the vector defining the second edge of the periodic box
        @param c      the vector defining the third edge of the periodic box

        """
        # TODO: Argument checking.
        self.periodicBoxVectors = [a,b,c]
    
    #==========================================================================
    # CONTAINERS
    #==========================================================================

    class ConstraintInfo(object):
        """
        Distance constraint information for particles in System object.

        """

        @accepts_compatible_units(None, None, units.nanometers)
        def __init__(self, particle1, particle2, distance):
            self.particle1 = particle1
            self.particle2 = particle2
            self.distance = distance
            return

    #==========================================================================
    # PYTHONIC EXTENSIONS
    #==========================================================================    
    
    @property    
    def nparticles(self):
        """
        The number of particles.

        """
        return len(self.masses)

    @property
    def nforces(self):
        """
        The number of force objects in the system.

        """
        return len(self.forces)
        
    @property
    def nconstraints(self):
        """
        The number of interparticle distance constraints defined.

        """
        return len(self.constraints)

    @property    
    def is_periodic(self):
        """
        True if system is periodic (contains a force with 'nonbondedMethod' set to 'CutoffPeriodic', 'Ewald', or 'PME').  False otherwise.

        """
        for force in self.forces:
            if hasattr(force, 'nonbondedMethod'):
                # Make a list of potential periodic methods.
                periodic_methods = list()
                for method in ['CutoffPeriodic', 'Ewald', 'PME']:
                    if hasattr(force, method): periodic_methods.append(getattr(force, method))
                if force.nonbondedMethod in periodic_methods: return True
            
        return False

    def __str__(self):
        """
        Return an 'informal' human-readable string representation of the System object.

        """
        
        r = ""
        r += "System object containing %d particles\n" % self.nparticles

        # Show particles.
        r += "Particle masses:\n"
        r += "%8s %24s\n" % ("particle", "mass")
        for index in range(self.getNumParticles()):
            mass = self.getParticleMass(index)
            r += "%8d %24s\n" % (index, str(mass))
        r += "\n"        

        # Show constraints.
        r += "Constraints:\n"
        r += "%8s %8s %16s\n" % ("particle1", "particle2", "distance")
        for index in range(self.getNumConstraints()):
            (particle1, particle2, distance) = self.getConstraintParameters(index)
            r += "%8d %8d %s\n" % (particle1, particle2, str(distance))
        r += "\n"
        
        # Show forces.
        r += "Forces:\n"
        for force in self.forces:
            r += str(force)
            r += "\n"
            
        return r

    def __add__(self, other):
        """
        Binary concatenation of two systems.

        The atoms of the second system appear, in order, after the atoms of the first system
        in the new combined system.

        USAGE

        combined_system = system1 + system2

        NOTES

        Both systems must have identical ordering of Force terms.
        Any non-particle settings from the first System override those of the second, if they differ.

        EXAMPLES

        Concatenate two systems in vacuum.

        >>> import testsystems
        >>> [system1, coordinates1] = testsystems.LennardJonesFluid()
        >>> system1 = System(system1) # convert to pure Python system
        >>> [system2, coordinates2] = testsystems.LennardJonesFluid()
        >>> system2 = System(system2) # convert to pure Python system        
        >>> combined_system = system1 + system2

        Check number of particles in each system.

        >>> print system1.nparticles
        216
        >>> print combined_system.nparticles
        432

        Make sure the appended Force terms are independent of the source system.
        >>> print "%.1f" % (combined_system.forces[0].particles[216].sigma / units.angstroms)
        3.4
        >>> system2.forces[0].particles[0].sigma *= -1
        >>> print "%.1f" % (combined_system.forces[0].particles[216].sigma / units.angstroms)
        3.4

        """        
        system = System(self)
        system += other
        return system
        
    def __iadd__(self, other):
        """
        Append specified system.

        USAGE

        system += additional_system

        NOTES

        Both systems must have identical ordering of Force terms.
        Any non-particle settings from the first System override those of the second, if they differ.
        
        EXAMPLES

        Append atoms from a second system to the first.
    
        >>> import testsystems
        >>> [system1, coordinates1] = testsystems.LennardJonesFluid()
        >>> system1 = System(system1) # convert to pure Python system
        >>> [system2, coordinates2] = testsystems.LennardJonesFluid()
        >>> system2 = System(system2) # convert to pure Python system        
        >>> system1 += system2
        
        Check the number of particles in each system.
        
        >>> print system2.nparticles
        216
        >>> print system1.nparticles
        432
        
        Make sure the appended Force terms are independent of the source system.
        >>> print "%.1f" % (system1.forces[0].particles[216].sigma / units.angstroms)
        3.4
        >>> system2.forces[0].particles[0].sigma *= -1
        >>> print "%.1f" % (system1.forces[0].particles[216].sigma / units.angstroms)
        3.4
        
        """
        # Make sure other system is Pythonic instance by making a deep copy.
        other = System(other)

        # Check to make sure both systems are compatible.
        if not isinstance(other, System):
            raise ValueError("both arguments must be System objects")
        if (self.nforces != other.nforces):
            raise ValueError("both System objects must have identical number of Force classes")
        for (force1,force2) in zip(self.forces, other.forces):
            if type(force1) != type(force2):
                raise ValueError("both System objects must have identical ordering of Force classes")

        # Combine systems.
        offset = self.nparticles
        self.masses += other.masses
        for constraint in other.constraints:            
            self.addConstraint(constraint.particle1+offset, constraint.particle2+offset, constraint.distance)
        for (force1, force2) in zip(self.forces, other.forces):
            force1._appendForce(force2, offset)
                
        return self
      
    def removeForce(self, force):
        """
        Remove the specified force from the System.

        ARGUMENTS

        force (Force) - the Force object to be removed

        EXAMPLES

        >>> system = System()        
        >>> nonbonded_force = NonbondedForce()
        >>> bond_force = HarmonicBondForce()
        >>> system.addForce(nonbonded_force)
        0
        >>> system.addForce(bond_force)
        1
        >>> system.removeForce(nonbonded_force)
        >>> system.removeForce(bond_force)

        """
        self.forces.remove(force)
        return

