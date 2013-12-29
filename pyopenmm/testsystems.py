#!/usr/local/bin/env python

"""
Module to generate Systems and positions for simple reference molecular systems for testing.

DESCRIPTION

This module provides functions for building a number of test systems of varying complexity,
useful for testing both OpenMM and various codes based on pyopenmm.

Note that the PYOPENMM_SOURCE_DIR must be set to point to where the PyOpenMM package is unpacked.

EXAMPLES

Create a 3D harmonic oscillator.

>>> import testsystems
>>> ho = testsystems.HarmonicOscillator()
>>> system, positions = ho.system, ho.positions

See list of methods for a complete list of provided test systems.

COPYRIGHT

@author Randall J. Radmer <radmer@stanford.edu>
@author John D. Chodera <jchodera@gmail.com>

All code in this repository is released under the GNU General Public License.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

TODO

* Add units checking code to check arguments.
* Change default arguments to Quantity objects, rather than None?

"""

import os
import os.path
import sys
import numpy as np
import numpy.random
import math
import copy
import scipy.special

import simtk
import simtk.openmm as mm
import simtk.unit as units
import simtk.openmm.app as app

from repex.utils import get_data_filename

kB = units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA 

#=============================================================================================
# Thermodynamic state description
#=============================================================================================

class ThermodynamicState(object):
    """
    Data specifying a thermodynamic state obeying Boltzmann statistics.

    EXAMPLES

    Specify an NVT state for a water box at 298 K.

    >>> import simtk.unit as u
    >>> state = ThermodynamicState(temperature=298.0*u.kelvin)

    Specify an NPT state at 298 K and 1 atm pressure.

    >>> state = ThermodynamicState(temperature=298.0*u.kelvin, pressure=1.0*u.atmospheres)
    
    Note that the pressure is only relevant for periodic systems.

    """
    
    def __init__(self, temperature=None, pressure=None):
        """
        Initialize the thermodynamic state.

        OPTIONAL ARGUMENTS

        system (simtk.openmm.System) - a System object describing the potential energy function for the system (default: None)
        temperature (simtk.unit.Quantity compatible with 'kelvin') - the temperature for a system with constant temperature (default: None)
        pressure (simtk.unit.Quantity compatible with 'atmospheres') - the pressure for constant-pressure systems (default: None)

        mm (simtk.openmm API) - OpenMM API implementation to use
        cache_context (boolean) - if True, will try to cache Context objects

        """

        # Initialize.
        self.temperature = temperature
        self.pressure = pressure

        return

#=============================================================================================
# Abstract base class for test systems
#=============================================================================================

class TestSystem(object):
    """Abstract base class for test systems, demonstrating how to implement a test system.

    Parameters
    ----------
    
    Attributes
    ----------
    system : simtk.openmm.System
        Openmm system with the harmonic oscillator
    positions : list
        positions of harmonic oscillator

    Notes
    -----

    Unimplemented methods will default to the base class methods, which raise a NotImplementedException.

    Examples
    --------

    Create a test system.

    >>> testsystem = TestSystem()
    
    Retrieve System object.

    >>> system = testsystem.system

    Retrieve the positions.
    
    >>> positions = testsystem.positions

    Serialize system and positions to XML (to aid in debugging).

    >>> (system_xml, positions_xml) = testsystem.serialize()

    """
    def __init__(self, temperature=None, pressure=None):
        """Abstract base class for test system.

        Parameters
        ----------

        temperature : simtk.unit.Quantity, optional, units compatible with simtk.unit.kelvin
            The temperature of the system.

        pressure : simtk.unit.Quantity, optional, units compatible with simtk.unit.atmospheres
            The pressure of the system.

        """
        
        # Create an empty system object.
        self._system = mm.System()

        # Store positions.
        self._positions = units.Quantity(np.zeros([0,3], np.float), units.nanometers)

        # Store thermodynamic parameters.
        self._temperature = temperature
        self._pressure = pressure
        
        return

    @property
    def system(self):
        """The simtk.openmm.System object corresponding to the test system."""
        return copy.deepcopy(self._system)

    @system.setter
    def system(self, value):
        self._system = value

    @system.deleter
    def system(self):
        del self._system

    @property
    def positions(self):
        """The simtk.unit.Quantity object containing the particle positions, with units compatible with simtk.unit.nanometers."""
        return copy.deepcopy(self._positions)

    @positions.setter
    def positions(self, value):
        self._positions = value
    
    @positions.deleter
    def positions(self):
        del self._positions

    @property
    def analytical_properties(self):
        """A list of available analytical properties, accessible via 'get_propertyname(thermodynamic_state)' calls."""
        return [ method[4:] for method in dir(self) if (method[0:4]=='get_') ]

    def reduced_potential_expectation(self, state_sampled_from, state_evaluated_in):
        """Calculate the expected potential energy in state_sampled_from, divided by kB * T in state_evaluated_in.
        
        Notes
        -----
        
        This is not called get_reduced_potential_expectation because this function
        requires two, not one, inputs.
        """
        
        if hasattr(self, "get_potential_expectation"):
            U = self.get_potential_expectation(state_sampled_from)
            U_red = U / (kB * state_evaluated_in.temperature)
            return U_red
        else:
            raise AttributeError("Cannot return reduced potential energy because system lacks get_potential_expectation")

    def serialize(self):
        """Return the System and positions in serialized XML form.

        Returns
        -------
        
        system_xml : str
            Serialized XML form of System object.
            
        state_xml : str
            Serialized XML form of State object containing particle positions.

        """

        from simtk.openmm import XmlSerializer
        
        # Serialize System.
        system_xml = XmlSerializer.serialize(self._system)

        # Serialize positions via State.
        if self._system.getNumParticles() == 0:
        # Cannot serialize the State of a system with no particles.
            state_xml = None
        else:
            platform = mm.Platform.getPlatformByName('Reference')
            integrator = mm.VerletIntegrator(1.0 * units.femtoseconds)
            context = mm.Context(self._system, integrator, platform)
            context.setPositions(self._positions)
            state = context.getState(getPositions=True)
            del context, integrator
            state_xml = XmlSerializer.serialize(state)

        return (system_xml, state_xml)

    @property
    def name(self):
        """The name of the test system."""
        return self.__class__.__name__

#=============================================================================================
# 3D harmonic oscillator
#=============================================================================================

class HarmonicOscillator(TestSystem):
    """Create a 3D harmonic oscillator, with a single particle confined in an isotropic harmonic well.

    Parameters
    ----------
    K : simtk.unit.Quantity, optional, default=100.0 * units.kilocalories_per_mole/units.angstrom**2
        harmonic restraining potential
    mass : simtk.unit.Quantity, optional, default=39.948 * units.amu
        particle mass
    
    Attributes
    ----------
    system : simtk.openmm.System
        Openmm system with the harmonic oscillator
    positions : list
        positions of harmonic oscillator

    Notes
    -----

    The natural period of a harmonic oscillator is T = sqrt(m/K), so you will want to use an
    integration timestep smaller than ~ T/10.

    The standard deviation in position in each dimension is sigma = (kT / K)^(1/2)

    The expectation and standard deviation of the potential energy of a 3D harmonic oscillator is (3/2)kT.

    Examples
    --------

    Create a 3D harmonic oscillator with default parameters:

    >>> ho = HarmonicOscillator()
    >>> (system, positions) = ho.system, ho.positions

    Create a harmonic oscillator with specified mass and spring constant:

    >>> mass = 12.0 * units.amu
    >>> K = 1.0 * units.kilocalories_per_mole / units.angstroms**2
    >>> ho = HarmonicOscillator(K=K, mass=mass)
    >>> (system, positions) = ho.system, ho.positions

    Get a list of the available analytically-computed properties.

    >>> print ho.analytical_properties
    ['potential_expectation', 'potential_standard_deviation']

    Compute the potential expectation and standard deviation

    >>> import simtk.unit as u
    >>> thermodynamic_state = ThermodynamicState(temperature=298.0*u.kelvin)
    >>> potential_mean = ho.get_potential_expectation(thermodynamic_state)
    >>> potential_stddev = ho.get_potential_standard_deviation(thermodynamic_state)
    
    """
    
    def __init__(self, K=100.0 * units.kilocalories_per_mole / units.angstroms**2, mass=39.948 * units.amu, **kwargs):

        TestSystem.__init__(self, kwargs)

        # Create an empty system object.
        system = mm.System()

        # Add the particle to the system.
        system.addParticle(mass)

        # Set the positions.
        positions = units.Quantity(np.zeros([1,3], np.float32), units.angstroms)

        # Add a restrining potential centered at the origin.
        force = mm.CustomExternalForce('(K/2.0) * (x^2 + y^2 + z^2)')
        force.addGlobalParameter('K', K)
        force.addParticle(0, [])
        system.addForce(force)
        
        self.K, self.mass = K, mass
        self.system, self.positions = system, positions
        
        # Number of degrees of freedom.
        self.ndof = 3

    def get_potential_expectation(self, state):
        """Return the expectation of the potential energy, computed analytically or numerically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.
        
        Returns
        -------
        
        potential_mean : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            The expectation of the potential energy.
        
        """

        return (3./2.) * kB * state.temperature
        
    def get_potential_standard_deviation(self, state):
        """Return the standard deviation of the potential energy, computed analytically or numerically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.

        Returns
        -------
        
        potential_stddev : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            potential energy standard deviation if implemented, or else None
        
        """

        return (3./2.) * kB * state.temperature

class PowerOscillator(TestSystem):
    """Create a 3D Power oscillator, with a single particle confined in an isotropic x^b well.

    Parameters
    ----------
    K : simtk.unit.Quantity, optional, default=100.0
        harmonic restraining potential.  The units depend on the power, 
        so we accept unitless inputs and add units of the form 
        units.kilocalories_per_mole / units.angstrom ** b
    mass : simtk.unit.Quantity, optional, default=39.948 * units.amu
        particle mass
    
    Attributes
    ----------
    system : simtk.openmm.System
        Openmm system with the harmonic oscillator
    positions : list
        positions of harmonic oscillator

    Notes
    -----

    Here we assume a potential energy of the form U(x) = k * x^b.  

    By the generalized equipartition theorem, the expectation of the 
    potential energy is 3 kT / b.
    
    """
    
    def __init__(self, K=100.0, b=2.0, mass=39.948 * units.amu, **kwargs):

        TestSystem.__init__(self, kwargs)
        
        K = K * units.kilocalories_per_mole / units.angstroms ** b

        # Create an empty system object.
        system = mm.System()

        # Add the particle to the system.
        system.addParticle(mass)

        # Set the positions.
        positions = units.Quantity(np.zeros([1,3], np.float32), units.angstroms)

        # Add a restrining potential centered at the origin.
        force = mm.CustomExternalForce('(K) * (x^%d + y^%d + z^%d)' %(b, b, b))
        force.addGlobalParameter('K', K)
        force.addParticle(0, [])
        system.addForce(force)
        
        self.K, self.mass = K, mass
        self.b = b
        self.system, self.positions = system, positions
        
        # Number of degrees of freedom.
        self.ndof = 3

    def get_potential_expectation(self, state):
        """Return the expectation of the potential energy, computed analytically or numerically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.
        
        Returns
        -------
        
        potential_mean : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            The expectation of the potential energy.
        
        """

        return (3.) * kB * state.temperature / self.b

    def _get_power_expectation(self, state, n):
        """Return the power of x^n.  Not currently used"""
        b = 1.0 * self.b
        beta = (1.0 * kB * state.temperature) ** -1.
        gamma = scipy.special.gamma
        return (self.K * beta) ** (-n / b) * gamma((n + 1.) / b) / gamma(1. / b)

    @classmethod
    def reduced_potential(cls, beta, a, b, a2, b2):
        gamma = scipy.special.gamma
        reduced_u = 3 * a2 * (a * beta) ** (-b2 / b) * gamma((b2 + 1.) / b) / gamma(1. / b) * beta
        return reduced_u

#=============================================================================================
# Diatomic molecule
#=============================================================================================

class Diatom(TestSystem):
    """Create a free diatomic molecule with a single harmonic bond between the two atoms.

    Parameters
    ----------
    K : simtk.unit.Quantity, optional, default=290.1 * units.kilocalories_per_mole / units.angstrom**2
        harmonic bond potential.  default is GAFF c-c bond
    r0 : simtk.unit.Quantity, optional, default=1.550 * units.amu
        bond length.  Default is Amber GAFF c-c bond.
    constraint : bool, default=False
        if True, the bond length will be constrained
    m1 : simtk.unit.Quantity, optional, default=12.01 * units.amu
        particle1 mass
    m2 : simtk.unit.Quantity, optional, default=12.01 * units.amu
        particle2 mass
    use_central_potential : bool, optional, default=False
        if True, a soft central potential will also be added to keep the system from drifting away        
    

    Notes
    -----

    The natural period of a harmonic oscillator is T = sqrt(m/K), so you will want to use an
    integration timestep smaller than ~ T/10.

    Examples
    --------

    Create a Diatom:

    >>> diatom = Diatom()
    >>> system, positions = diatom.system, diatom.positions

    Create a Diatom with constraint in a central potential
    >>> diatom = Diatom(constraint=True, use_central_potential=True)
    >>> system, positions = diatom.system, diatom.positions

    """

    def __init__(self, 
        K=290.1 * units.kilocalories_per_mole / units.angstrom**2,
        r0=1.550 * units.angstroms, 
        m1=39.948 * units.amu,
        m2=39.948 * units.amu,
        constraint=False,
        use_central_potential=False):

        # Create an empty system object.
        system = mm.System()

        # Add two particles to the system.
        system.addParticle(m1)
        system.addParticle(m2)

        # Add a harmonic bond.
        force = mm.HarmonicBondForce()
        force.addBond(0, 1, r0, K)
        system.addForce(force)

        if constraint:
            # Add constraint between particles.
            system.addConstraint(0, 1, r0)
        
        # Set the positions.
        positions = units.Quantity(np.zeros([2,3], np.float32), units.angstroms)
        positions[1,0] = r0

        if use_central_potential:
            # Add a central restraining potential.
            Kcentral = 1.0 * units.kilocalories_per_mole / units.nanometer**2
            force = mm.CustomExternalForce('(Kcentral/2.0) * (x^2 + y^2 + z^2)')
            force.addGlobalParameter('Kcentral', Kcentral)
            force.addParticle(0, [])
            force.addParticle(1, [])    
            system.addForce(force)

        self.system, self.positions = system, positions
        self.K, self.r0, self.m1, self.m2, self.constraint, self.use_central_potential = K, r0, m1, m2, constraint, use_central_potential
        
        # Store number of degrees of freedom.
        self.ndof = 6 - 1*constraint

    def get_potential_expectation(self, state):
        """Return the expectation of the potential energy, computed analytically or numerically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.
        
        Returns
        -------
        
        potential_mean : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            The expectation of the potential energy.
        
        """

        return (self.ndof/2.) * kB * state.temperature
        
#=============================================================================================
# Constraint-coupled harmonic oscillator
#=============================================================================================

class ConstraintCoupledHarmonicOscillator(TestSystem):
    """Create a pair of particles in 3D harmonic oscillator wells, coupled by a constraint.

    Parameters
    ----------
    K : simtk.unit.Quantity, optional, default=1.0 * units.kilojoules_per_mole / units.nanometer**2
        harmonic restraining potential
    d : simtk.unit.Quantity, optional, default=1.0 * units.nanometer
        distance between harmonic oscillators.  Default is Amber GAFF c-c bond.
    mass : simtk.unit.Quantity, default=39.948 * units.amu
        particle mass
    
    Attributes
    ----------
    system : simtk.openmm.System
    positions : list

    Notes
    -----

    The natural period of a harmonic oscillator is T = sqrt(m/K), so you will want to use an
    integration timestep smaller than ~ T/10.

    Examples
    --------

    Create a constraint-coupled harmonic oscillator with specified mass, distance, and spring constant.

    >>> ccho = ConstraintCoupledHarmonicOscillator()
    >>> mass = 12.0 * units.amu
    >>> d = 5.0 * units.angstroms
    >>> K = 1.0 * units.kilocalories_per_mole / units.angstroms**2
    >>> ccho = ConstraintCoupledHarmonicOscillator(K=K, d=d, mass=mass)
    >>> system, positions = ccho.system, ccho.positions
    """

    def __init__(self, 
        K=1.0 * units.kilojoules_per_mole/units.nanometer**2, 
        d=1.0 * units.nanometer, 
        mass=39.948 * units.amu):

        # Create an empty system object.
        system = mm.System()

        # Add particles to the system.
        system.addParticle(mass)
        system.addParticle(mass)    

        # Set the positions.
        positions = units.Quantity(np.zeros([2,3], np.float32), units.angstroms)
        positions[1,0] = d

        # Add a restrining potential centered at the origin.
        force = mm.CustomExternalForce('(K/2.0) * ((x-d)^2 + y^2 + z^2)')
        force.addGlobalParameter('K', K)
        force.addPerParticleParameter('d')    
        force.addParticle(0, [0.0])
        force.addParticle(1, [d / units.nanometers])    
        system.addForce(force)

        # Add constraint between particles.
        system.addConstraint(0, 1, d)   

        # Add a harmonic bond force as well so minimization will roughly satisfy constraints.
        force = mm.HarmonicBondForce()
        K = 10.0 * units.kilocalories_per_mole / units.angstrom**2 # force constant
        force.addBond(0, 1, d, K)
        system.addForce(force)
        
        self.system, self.positions = system, positions
        self.K, self.d, self.mass = K, d, mass

#=============================================================================================
# Harmonic oscillator array
#=============================================================================================

class HarmonicOscillatorArray(TestSystem):
    """Create a 1D array of noninteracting particles in 3D harmonic oscillator wells.
    
    Parameters
    ----------
    K : simtk.unit.Quantity, optional, default=90.0 * units.kilocalories_per_mole/units.angstroms**2
        harmonic restraining potential
    d : simtk.unit.Quantity, optional, default=1.0 * units.nanometer
        distance between harmonic oscillators.  Default is Amber GAFF c-c bond.
    mass : simtk.unit.Quantity, default=39.948 * units.amu
        particle mass
    N : int, optional, default=5
        Number of harmonic oscillators
    
    Attributes
    ----------
    system : simtk.openmm.System
    positions : list

    Notes
    -----

    The natural period of a harmonic oscillator is T = sqrt(m/K), so you will want to use an
    integration timestep smaller than ~ T/10.

    Examples
    --------

    Create a constraint-coupled 3D harmonic oscillator with default parameters.

    >>> ho_array = HarmonicOscillatorArray()
    >>> mass = 12.0 * units.amu
    >>> d = 5.0 * units.angstroms
    >>> K = 1.0 * units.kilocalories_per_mole / units.angstroms**2
    >>> ccho = HarmonicOscillatorArray(K=K, d=d, mass=mass)
    >>> system, positions = ccho.system, ccho.positions
    """

    def __init__(self, K=90.0 * units.kilocalories_per_mole/units.angstroms**2,
        d=1.0 * units.nanometer,
        mass=39.948 * units.amu ,
        N=5):        

        # Create an empty system object.
        system = mm.System()

        # Add particles to the system.
        for n in range(N):
            system.addParticle(mass)

        # Set the positions for a 1D array of particles spaced d apart along the x-axis.
        positions = units.Quantity(np.zeros([N,3], np.float32), units.angstroms)
        for n in range(N):
            positions[n,0] = n*d

        # Add a restrining potential for each oscillator.
        force = mm.CustomExternalForce('(K/2.0) * ((x-x0)^2 + y^2 + z^2)')
        force.addGlobalParameter('K', K)
        force.addPerParticleParameter('x0')
        for n in range(N):
            parameters = (d*n / units.nanometers, )
            force.addParticle(n, parameters)
        system.addForce(force)

        self.system, self.positions = system, positions
        self.K, self.d, self.mass, self.N = K, d, mass, N
        self.ndof = 3*N

    def get_potential_expectation(self, state):
        """Return the expectation of the potential energy, computed analytically or numerically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.
        
        Returns
        -------
        
        potential_mean : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            The expectation of the potential energy.
        
        """

        return (self.ndof/2.) * kB * state.temperature
        
    def get_potential_standard_deviation(self, state):
        """Return the standard deviation of the potential energy, computed analytically or numerically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.

        Returns
        -------
        
        potential_stddev : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            potential energy standard deviation if implemented, or else None
        
        """

        return (self.ndof/2.) * kB * state.temperature

#=============================================================================================
# Sodium chloride FCC crystal.
#=============================================================================================

class SodiumChlorideCrystal(TestSystem):
    """Create an FCC crystal of sodium chloride.

    Each atom is represented by a charged Lennard-Jones sphere in an Ewald lattice.


    Notes
    -----

    TODO

    * Lennard-Jones interactions aren't correctly being included now, due to LJ cutoff.  Fix this by hard-coding LJ interactions?
    * Add nx, ny, nz arguments to allow user to specify replication of crystal unit in x,y,z.
    * Choose more appropriate lattice parameters and lattice spacing.

    Examples
    --------

    Create sodium chloride crystal unit.
    
    >>> crystal = SodiumChlorideCrystal()
    >>> system, positions = crystal.system, crystal.positions
    """
    def __init__(self):
        # Set default parameters (from Tinker).
        mass_Na     = 22.990 * units.amu
        mass_Cl     = 35.453 * units.amu
        q_Na        = 1.0 * units.elementary_charge 
        q_Cl        =-1.0 * units.elementary_charge 
        sigma_Na    = 3.330445 * units.angstrom
        sigma_Cl    = 4.41724 * units.angstrom
        epsilon_Na  = 0.002772 * units.kilocalorie_per_mole 
        epsilon_Cl  = 0.118 * units.kilocalorie_per_mole 

        # Create system
        system = mm.System()

        # Set box vectors.
        box_size = 5.628 * units.angstroms # box width
        a = units.Quantity(np.zeros([3]), units.nanometers); a[0] = box_size
        b = units.Quantity(np.zeros([3]), units.nanometers); b[1] = box_size
        c = units.Quantity(np.zeros([3]), units.nanometers); c[2] = box_size
        system.setDefaultPeriodicBoxVectors(a, b, c)

        # Create nonbonded force term.
        force = mm.NonbondedForce()

        # Set interactions to be periodic Ewald.
        force.setNonbondedMethod(mm.NonbondedForce.Ewald)

        # Set cutoff to be less than one half the box length.
        cutoff = box_size / 2.0 * 0.99
        force.setCutoffDistance(cutoff)
        
        # Allocate storage for positions.
        natoms = 2
        positions = units.Quantity(np.zeros([natoms,3], np.float32), units.angstroms)

        # Add sodium ion.
        system.addParticle(mass_Na)
        force.addParticle(q_Na, sigma_Na, epsilon_Na)
        positions[0,0] = 0.0 * units.angstrom
        positions[0,1] = 0.0 * units.angstrom
        positions[0,2] = 0.0 * units.angstrom
        
        # Add chloride atom.
        system.addParticle(mass_Cl)
        force.addParticle(q_Cl, sigma_Cl, epsilon_Cl)
        positions[1,0] = 2.814 * units.angstrom
        positions[1,1] = 2.814 * units.angstrom
        positions[1,2] = 2.814 * units.angstrom

        # Add nonbonded force term to the system.
        system.addForce(force)
           
        self.system, self.positions = system, positions

#=============================================================================================
# Lennard-Jones cluster
#=============================================================================================

class LennardJonesCluster(TestSystem):
    """Create a non-periodic rectilinear grid of Lennard-Jones particles in a harmonic restraining potential.


    Parameters
    ----------
    nx : int, optional, default=3
        number of particles in the x direction
    ny : int, optional, default=3
        number of particles in the y direction
    nz : int, optional, default=3
        number of particles in the z direction        
    K : simtk.unit.Quantity, optional, default=1.0 * units.kilojoules_per_mole/units.nanometer**2
        harmonic restraining potential


    Examples
    --------

    Create Lennard-Jones cluster.
    
    >>> cluster = LennardJonesCluster()
    >>> system, positions = cluster.system, cluster.positions

    Create default 3x3x3 Lennard-Jones cluster in a harmonic restraining potential.

    >>> cluster = LennardJonesCluster(nx=10, ny=10, nz=10)
    >>> system, positions = cluster.system, cluster.positions
    """
    def __init__(self, nx=3, ny=3, nz=3, K=1.0 * units.kilojoules_per_mole/units.nanometer**2):        

        # Default parameters
        mass_Ar     = 39.9 * units.amu
        q_Ar        = 0.0 * units.elementary_charge
        sigma_Ar    = 3.350 * units.angstrom
        epsilon_Ar  = 0.001603 * units.kilojoule_per_mole

        scaleStepSizeX = 1.0
        scaleStepSizeY = 1.0
        scaleStepSizeZ = 1.0

        # Determine total number of atoms.
        natoms = nx * ny * nz

        # Create an empty system object.
        system = mm.System()

        # Create a NonbondedForce object with no cutoff.
        nb = mm.NonbondedForce()
        nb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)

        positions = units.Quantity(np.zeros([natoms,3],np.float32), units.angstrom)

        atom_index = 0
        for ii in range(nx):
            for jj in range(ny):
                for kk in range(nz):
                    system.addParticle(mass_Ar)
                    nb.addParticle(q_Ar, sigma_Ar, epsilon_Ar)
                    x = sigma_Ar*scaleStepSizeX*(ii - nx/2.0)
                    y = sigma_Ar*scaleStepSizeY*(jj - ny/2.0)
                    z = sigma_Ar*scaleStepSizeZ*(kk - nz/2.0)

                    positions[atom_index,0] = x
                    positions[atom_index,1] = y
                    positions[atom_index,2] = z
                    atom_index += 1

        # Add the nonbonded force.
        system.addForce(nb)

        # Add a restrining potential centered at the origin.
        force = mm.CustomExternalForce('(K/2.0) * (x^2 + y^2 + z^2)')
        force.addGlobalParameter('K', K)
        for particle_index in range(natoms):
            force.addParticle(particle_index, [])
        system.addForce(force)

        self.system, self.positions = system, positions

#=============================================================================================
# Lennard-Jones fluid
#=============================================================================================

class LennardJonesFluid(TestSystem):
    """Create a periodic rectilinear grid of Lennard-Jones particles.    
    Parameters for argon are used by default. Cutoff is set to 3 sigma by default.
    
    Parameters
    ----------
    nx : int, optional, default=6
        number of particles in the x direction
    ny : int, optional, default=6
        number of particles in the y direction
    nz : int, optional, default=6
        number of particles in the z direction        
    mass : simtk.unit.Quantity, optional, default=39.9 * units.amu
        mass of each particle.
    sigma : simtk.unit.Quantity, optional, default=3.4 * units.angstrom
        Lennard-Jones sigma parameter
    epsilon : simtk.unit.Quantity, optional, default=0.238 * units.kilocalories_per_mole
        Lennard-Jones well depth
    cutoff : simtk.unit.Quantity, optional, default=None
        Cutoff for nonbonded interactions.  If None, defaults to 2.5 * sigma
    switch : simtk.unit.Quantity, optional, default=1.0 * units.kilojoules_per_mole/units.nanometer**2
        if specified, the switching function will be turned on at this distance (default: None)
    dispersion_correction : bool, optional, default=True
        if True, will use analytical dispersion correction (if not using switching function)

    Examples
    --------

    Create default-size Lennard-Jones fluid.

    >>> fluid = LennardJonesFluid()
    >>> system, positions = fluid.system, fluid.positions

    Create a larger 10x8x5 box of Lennard-Jones particles.

    >>> fluid = LennardJonesFluid(nx=10, ny=8, nz=5)
    >>> system, positions = fluid.system, fluid.positions

    Create Lennard-Jones fluid using switched particle interactions (switched off betwee 7 and 9 A) and more particles.

    >>> fluid = LennardJonesFluid(nx=10, ny=10, nz=10, switch=7.0*units.angstroms, cutoff=9.0*units.angstroms)
    >>> system, positions = fluid.system, fluid.positions
    """

    def __init__(self, nx=6, ny=6, nz=6, 
        mass=39.9 * units.amu, # argon
        sigma=3.4 * units.angstrom, # argon, 
        epsilon=0.238 * units.kilocalories_per_mole, # argon, 
        cutoff=None, 
        switch=False, 
        dispersion_correction=True):        

        if cutoff is None:
            cutoff = 2.5 * sigma

        charge        = 0.0 * units.elementary_charge

        scaleStepSizeX = 1.0
        scaleStepSizeY = 1.0
        scaleStepSizeZ = 1.0

        # Determine total number of atoms.
        natoms = nx * ny * nz

        # Create an empty system object.
        system = mm.System()

        # Set up periodic nonbonded interactions with a cutoff.
        if switch:
            energy_expression = "LJ * S;"
            energy_expression += "LJ = 4*epsilon*((sigma/r)^12 - (sigma/r)^6);"
            #energy_expression += "sigma = 0.5 * (sigma1 + sigma2);"
            #energy_expression += "epsilon = sqrt(epsilon1*epsilon2);"
            energy_expression += "S = (cutoff^2 - r^2)^2 * (cutoff^2 + 2*r^2 - 3*switch^2) / (cutoff^2 - switch^2)^3;"
            nb = mm.CustomNonbondedForce(energy_expression)
            nb.addGlobalParameter('switch', switch)
            nb.addGlobalParameter('cutoff', cutoff)
            nb.addGlobalParameter('sigma', sigma)
            nb.addGlobalParameter('epsilon', epsilon)
            nb.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
            nb.setCutoffDistance(cutoff)        
        else:
            nb = mm.NonbondedForce()
            nb.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)    
            nb.setCutoffDistance(cutoff)
            nb.setUseDispersionCorrection(dispersion_correction)
            
        positions = units.Quantity(np.zeros([natoms,3],np.float32), units.angstrom)

        maxX = 0.0 * units.angstrom
        maxY = 0.0 * units.angstrom
        maxZ = 0.0 * units.angstrom

        atom_index = 0
        for ii in range(nx):
            for jj in range(ny):
                for kk in range(nz):
                    system.addParticle(mass)
                    if switch:
                        nb.addParticle([])                    
                    else:
                        nb.addParticle(charge, sigma, epsilon)
                    x = sigma*scaleStepSizeX*ii
                    y = sigma*scaleStepSizeY*jj
                    z = sigma*scaleStepSizeZ*kk

                    positions[atom_index,0] = x
                    positions[atom_index,1] = y
                    positions[atom_index,2] = z
                    atom_index += 1
                    
                    # Wrap positions as needed.
                    if x>maxX: maxX = x
                    if y>maxY: maxY = y
                    if z>maxZ: maxZ = z
                    
        # Set periodic box vectors.
        x = maxX+2*sigma*scaleStepSizeX
        y = maxY+2*sigma*scaleStepSizeY
        z = maxZ+2*sigma*scaleStepSizeZ

        a = units.Quantity((x,                0*units.angstrom, 0*units.angstrom))
        b = units.Quantity((0*units.angstrom,                y, 0*units.angstrom))
        c = units.Quantity((0*units.angstrom, 0*units.angstrom, z))
        system.setDefaultPeriodicBoxVectors(a, b, c)

        # Add the nonbonded force.
        system.addForce(nb)

        self.system, self.positions = system, positions

#=============================================================================================
# Custom Lennard-Jones fluid
#=============================================================================================

class CustomLennardJonesFluid(TestSystem):
    """Create a periodic rectilinear grid of Lennard-Jones particled, but implemented via CustomBondForce rather than NonbondedForce.
    Parameters for argon are used by default. Cutoff is set to 3 sigma by default.
    
    Parameters
    ----------
    nx : int, optional, default=6
        number of particles in the x direction
    ny : int, optional, default=6
        number of particles in the y direction
    nz : int, optional, default=6
        number of particles in the z direction        
    mass : simtk.unit.Quantity, optional, default=39.9 * units.amu
        mass of each particle.
    sigma : simtk.unit.Quantity, optional, default=3.4 * units.angstrom
        Lennard-Jones sigma parameter
    epsilon : simtk.unit.Quantity, optional, default=0.238 * units.kilocalories_per_mole
        Lennard-Jones well depth
    cutoff : simtk.unit.Quantity, optional, default=None
        Cutoff for nonbonded interactions.  If None, defaults to 2.5 * sigma
    switch : simtk.unit.Quantity, optional, default=1.0 * units.kilojoules_per_mole/units.nanometer**2
        if specified, the switching function will be turned on at this distance (default: None)
    dispersion_correction : bool, optional, default=True
        if True, will use analytical dispersion correction (if not using switching function)

    Notes
    -----

    No analytical dispersion correction is included here.

    Examples
    --------

    Create default-size Lennard-Jones fluid.

    >>> fluid = CustomLennardJonesFluid()
    >>> system, positions = fluid.system, fluid.positions

    Create a larger 10x8x5 box of Lennard-Jones particles.

    >>> fluid = CustomLennardJonesFluid(nx=10, ny=8, nz=5)
    >>> system, positions = fluid.system, fluid.positions

    Create Lennard-Jones fluid using switched particle interactions (switched off betwee 7 and 9 A) and more particles.

    >>> fluid = CustomLennardJonesFluid(nx=10, ny=10, nz=10, switch=7.0*units.angstroms, cutoff=9.0*units.angstroms)
    >>> system, positions = fluid.system, fluid.positions
    """

    def __init__(self, nx=6, ny=6, nz=6, 
        mass=39.9 * units.amu, # argon
        sigma=3.4 * units.angstrom, # argon, 
        epsilon=0.238 * units.kilocalories_per_mole, # argon, 
        cutoff=None, 
        switch=False, 
        dispersion_correction=True):        

        if cutoff is None:
            cutoff = 2.5 * sigma
            
        charge        = 0.0 * units.elementary_charge            
        scaleStepSizeX = 1.0
        scaleStepSizeY = 1.0
        scaleStepSizeZ = 1.0

        # Determine total number of atoms.
        natoms = nx * ny * nz

        # Create an empty system object.
        system = mm.System()

        # Set up periodic nonbonded interactions with a cutoff.
        if switch:
            energy_expression = "LJ * S;"
            energy_expression += "LJ = 4*epsilon*((sigma/r)^12 - (sigma/r)^6);"
            energy_expression += "S = (cutoff^2 - r^2)^2 * (cutoff^2 + 2*r^2 - 3*switch^2) / (cutoff^2 - switch^2)^3;"
            nb = mm.CustomNonbondedForce(energy_expression)
            nb.addGlobalParameter('switch', switch)
            nb.addGlobalParameter('cutoff', cutoff)
            nb.addGlobalParameter('sigma', sigma)
            nb.addGlobalParameter('epsilon', epsilon)
            nb.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
            nb.setCutoffDistance(cutoff)        
        else:
            energy_expression = "4*epsilon*((sigma/r)^12 - (sigma/r)^6);"
            nb = mm.CustomNonbondedForce(energy_expression)
            nb.addGlobalParameter('sigma', sigma)
            nb.addGlobalParameter('epsilon', epsilon)
            nb.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
            nb.setCutoffDistance(cutoff)        
            
        positions = units.Quantity(np.zeros([natoms,3],np.float32), units.angstrom)

        maxX = 0.0 * units.angstrom
        maxY = 0.0 * units.angstrom
        maxZ = 0.0 * units.angstrom

        atom_index = 0
        for ii in range(nx):
            for jj in range(ny):
                for kk in range(nz):
                    system.addParticle(mass)
                    nb.addParticle([])                    
                    x = sigma*scaleStepSizeX*ii
                    y = sigma*scaleStepSizeY*jj
                    z = sigma*scaleStepSizeZ*kk

                    positions[atom_index,0] = x
                    positions[atom_index,1] = y
                    positions[atom_index,2] = z
                    atom_index += 1
                    
                    # Wrap positions as needed.
                    if x>maxX: maxX = x
                    if y>maxY: maxY = y
                    if z>maxZ: maxZ = z
                    
        # Set periodic box vectors.
        x = maxX+2*sigma*scaleStepSizeX
        y = maxY+2*sigma*scaleStepSizeY
        z = maxZ+2*sigma*scaleStepSizeZ

        a = units.Quantity((x,                0*units.angstrom, 0*units.angstrom))
        b = units.Quantity((0*units.angstrom,                y, 0*units.angstrom))
        c = units.Quantity((0*units.angstrom, 0*units.angstrom, z))
        system.setDefaultPeriodicBoxVectors(a, b, c)

        # Add the nonbonded force.
        system.addForce(nb)

        # Add long-range correction.
        if switch:
            # TODO
            pass
        else:
            volume = x*y*z
            density = natoms / volume
            per_particle_dispersion_energy = -(8./3.)*math.pi*epsilon*(sigma**6)/(cutoff**3)*density  # attraction
            per_particle_dispersion_energy += (8./9.)*math.pi*epsilon*(sigma**12)/(cutoff**9)*density  # repulsion
            energy_expression = "%f" % (per_particle_dispersion_energy / units.kilojoules_per_mole)
            force = mm.CustomExternalForce(energy_expression)
            for i in range(natoms):
                force.addParticle(i, [])
            system.addForce(force)
        
        self.system, self.positions = system, positions

#=============================================================================================
# Ideal gas
#=============================================================================================

class IdealGas(TestSystem):
    """Create an 'ideal gas' of noninteracting particles in a periodic box.

    Parameters
    ----------
    nparticles : int, optional, default=216
        number of particles
    mass : int, optional, default=39.9 * units.amu
    temperature : int, optional, default=298.0 * units.kelvin
    pressure : int, optional, default=1.0 * units.atmosphere
    volume : None
        if None, defaults to (nparticles * temperature * units.BOLTZMANN_CONSTANT_kB / pressure).in_units_of(units.nanometers**3)

    Examples
    --------

    Create an ideal gas system.

    >>> gas = IdealGas()
    >>> system, positions = gas.system, gas.positions

    Create a smaller ideal gas system containing 64 particles.

    >>> gas = IdealGas(nparticles=64)
    >>> system, positions = gas.system, gas.positions

    """

    def __init__(self, nparticles=216, mass=39.9 * units.amu, temperature=298.0 * units.kelvin, pressure=1.0 * units.atmosphere, volume=None):

        if volume is None: 
            volume = (nparticles * temperature * units.BOLTZMANN_CONSTANT_kB / pressure).in_units_of(units.nanometers**3)

        charge   = 0.0 * units.elementary_charge
        sigma    = 3.350 * units.angstrom # argon LJ 
        epsilon  = 0.0 * units.kilojoule_per_mole # zero interaction

        # Create an empty system object.
        system = mm.System()
        
        # Compute box size.
        length = volume**(1.0/3.0)
        a = units.Quantity((length,           0*units.nanometer, 0*units.nanometer))
        b = units.Quantity((0*units.nanometer,           length, 0*units.nanometer))
        c = units.Quantity((0*units.nanometer, 0*units.nanometer, length))
        system.setDefaultPeriodicBoxVectors(a, b, c)

        # Add particles.
        for index in range(nparticles):
            system.addParticle(mass)
        
        # Place particles at random positions within the box.
        # TODO: Use reproducible seed.
        # NOTE: This may not be thread-safe.

        state = np.random.get_state()
        np.random.seed(0)
        positions = units.Quantity((length/units.nanometer) * np.random.rand(nparticles,3), units.nanometer)
        np.random.set_state(state)

        self.system, self.positions = system, positions
        self.ndof = 3 * nparticles

    def get_potential_expectation(self, state):
        """Return the expectation of the potential energy, computed analytically or numerically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.
        
        Returns
        -------
        
        potential_mean : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            The expectation of the potential energy.
        
        """

        return 0.0 * units.kilojoules_per_mole
        
    def get_potential_standard_deviation(self, state):
        """Return the standard deviation of the potential energy, computed analytically or numerically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.

        Returns
        -------
        
        potential_stddev : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            potential energy standard deviation if implemented, or else None
        
        """

        return 0.0 * units.kilojoules_per_mole

    def get_kinetic_expectation(self, state):
        """Return the expectation of the kinetic energy, computed analytically or numerically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.
        
        Returns
        -------
        
        potential_mean : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            The expectation of the potential energy.
        
        """

        return (3./2.) * kB * state.temperature 
        
    def get_kinetic_standard_deviation(self, state):
        """Return the standard deviation of the kinetic energy, computed analytically or numerically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.

        Returns
        -------
        
        potential_stddev : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            potential energy standard deviation if implemented, or else None
        
        """

        return (3./2.) * kB * state.temperature 

    def get_volume_expectation(self, state):
        """Return the expectation of the volume, computed analytically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature and pressure defined
            The thermodynamic state at which the property is to be computed.
        
        Returns
        -------
        
        volume_mean : simtk.unit.Quantity compatible with simtk.unit.nanometers**3
            The expectation of the volume at equilibrium.
        
        Notes
        -----
        
        The true mean volume is used, rather than the large-N limit.

        """
        
        if not state.pressure:
            box_vectors = self.system.getDefaultPeriodicBoxVectors()
            volume = box_vectors[0][0] * box_vectors[1][1] * box_vectors[2][2] 
            return volume

        N = self._system.getNumParticles()
        return ((N+1) * units.BOLTZMANN_CONSTANT_kB * state.temperature / state.pressure).in_units_of(units.nanometers**3)
        
    def get_volume_standard_deviation(self, state):
        """Return the standard deviation of the volume, computed analytically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature and pressure defined
            The thermodynamic state at which the property is to be computed.

        Returns
        -------
        
        volume_stddev : simtk.unit.Quantity compatible with simtk.unit.nanometers**3
            The standard deviation of the volume at equilibrium.
        
        Notes
        -----
        
        The true mean volume is used, rather than the large-N limit.

        """
        
        if not state.pressure:
            return 0.0 * units.nanometers**3

        N = self._system.getNumParticles()
        return (numpy.sqrt(N+1) * units.BOLTZMANN_CONSTANT_kB * state.temperature / state.pressure).in_units_of(units.nanometers**3)
    
#=============================================================================================
# Water box
#=============================================================================================

class WaterBox(TestSystem):
    """Create a test system containing a periodic box of TIP3P water.

    Flexible bonds and angles are always added, and constraints are optional (but on by default).
    Addition of flexible bond and angle terms doesn't affect constrained dynamics, but allows for minimization to work properly.

    Parameters
    ----------
    constrain : bool, optional, default=True
        if True, will also constrain OH and HH bonds in water (default: True)    
    flexible : bool, optional, default=True,
        if True, will add harmonic OH bonds and HOH angle between
    cutoff : Quantity, optional, default=None,
        If None, defaults to box_size / 2.0 * 0.999
    nonbonded_method : default=None
    filename : str, optional, default="watbox216.pdb"
        name of file containing water positions    
    charges : bool, optional, default=True        
    
    Examples
    --------

    Create a 216-water system.
    
    >>> water_box = WaterBox()
    >>> (system, positions) = water_box.system, water_box.positions

    TODO
    ----

    * Allow size of box (either dimensions or number of waters) to be specified, replicating equilibrated waterbox to fill these dimensions.

    """

    def __init__(self, constrain=True, flexible=True, cutoff=None, nonbonded_method=None, filename=None, charges=True):
        

        # Construct filename
        if filename is None:
            filename = os.path.join(os.path.dirname(__file__), 'data', 'waterbox', 'watbox216.pdb')    
        
        # Partial atomic charges for water
        massO  = 16.0 * units.amu
        massH  =  1.0 * units.amu

        # Partial atomic charges for TIP3P water
        if charges:
            qO  = -0.8340 * units.elementary_charge
            qH  =  0.4170 * units.elementary_charge
        else:
            qO  = 0.0 * units.elementary_charge
            qH  = 0.0 * units.elementary_charge        
        
        # Lennard-Jones parameters for oxygen-oxygen interactions
        sigma   = 3.15061 * units.angstrom
        epsilon = 0.6364  * units.kilojoule_per_mole

        # Water bond and angle values
        rOH   = 0.9572  * units.angstrom
        aHOH  = 104.52  * units.degree

        # Water bond and angle spring constants.
        kOH   = 553.0 * units.kilocalories_per_mole / units.angstrom**2 # from AMBER parm96
        kHOH  = 100.0 * units.kilocalories_per_mole / units.radian**2 # from AMBER parm96

        # Distance between the two H atoms in water
        rHH = 2*rOH*units.sin(aHOH/2.0)

        def loadCoordsHOH(infile):
            """
            Load water positions from a PDB file.
            
            """
            pdbData = []
            atomNum = 0
            resNum = 0
            for line in infile:
                if line.find('HETATM')==0 or line.find('ATOM  ')==0:
                    resName=line[17:20]
                    if resName in ['HOH', 'WAT', 'SOL']:
                        atomNum+=1
                        atomName=line[12:16].strip()
                        if atomName in ['O', 'OW']:
                            resNum+=1
                        try:
                            if atomName==pdbData[-1][1] or atomName==pdbData[-2][1]:
                                raise Exception("bad water molecule near %s..." % line[:27])
                        except IndexError:
                            pass
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        pdbData.append( (atomNum, atomName, resName, resNum, ((x, y, z) * units.angstrom) ) )
                        
            return pdbData

        # Load waters from input pdb file
        infile = open(filename, 'r')
        pdbData = loadCoordsHOH(infile)
        infile.close()

        # Determine number of atoms.
        natoms = len(pdbData)
        
        # Create water system, which will inlcude force field
        # and general system parameters
        system = mm.System()

        # Create nonbonded, harmonic bond, and harmonic angle forces.
        nb = mm.NonbondedForce()
        bond = mm.HarmonicBondForce()
        angle = mm.HarmonicAngleForce()

        # Add water molecules to system
        # Note that no bond forces are used. Bond lenths are rigid
        count = 0
        positions = units.Quantity(np.zeros([natoms,3], np.float32), units.nanometer)
        for atomNum, atomName, resName, resNum, xyz in pdbData:
            if atomName in ['O', 'OW']:
                # Add an oxygen atom
                system.addParticle(massO)
                nb.addParticle(qO, sigma, epsilon)
                lastOxygen = count
            elif atomName in ['H1', 'H2', 'HW1', 'HW2']:            
                # Add an hydrogen atom
                system.addParticle(massH)
                zero_chargeprod = 0.0 * units.elementary_charge**2
                unit_sigma = 1.0 * units.angstroms
                zero_epsilon = 0.0 * units.kilocalories_per_mole            
                nb.addParticle(qH, unit_sigma, zero_epsilon)
                if (count == lastOxygen+1):
                    # For the last oxygen and hydrogen number 1:
                    if constrain: system.addConstraint(lastOxygen, count, rOH)        #O-H1
                    # Add harmonic bond force for this bond.
                    bond.addBond(lastOxygen, count, rOH, kOH)
                    # Exception: chargeProd=0.0, sigma=1.0, epsilon=0.0
                    nb.addException(lastOxygen, count, zero_chargeprod, unit_sigma, zero_epsilon)   #O-H1
                elif (count == lastOxygen+2):
                    # For the last oxygen and hydrogen number 2:
                    if constrain: system.addConstraint(lastOxygen, count, rOH)        #O-H2
                    # Add harmonic bond force for this bond.
                    bond.addBond(lastOxygen, count, rOH, kOH)
                    # Exception: chargeProd=0.0, sigma=1.0, epsilon=0.0
                    nb.addException(lastOxygen, count, zero_chargeprod, unit_sigma, zero_epsilon)   #O-H2

                    # For hydrogen number 1 and hydrogen number 2
                    if constrain: system.addConstraint(count-1, count, rHH)           #H1-H2
                    # Add harmonic angle bend.
                    angle.addAngle(count-1, lastOxygen, count, aHOH, kHOH)
                    # Exception: chargeProd=0.0, sigma=1.0, epsilon=0.0
                    nb.addException(count-1, count, zero_chargeprod, unit_sigma, zero_epsilon)      #H1-H2
                else:
                    s = "too many hydrogens:"
                    s += " atomNum=%d, resNum=%d, resName=%s, atomName=%s" % (atomNum, resNum, resName, atomName)
                    raise Exception(s)
            else:
                raise Exception("bad atom : %s" % atomName)
            for k in range(3):
                positions[count,k] = xyz[k]        
            count += 1

        # Determine box size from maximum extent.
        box_extents = units.Quantity(np.zeros([3]), units.nanometers)
        for k in range(3):
            box_extents[k] = (positions[:,k] / units.nanometers).max() * units.nanometers - (positions[:,k] / units.nanometers).min() * units.nanometers
        box_size = (box_extents / units.nanometers).max() * units.nanometers 
        
        # Set box vectors.
        a = units.Quantity(np.zeros([3]), units.nanometers); a[0] = box_size
        b = units.Quantity(np.zeros([3]), units.nanometers); b[1] = box_size
        c = units.Quantity(np.zeros([3]), units.nanometers); c[2] = box_size
        system.setDefaultPeriodicBoxVectors(a, b, c)

        # Set nonbonded cutoff.
        nb.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
        if (nonbonded_method is not None):
            nb.setNonbondedMethod(nonbonded_method)        
        if (cutoff is None) or (cutoff >= box_size / 2.0):
            cutoff = box_size / 2.0 * 0.999 # cutoff should be smaller than half the box length    
        nb.setCutoffDistance(cutoff)

        # Add force terms to system.
        system.addForce(nb)
        if flexible:
            system.addForce(bond)
            system.addForce(angle)

        self.system, self.positions = system, positions

#=============================================================================================
# Alanine dipeptide in vacuum.
#=============================================================================================

class AlanineDipeptideVacuum(TestSystem):
    """Alanine dipeptide ff96 in vacuum.
    
    Parameters
    ----------
    constraints : optional, default=simtk.openmm.app.HBonds
    
    Examples
    --------
    
    Create alanine dipeptide with constraints on bonds to hydrogen
    >>> alanine = AlanineDipeptideVacuum()
    >>> (system, positions) = alanine.system, alanine.positions
    """

    def __init__(self, constraints=app.HBonds):

        prmtop_filename = get_data_filename("data/alanine-dipeptide-gbsa/alanine-dipeptide.prmtop")
        crd_filename = get_data_filename("data/alanine-dipeptide-gbsa/alanine-dipeptide.crd")

        prmtop = app.AmberPrmtopFile(prmtop_filename)
        system = prmtop.createSystem(implicitSolvent=None, constraints=constraints, nonbondedCutoff=None)

        # Read positions.
        inpcrd = app.AmberInpcrdFile(crd_filename)
        positions = inpcrd.getPositions(asNumpy=True)

        self.system, self.positions = system, positions
#=============================================================================================
# Alanine dipeptide in implicit solvent.
#=============================================================================================

class AlanineDipeptideImplicit(TestSystem):
    """Alanine dipeptide ff96 in OBC GBSA implicit solvent.
    
    Parameters
    ----------
    constraints : optional, default=simtk.openmm.app.HBonds
    
    Examples
    --------
    
    Create alanine dipeptide with constraints on bonds to hydrogen
    >>> alanine = AlanineDipeptideImplicit()
    >>> (system, positions) = alanine.system, alanine.positions
    """

    def __init__(self, constraints=app.HBonds):

        prmtop_filename = get_data_filename("data/alanine-dipeptide-gbsa/alanine-dipeptide.prmtop")
        crd_filename = get_data_filename("data/alanine-dipeptide-gbsa/alanine-dipeptide.crd")        

        # Initialize system.
        
        prmtop = app.AmberPrmtopFile(prmtop_filename)
        system = prmtop.createSystem(implicitSolvent=app.OBC1, constraints=constraints, nonbondedCutoff=None)

        # Read positions.
        inpcrd = app.AmberInpcrdFile(crd_filename)
        positions = inpcrd.getPositions(asNumpy=True)

        self.system, self.positions = system, positions

#=============================================================================================
# Alanine dipeptide in explicit solvent
#=============================================================================================

class AlanineDipeptideExplicit(TestSystem):
    """Alanine dipeptide ff96 in TIP3P explicit solvent with PME electrostatics.
    
    Parameters
    ----------
    constraints : optional, default=simtk.openmm.app.HBonds
    rigid_water : bool, optional, default=True
    nonbondedCutoff : Quantity, optional, default=9.0 * units.angstroms
    use_dispersion_correction : bool, optional, default=True
        If True, the long-range disperson correction will be used.
    
    Examples
    --------
    
    >>> alanine = AlanineDipeptideExplicit()
    >>> (system, positions) = alanine.system, alanine.positions
    """

    def __init__(self, constraints=app.HBonds, rigid_water=True, nonbondedCutoff=9.0 * units.angstroms, use_dispersion_correction=True):

        prmtop_filename = get_data_filename("data/alanine-dipeptide-explicit/alanine-dipeptide.prmtop")
        crd_filename = get_data_filename("data/alanine-dipeptide-explicit/alanine-dipeptide.crd")        

        # Initialize system.
        
        prmtop = app.AmberPrmtopFile(prmtop_filename)
        system = prmtop.createSystem(constraints=constraints, nonbondedMethod=app.PME, rigidWater=rigid_water, nonbondedCutoff=0.9*units.nanometer)

        # Set dispersion correction use.
        forces = { system.getForce(index).__class__.__name__ : system.getForce(index) for index in range(system.getNumForces()) }
        forces['NonbondedForce'].setUseDispersionCorrection(use_dispersion_correction)

        # Read positions.
        inpcrd = app.AmberInpcrdFile(crd_filename, loadBoxVectors=True)
        positions = inpcrd.getPositions(asNumpy=True)

        # Set box vectors.
        box_vectors = inpcrd.getBoxVectors(asNumpy=True)
        system.setDefaultPeriodicBoxVectors(box_vectors[0], box_vectors[1], box_vectors[2])
        
        self.system, self.positions = system, positions

#=============================================================================================
# T4 lysozyme L99A mutant with p-xylene ligand.
#=============================================================================================

class LysozymeImplicit(TestSystem):
    """T4 lysozyme L99A (AMBER ff96) with p-xylene ligand (GAFF + AM1-BCC) in implicit OBC GBSA solvent.

    Parameters
    ----------
    flexibleConstraints : bool, optional, default=True
    shake : string, optional, default="h-bonds"
    
    Examples
    --------
    
    >>> lysozyme = LysozymeImplicit()
    >>> (system, positions) = lysozyme.system, lysozyme.positions
    """

    def __init__(self, flexibleConstraints=True, shake='h-bonds'):

        prmtop_filename = get_data_filename("data/T4-lysozyme-L99A-implicit/complex.prmtop")
        crd_filename = get_data_filename("data/T4-lysozyme-L99A-implicit/complex.crd")

        # Initialize system.
        
        prmtop = app.AmberPrmtopFile(prmtop_filename)
        system = prmtop.createSystem(implicitSolvent=app.OBC1, constraints=app.HBonds, nonbondedCutoff=None)

        # Read positions.
        inpcrd = app.AmberInpcrdFile(crd_filename)
        positions = inpcrd.getPositions(asNumpy=True)
        
        self.system, self.positions = system, positions


class SrcImplicit(TestSystem):
    """Src kinase in implicit AMBER 99sb-ildn with OBC GBSA solvent.

    Examples
    --------
    >>> src = SrcImplicit()
    >>> system, positions = src.system, src.positions
    """
    
    def __init__(self):

        pdb_filename = get_data_filename("data/src-implicit/implicit-refined.pdb")
        pdbfile = app.PDBFile(pdb_filename)

        # Construct system.
        forcefields_to_use = ['amber99sbildn.xml', 'amber99_obc.xml'] # list of forcefields to use in parameterization
        forcefield = app.ForceField(*forcefields_to_use)
        system = forcefield.createSystem(pdbfile.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)

        # Get positions.
        positions = pdbfile.getPositions()
        
        self.system, self.positions = system, positions

#=============================================================================================
# Src kinase in explicit solvent.
#=============================================================================================

class SrcExplicit(TestSystem):
    """Src kinase (AMBER 99sb-ildn) in explicit TIP3P solvent.

    Examples
    --------
    >>> src = SrcExplicit()
    >>> system, positions = src.system, src.positions

    """
    def __init__(self):

        system_xml_filename = get_data_filename("data/src-explicit/system.xml")
        state_xml_filename = get_data_filename("data/src-explicit/state.xml")        

        # Read system.
        infile = open(system_xml_filename, 'r')
        system = mm.XmlSerializer.deserialize(infile.read())
        infile.close()

        # Read state.
        infile = open(state_xml_filename, 'r')
        serialized_state = mm.XmlSerializer.deserialize(infile.read())
        infile.close()

        positions = serialized_state.getPositions()
        box_vectors = serialized_state.getPeriodicBoxVectors()
        system.setDefaultPeriodicBoxVectors(*box_vectors)
        
        self.system, self.positions = system, positions

#=============================================================================================
# Methanol box.
#=============================================================================================

class MethanolBox(TestSystem):
    """Methanol box.

    Parameters
    ----------
    flexibleConstraints : bool, optional, default=True
    shake : string, optional, default="h-bonds"
    nonbondedCutoff : Quantity, optional, default=7.0 * units.angstroms
    nonbondedMethod : str, optional, default="CutoffPeriodic"

    Examples
    --------
    
    >>> methanol_box = MethanolBox()
    >>> system, positions = methanol_box.system, methanol_box.positions
    """

    def __init__(self, flexibleConstraints=True, shake='h-bonds', nonbondedCutoff=7.0 * units.angstroms, nonbondedMethod='CutoffPeriodic'):

        system_name = 'methanol-box'
        prmtop_filename = get_data_filename("data/%s/%s.prmtop" % (system_name, system_name))
        crd_filename = get_data_filename("data/%s/%s.crd" % (system_name, system_name))
        
        # Initialize system.
        
        prmtop = app.AmberPrmtopFile(prmtop_filename)
        system = prmtop.createSystem(constraints=app.HBonds, nonbondedMethod=app.PME, rigidWater=True, nonbondedCutoff=0.9*units.nanometer)

        # Read positions.
        inpcrd = app.AmberInpcrdFile(crd_filename, loadBoxVectors=True)
        positions = inpcrd.getPositions(asNumpy=True)

        # Set box vectors.
        box_vectors = inpcrd.getBoxVectors(asNumpy=True)
        system.setDefaultPeriodicBoxVectors(box_vectors[0], box_vectors[1], box_vectors[2])
        
        self.system, self.positions = system, positions

#=============================================================================================
# Molecular ideal gas (methanol box).
#=============================================================================================

class MolecularIdealGas(TestSystem):
    """Molecular ideal gas (methanol box).

    Parameters
    ----------
    flexibleConstraints : bool, optional, default=True
    shake : string, optional, default=None
    nonbondedCutoff : Quantity, optional, default=7.0 * units.angstroms
    nonbondedMethod : str, optional, default="CutoffPeriodic"

    Examples
    --------
    
    >>> methanol_box = MolecularIdealGas()
    >>> system, positions = methanol_box.system, methanol_box.positions
    """

    def __init__(self, flexibleConstraints=True, shake=None, nonbondedCutoff=7.0 * units.angstroms, nonbondedMethod='CutoffPeriodic'):

        system_name = 'methanol-box'
        prmtop_filename = get_data_filename("data/%s/%s.prmtop" % (system_name, system_name))
        crd_filename = get_data_filename("data/%s/%s.crd" % (system_name, system_name))

        # Initialize system.
        
        prmtop = app.AmberPrmtopFile(prmtop_filename)
        reference_system = prmtop.createSystem(constraints=app.HBonds, nonbondedMethod=app.PME, rigidWater=True, nonbondedCutoff=0.9*units.nanometer)

        # Make a new system that contains no intermolecular interactions.
        system = mm.System()
            
        # Add atoms.
        for atom_index in range(reference_system.getNumParticles()):
            mass = reference_system.getParticleMass(atom_index)
            system.addParticle(mass)

        # Add constraints
        for constraint_index in range(reference_system.getNumConstraints()):
            [iatom, jatom, r0] = reference_system.getConstraintParameters(constraint_index)
            system.addConstraint(iatom, jatom, r0)

        # Copy only intramolecular forces.
        nforces = reference_system.getNumForces()
        for force_index in range(nforces):
            reference_force = reference_system.getForce(force_index)
            if isinstance(reference_force, mm.HarmonicBondForce):
                # HarmonicBondForce
                force = mm.HarmonicBondForce()
                for bond_index in range(reference_force.getNumBonds()):
                    [iatom, jatom, r0, K] = reference_force.getBondParameters(bond_index)
                    force.addBond(iatom, jatom, r0, K)
                system.addForce(force)
            elif isinstance(reference_force, mm.HarmonicAngleForce):
                # HarmonicAngleForce
                force = mm.HarmonicAngleForce()
                for angle_index in range(reference_force.getNumAngles()):
                    [iatom, jatom, katom, theta0, Ktheta] = reference_force.getAngleParameters(angle_index)
                    force.addAngle(iatom, jatom, katom, theta0, Ktheta)
                system.addForce(force)
            elif isinstance(reference_force, mm.PeriodicTorsionForce):
                # PeriodicTorsionForce
                force = mm.PeriodicTorsionForce()
                for torsion_index in range(reference_force.getNumTorsions()):
                    [particle1, particle2, particle3, particle4, periodicity, phase, k] = reference_force.getTorsionParameters(torsion_index)
                    force.addTorsion(particle1, particle2, particle3, particle4, periodicity, phase, k)
                system.addForce(force)
            else:
                # Don't add any other forces.
                pass

        # Read positions.
        inpcrd = app.AmberInpcrdFile(crd_filename, loadBoxVectors=True)
        positions = inpcrd.getPositions(asNumpy=True)

        # Set box vectors.
        box_vectors = inpcrd.getBoxVectors(asNumpy=True)
        system.setDefaultPeriodicBoxVectors(box_vectors[0], box_vectors[1], box_vectors[2])
        
        self.system, self.positions = system, positions

#=============================================================================================
# System of particles with CustomGBForce
#=============================================================================================

class CustomGBForceSystem(TestSystem):
    """A system of particles with a CustomGBForce.

    Notes
    -----

    This example comes from TestReferenceCustomGBForce.cpp from the OpenMM distribution.
    
    Examples
    --------
    
    >>> gb_system = CustomGBForceSystem()
    >>> system, positions = gb_system.system, gb_system.positions
    """

    def __init__(self):

        numMolecules = 70
        numParticles = numMolecules*2
        boxSize = 10.0 * units.nanometers

        # Default parameters
        mass     = 39.9 * units.amu
        sigma    = 3.350 * units.angstrom
        epsilon  = 0.001603 * units.kilojoule_per_mole
        cutoff   = 2.0 * units.nanometers
        
        system = mm.System()
        for i in range(numParticles):
            system.addParticle(mass)

        system.setDefaultPeriodicBoxVectors(mm.Vec3(boxSize, 0.0, 0.0), mm.Vec3(0.0, boxSize, 0.0), mm.Vec3(0.0, 0.0, boxSize))

        # Create NonbondedForce.
        nonbonded = mm.NonbondedForce()    
        nonbonded.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
        nonbonded.setCutoffDistance(cutoff)

        # Create CustomGBForce.
        custom = mm.CustomGBForce()
        custom.setNonbondedMethod(mm.CustomGBForce.CutoffPeriodic)
        custom.setCutoffDistance(cutoff)
        
        custom.addPerParticleParameter("q")
        custom.addPerParticleParameter("radius")
        custom.addPerParticleParameter("scale")

        custom.addGlobalParameter("solventDielectric", 80.0)
        custom.addGlobalParameter("soluteDielectric", 1.0)
        custom.addComputedValue("I", "step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(1/U^2-1/L^2)*(r-sr2*sr2/r)+0.5*log(L/U)/r+C);"
                                      "U=r+sr2;"
                                      "C=2*(1/or1-1/L)*step(sr2-r-or1);"
                                      "L=max(or1, D);"
                                      "D=abs(r-sr2);"
                                      "sr2 = scale2*or2;"
                                      "or1 = radius1-0.009; or2 = radius2-0.009", mm.CustomGBForce.ParticlePairNoExclusions);
        custom.addComputedValue("B", "1/(1/or-tanh(1*psi-0.8*psi^2+4.85*psi^3)/radius);"
                                      "psi=I*or; or=radius-0.009", mm.CustomGBForce.SingleParticle);
        custom.addEnergyTerm("28.3919551*(radius+0.14)^2*(radius/B)^6-0.5*138.935485*(1/soluteDielectric-1/solventDielectric)*q^2/B", mm.CustomGBForce.SingleParticle);
        custom.addEnergyTerm("-138.935485*(1/soluteDielectric-1/solventDielectric)*q1*q2/f;"
                              "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))", mm.CustomGBForce.ParticlePairNoExclusions);

        # Add particles.
        for i in range(numMolecules):
            if (i < numMolecules/2):
                charge = 1.0 * units.elementary_charge
                radius = 0.2 * units.nanometers
                scale = 0.5
                nonbonded.addParticle(charge, sigma, epsilon)
                custom.addParticle([charge, radius, scale])

                charge = -1.0 * units.elementary_charge
                radius = 0.1 * units.nanometers
                scale = 0.5
                nonbonded.addParticle(charge, sigma, epsilon)            
                custom.addParticle([charge, radius, scale]);
            else:
                charge = 1.0 * units.elementary_charge
                radius = 0.2 * units.nanometers
                scale = 0.8
                nonbonded.addParticle(charge, sigma, epsilon)
                custom.addParticle([charge, radius, scale])

                charge = -1.0 * units.elementary_charge
                radius = 0.1 * units.nanometers
                scale = 0.8
                nonbonded.addParticle(charge, sigma, epsilon)            
                custom.addParticle([charge, radius, scale]);

        system.addForce(nonbonded)
        system.addForce(custom)    

        # Place particles at random positions within the box.
        # TODO: Use reproducible random number seed.
        # NOTE: This may not be thread-safe.
        
        state = np.random.get_state()
        np.random.seed(0)
        positions = units.Quantity((boxSize/units.nanometer) * np.random.rand(numParticles,3), units.nanometer)
        np.random.set_state(state)

        self.system, self.positions = system, positions

