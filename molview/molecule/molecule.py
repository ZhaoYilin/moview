from moha.molecule import periodic
from moha.molecule.atom import Atom
from moha.basis.basis_set import BasisSet
from moha.io.log import log

import numpy as np

__all__ = ['Molecule']

class Molecule(list):
    """Molecule Class.

    Attributes
    ----------
    pg: bool or str
        Point group symmetry of the molecule.

    Property
    --------
    E_nuc(self)
        Nucelar repulsion energy.

    center_of_mass(self)
        Center of mass for molecular.

    Methods
    -------
    assign_point_group(self,pg)
        Assign point group symmetry to the molecule.

    bond_length(self,i,j)
        Calculate the the distances between atom i and j.

    bond_angle(self,i,j,k)
        Calculate the bond angle between atoms i,j,k.
 
    dihedral_angle(self,A,B,C,D)
        Calculate dihedral angle for atom connectivit A,B,C,D.
    """
    def __new__(cls):
        """Generate new molecule object.
        """
        obj = list().__new__(cls)

        return obj
 
    @classmethod
    def build(cls,geometry,pg=False):
        """Generate N electron basis set.

        Parameters
        ----------
        geometry : str or iterator
            Geometry object.

        pg: bool or str
            Point group symmetry of the molecule.

        Raises
        ------
        TypeError 
            If geometry is unknow format.
        """
        log('Molecule Section'.format())
        log.hline()

        molecule = Molecule()
        molecule.assign_point_group(pg)
        if isinstance(geometry, str):
            if geometry.endswith('.xyz'):
                from moha.molecule import xyz
                atoms = xyz.load(geometry)

        elif isinstance(geometry,(list,tuple,dict)):
            from moha.molecule import iterator
            atoms = iterator.load(geometry)

        else:
            raise ValueError("Unknown format")

        log('{0:<20s} {1:<20s} {2:<20s} {3:<20s}'.format('Element', 'X', 'Y','Z'))
        log.hline()
        for atom in atoms:
            molecule.append(atom)
            log('{0:<20s} {1:<20f} {2:<20f} {3:<20f}'.format(atom.symbol, 
            atom.coordinate[0], atom.coordinate[1],atom.coordinate[2]))
        log.blank()


        return molecule

    def append(self,atom):
        """Add atom to the molecule.

        Parameters
        ----------
        atom : Atom
            Atom instance.
        
        Raises
        ------
        TypeError
            If atom is not an Atom instance.
        """
        if not isinstance(atom, Atom):
            raise TypeError("Atom must be a Atom instance")
        super().append(atom)

    def assign_point_group(self,pg):
        """Assign point group symmetry to the molecule.

        Parameters
        ----------
        pg: bool or str
            Point group symmetry of the molecule.

        Raises
        ------
        TypeError
            If pg is not bool or str.
        """
        if not isinstance(pg, bool or str):
            raise TypeError("Parameters pg must be bool or str.")
        self.pg = pg
                                
    @property
    def e_nuc(self):
        """Nucelar repulsion energy.

        Return
        ------
        energy: float
            Nuclear repulsion energy.
        """
        energy = 0.
        for i,A in enumerate(self):
            for j,B in enumerate(self):
                if j>i:
                    energy += (A.number*B.number)/self.bond_length(i,j)
        return energy

    def bond_length(self,i,j):
        """Calculate the the distances between atom i and j.

        Parameters
        ----------
        i : int
            Index of first atom.
        
        j : int
            Index of second atom.
            
        Raises
        ------
        TypeError
            If index of atom is not int.
        ValueError
            If index of atom is not zero or positive number.            
            If index of atom is not samller than number of atoms.            
        """
        for item in (i,j):
            if not isinstance(item,int):
                raise TypeError("Index of atom must be int")
            if item < 0:
                raise ValueError("Index of atom must be none negative number")
            if item >= len(self):
                raise ValueError("Index of atom must be smaller than number of atoms")
        A = np.array(self[i].coordinate)
        B = np.array(self[j].coordinate)
        return np.linalg.norm(A-B)

    def bond_angle(self,i,j,k):
        """Calculate the bond angle between atoms i,j,k.
           Where atom j is the central atom.

        Parameters
        ----------
        i : int
            Index of first atom.
        
        j : int
            Index of second atom.
            
        k : int
            Index of third atom.

        Raises
        ------
        TypeError
            If index of atom is not int.
        ValueError
            If index of atom is not zero or positive number.            
            If index of atom is not samller than number of atoms.            
        """
        for item in (i,j,k):
            if not isinstance(item,int):
                raise TypeError("Index of atom must be int")
            if item < 0:
                raise ValueError("Index of atom must be none negative number")
            if item >= len(self):
                raise ValueError("Index of atom must be smaller than number of atoms")
        A = np.array(self[i].coordinate)
        B = np.array(self[j].coordinate)
        C = np.array(self[k].coordinate)
        return np.linalg.norm(A-B)

    def dihedral_angle(self,A,B,C,D):
        """Calculate dihedral angle for atom connectivit A,B,C,D.

        Parameters
        ----------
        i : int
            Index of first atom.
        
        j : int
            Index of second atom.
            
        k : int
            Index of third atom.

        l : int
            Index of fourth atom.

        Raises
        ------
        TypeError
            If index of atom is not int.
        ValueError
            If index of atom is not zero or positive number.            
            If index of atom is not samller than number of atoms.            
        """
        for item in (i,j,k):
            if not isinstance(item,int):
                raise TypeError("Index of atom must be int")
            if item < 0:
                raise ValueError("Index of atom must be none negative number")
            if item >= len(self):
                raise ValueError("Index of atom must be smaller than number of atoms")
        A = np.array(self[i].coordinate)
        B = np.array(self[j].coordinate)
        C = np.array(self[k].coordinate)
        D = np.array(self[l].coordinate)

    @property
    def center_of_mass(self):
        """  This function calculates the center of mass of a given data array.
        The format of the data array is a two dimensional array, which contains
        the x, y, z Cartesion coordinates and the corresponding mass values for
        each x, y, z coordinates.

        It returns a 1-D data array.
        """
        geometry = []
        mass = [] 
        for atom in self:
            geometry.append(atom.coordinate)
            mass.append(atom.mass)
        return np.average(geometry, axis=0, weights=mass)

    @property
    def moment_of_inertia(self):
        """  This function calculates the Elements of inertia tensor for the
        center-moved coordinates.
    
        The structure of the array is a two dimensional array, which contains
        the mass of elements (the first column) and their corresponded x, y, z
        Cartesion coordinates.

        It returns an inertia tensor.
        """
        geometry = []
        mass = [] 
        for atom in self:
            geometry.append(atom.coordinate)
            mass.append(atom.mass)
        geometry = np.array(geometry)
        mass = np.array(mass)
        com_coord = np.average(geometry, axis=0, weights=mass)
        new_coord = geometry - com_coord
        I_xx = (mass * np.sum(np.square(new_coord[:,1:3:1]),axis=1)).sum()
        I_yy = (mass * np.sum(np.square(new_coord[:,0:3:2]),axis=1)).sum()
        I_zz = (mass * np.sum(np.square(new_coord[:,0:2:1]),axis=1)).sum()
        I_xy = (-1 * mass * np.prod(new_coord[:,0:2:1],axis=1)).sum()
        I_yz = (-1 * mass * np.prod(new_coord[:,1:3:1],axis=1)).sum()
        I_xz = (-1 * mass * np.prod(new_coord[:,0:3:2],axis=1)).sum()
        I = np.array([[I_xx, I_xy, I_xz],
		    [I_xy, I_yy, I_yz],
		    [I_xz, I_yz, I_zz]])
        return I
