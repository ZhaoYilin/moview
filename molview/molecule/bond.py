__all__ = ['Bond']

class Bond(list):
    
    def __init__(self, atom1, atom2, order=0):
        """_summary_

        Parameters
        ----------
        element : Element
            _description_
        coordinate : list
            _description_
        """
        self.assign_atom(atom)

    def assign_atom(self,atom):
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

