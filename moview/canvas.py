import numpy as np
import matplotlib.pyplot as plt

__all__ = ['Canvas']

class Canvas(list):
    def __new__(cls):
        """Generate new molecule object.
        """
        obj = list().__new__(cls)

        return obj

    def __init__(self,background='pink', dimension=3, axis='off'):
        """Generate new molecule object.
        """
        self.background = background
        self.dimension = dimension
        self.axis = axis 

        fig = plt.figure()
        if dimension==2:
            ax = plt.axes()
        elif dimension==3:
            ax = plt.axes(projection='3d')
        else:
            raise ValueError('Dimension of canvas must be two or three.')

        ax.set_facecolor(color=background)

        self.plt = plt
        self.fig = fig
        self.ax = ax

    def append(self, obj):
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
        super().append(obj)
        #if isinstance(obj, Molecule):
        #    super().append(obj)
        #elif isinstance(obj, Orbital):
        #    super().append(obj)
        #else:
        #    raise ValueError("Unknown format")

    def plot(self):
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
        
        #for obj in self:
        #    if isinstance(obj, Molecule):
        #        self.plot_atom(obj)
        #    elif isinstance(obj, Orbital):
        #        self.plot_orbital(obj)
        #    else:
        #        raise ValueError("Unknown format")
        self.plot_atom(self[0])
        self.plot_orbital(self[1],0.2)

        # Plot
        self.plt.axis('off')
        self.plt.show()
            
    def plot_atom(self,mol):
        x = []
        y = []
        z = []
        size = []
        color = []
        for atom in mol:
            x.append(atom.coordinate[0])
            y.append(atom.coordinate[1])
            z.append(atom.coordinate[2])
            size.append(atom.radius_calculated*50)
            color.append([i/255 for i in atom.color_jmol])

        self.ax.scatter3D(x, y, z, s=size,c=color)


    def plot_bond(self,mol):
        pass

    def plot_molecule(self,mol):
        self.plot_atom(mol)
        self.plot_bond(mol)

    def plot_orbital(self,orb,iso_val):
        from skimage.measure import marching_cubes
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        x, y, z = np.mgrid[-10:10:31j, -10:10:31j, -10:10:31j]
        orb = orb(x,y,z)

        # Use marching cubes to obtain the surface mesh of these ellipsoids
        verts, faces, normals, values = marching_cubes(orb, iso_val*orb.max())
        verts_, faces_, normals_, values_ = marching_cubes(orb, iso_val*orb.min())

        # Display resulting triangular mesh using Matplotlib. This can also be done
        # with mayavi (see skimage.measure.marching_cubes_lewiner docstring).

        # Fancy indexing: `verts[faces]` to generate a collection of triangles
        mesh = Poly3DCollection(verts[faces])
        mesh_ = Poly3DCollection(verts_[faces_])

        mesh.set_edgecolor((0.0, 0.0, 1.0, 1.0))
        mesh.set_facecolor((0.0, 0.0, 1.0, 0.0))
        mesh_.set_edgecolor((1.0, 0.0, 0.0, 1.0))
        mesh_.set_facecolor((1.0, 0.0, 0.0, 0.0))

        self.ax.add_collection3d(mesh)
        self.ax.add_collection3d(mesh_)
        # Centering
        xmin = verts[:,0].min()
        xmax = verts[:,0].max()
        ymin = verts[:,1].min()
        ymax = verts[:,1].max()
        zmin = verts[:,2].min()
        zmax = verts[:,2].max()
        self.ax.set_xlim(-10, 10)
        self.ax.set_ylim(-10, 10)
        self.ax.set_zlim(-10, 10)

        # Plot
        self.plt.axis('off')
        self.plt.show()
                                                             
