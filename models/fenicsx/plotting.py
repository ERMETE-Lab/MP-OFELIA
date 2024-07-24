import dolfinx.plot
from dolfinx import geometry
from dolfinx.fem import Function
import numpy as np
import warnings
from matplotlib import cm
import pyvista as pv

def get_scalar_grid(fun: Function, varname: str, real = True):
    """
    This function extracts the dofs of a scalar function for the plot using pyvista.

    Parameters
    ----------
    fun : Function 
        Function from which the dofs are extracted.
    varname : str
        Name of the variable.
    real : boolean, optional (Default = True) 
        Real dofs are considered, if `False` imaginary dofs are used.
        
    Returns
    -------
    u_grid : pyvista.UnstructuredGrid
        Unstructured grid of values for the selected function.
    """
    topology, cells, geometry = dolfinx.plot.create_vtk_mesh(fun.function_space)
    u_grid = pv.UnstructuredGrid(topology, cells, geometry)

    if real:
        u_grid.point_data[varname] = fun.x.array[:].real
    else: 
        u_grid.point_data[varname] = fun.x.array[:].imag

    return u_grid

def extract_cells(domain: dolfinx.mesh.Mesh, points: np.ndarray):
    """
    This function can be used to extract data along a line defined by the variables `points`, crossing the domain.
 
    Parameters
    ----------
    domain  : dolfinx.mesh.Mesh
        Domain to extract data from.
    points : np.ndarray 
        Points listing the line from which data are extracted.

    Returns
    -------
    xPlot : np.ndarray 
        Coordinate denoting the cell from which data are extracted.
    cells : list
        List of cells of the mesh.
    """
    bb_tree = geometry.BoundingBoxTree(domain, domain.topology.dim)
    cells = []
    points_on_proc = []
    cell_candidates = geometry.compute_collisions(bb_tree, points.T)
    colliding_cells = geometry.compute_colliding_cells(domain, cell_candidates, points.T)
    for i, point in enumerate(points.T):
        if len(colliding_cells.links(i))>0:
            points_on_proc.append(point)
            cells.append(colliding_cells.links(i)[0])
    xPlot = np.array(points_on_proc, dtype=np.float64)

    return xPlot, cells

def PlotScalar(fun: Function, filename: str = None, format: str = 'png', varname: str = None,
               clim = None, colormap = cm.jet, resolution = [1080, 720], show=False):
    """
    Python function to plot a scalar field.

    Parameters
    ----------
    fun : Function 
        Field to plot.
    varname : str
        Name of the variable.
    filename : str
        Name of the file to save.
    clim : optional (Default = None)
        Colorbar limit, if `None` the mininum and maximum of `fun` are computed
    colormap : optional (Default = jet)
        Colormap for the plot
    resolution : list, optional (Default = [1080, 720])
        Resolution of the image

    """

    plotter = pv.Plotter(off_screen=True, border=False, window_size=resolution)
    lab_fontsize = 20
    title_fontsize = 25
    zoom = 1.
    
    if varname is None:
        varname = 'f'
    
    u_grid = get_scalar_grid(fun, varname)
    u_grid.set_active_scalars(varname)

    dict_cb = dict(title = varname, width = 0.75,
                    title_font_size=title_fontsize,
                    label_font_size=lab_fontsize,
                    color='k',
                    position_x=0.125, position_y=0.875,
                    shadow=True) 
    
    if clim is None:
        clim = [min(fun.x.array.real) * 0.975, max(fun.x.array.real)* 1.025]
        
    plotter.add_mesh(u_grid, cmap = colormap, clim = clim, show_edges=False, scalar_bar_args=dict_cb)
    plotter.view_xy()
    # plotter.add_title(varname, font_size=25, color ='k')
    plotter.camera.zoom(zoom)

    plotter.set_background('white', top='white')

    ## Save figure
    if filename is not None:
        if format == 'pdf':
            plotter.save_graphic(filename+'.pdf')
        elif format == 'svg':
            plotter.save_graphic(filename+'.svg')
        elif format == 'png':
            plotter.screenshot(filename+'.png', transparent_background = True,  window_size=resolution)
        else:
            warnings.warn("Available output format are 'pdf', 'svg' and 'png'. Saved 'png' screenshot.")
            plotter.screenshot(filename+'.png', transparent_background = True,  window_size=resolution)
        
    if show:
        plotter.show()