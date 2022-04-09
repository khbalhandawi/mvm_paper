# -*- coding: utf-8 -*-
"""
Created on Tue Aug 08 02:41:33 2017

@author: Lydia

"""
import vtk
from vtk.util import numpy_support as VN
import numpy as np
from numpy.linalg import norm
import os

def read_unstructured_grid(filepath, subcase):
    """
    Reads unstructured grid from a VTK file
    :param filepath: Path to *.vtk file, containing unstructured grid
    :return: instanse of vtkUnstructuredGridReader, containing pre-processed grid
    """
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(filepath)
    reader.ReadAllScalarsOn()
    reader.ReadAllFieldsOn()
    reader.ReadAllTCoordsOn()
    reader.ReadAllVectorsOn()
    
    field = reader.GetFieldDataNameInFile(subcase)
    reader.SetFieldDataName(field)
    reader.Update()
    
    reader.Update()
    return reader

def color_unstructured_grid(grid, r, g, b):
    """
    Colors all the cells of input unstructured grid with a desired color
    :param grid: input grid, instance of vtkUnstructuredGrid
    :param r: R part of RGB
    :param g: G part of RGB
    :param b: B part of RGB
    :return:
    """
    colors = vtk.vtkUnsignedCharArray()
    colors.SetName("Colors")
    colors.SetNumberOfComponents(3)
    colors.SetNumberOfTuples(grid.GetNumberOfCells())
 
    for i in range(0, grid.GetNumberOfCells()):
        colors.InsertTuple3(i, r, g, b)
 
    grid.GetCellData().SetScalars(colors)
    
    return grid

def yield_background_grid_actor(grid):
    grid = color_unstructured_grid(grid, 224, 224, 224)
 
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(grid)
    mapper.SetScalarModeToUseCellData()
    
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetOpacity(1)
    actor.GetProperty().LightingOff()
    actor.GetProperty().EdgeVisibilityOn()
    actor.GetProperty().SetEdgeColor(0.36, 0.36, 0.36)
    actor.GetProperty().SetLineWidth(15)
    return actor

def set_view(az,el):
    #Set viewpoint specifications
        camera = vtk.vtkCamera()
        camera.ParallelProjectionOn()
        camera.SetParallelScale(829.801)
        camera.SetViewUp(0, 1, -1)
        camera.SetPosition(1081.75, 3206.44, 0.130707)
        camera.SetFocalPoint(1081.75, 0.339233, 0.130707)
        # set a 3D view according to az and el
        r = camera.GetRoll()

        camera.Azimuth(az)
        camera.Elevation(el)
        return camera

def make_marker(interactor):
    #create a X,Y,Z axes to show 3d position:
    # create axes variable and load vtk axes actor
    axes = vtk.vtkAxesActor()
    marker = vtk.vtkOrientationMarkerWidget()
    
    # set the interactor. self.widget._Iren is inbuilt python mechanism for current renderer window.
    marker.SetInteractor(interactor)
    marker.SetOrientationMarker(axes)
    
    # set size and position of window (Xmin,Ymin,Xmax,Ymax)
    marker.SetViewport(0.75,0,1,0.25)
    
    #Allow user input
    marker.SetEnabled(1)

def renderer_win(target, actor, camera, resolution, scalar_bar):
        # Setup renderer and rendererwindow
        renderer = vtk.vtkRenderer()
        renderer.AddActor(actor)
        renderer.AddActor(scalar_bar)

        bkc = np.divide([169, 160, 191],255, dtype=float)
        renderer.SetBackground(bkc[0], bkc[1], bkc[2])
        
        renderer.SetActiveCamera(camera)
        renderer.ResetCamera()
        
        renderer_window = vtk.vtkRenderWindow()
        renderer_window.AddRenderer(renderer)
        renderer_window.SetSize(resolution) # Set window size
        renderer_window.Render()
        return renderer_window

def interactor_win(renderer_window):
        # Setup renderer and rendererwindow
        interactor = vtk.vtkRenderWindowInteractor()
        interactor.SetRenderWindow(renderer_window)
    
        interactor_style = vtk.vtkInteractorStyleTrackballCamera()  # ParaView-like interaction
        interactor.SetInteractorStyle(interactor_style)
        return interactor
    
def close_window(interactor):
    render_window = interactor.GetRenderWindow()
    render_window.Finalize()
    interactor.TerminateApp()

def scalarbar_act(lookuptable, sc_title):
        # create the scalar_bar
        scalar_bar = vtk.vtkScalarBarActor()
        scalar_bar.SetOrientationToHorizontal()
        scalar_bar.SetTitle(sc_title)
        
        title_property = scalar_bar.GetTitleTextProperty()
        title_property.SetFontSize(20)
        scalar_bar.SetTitleTextProperty(title_property)
        scalar_bar.SetTitleRatio(1)

        scalar_bar.SetLookupTable(lookuptable)
        
        return scalar_bar

def scalarbar_wid(scalar_bar, interactor):         
        # create the scalar_bar_widget
        scalar_bar_widget = vtk.vtkScalarBarWidget()
        scalar_bar_widget.SetInteractor(interactor)
        scalar_bar_widget.SetScalarBarActor(scalar_bar)
        scalar_bar_widget.On()

def write_outs(renderer_window, target, pngname, textname):        
        window_to_image_filter = vtk.vtkWindowToImageFilter()
        window_to_image_filter.SetInput(renderer_window)
        window_to_image_filter.SetScale(1)
        window_to_image_filter.Update()
         
        writer = vtk.vtkPNGWriter()
        writer.SetFileName(pngname)
        writer.SetInputConnection(window_to_image_filter.GetOutputPort())
        writer.Write()
        
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(textname);
        if vtk.VTK_MAJOR_VERSION <= 5:
            writer.SetInput(target)
        else:
            writer.SetInputData(target)
        writer.Write()

def evaluate_mapper_on_defogrid(grid, isubcase, method): 
    
    point_cond = 'POINT'
    cell_cond = 'CELL'
    itime = 0
    # Join relevant strings
    
    name = 'T1'
    disp_comp_x = ''.join((name, '_isubcase_', str(isubcase), 
                            '_itime_', str(itime))) # Displacement component name
    name = 'T2'
    disp_comp_y = ''.join((name, '_isubcase_', str(isubcase), 
                            '_itime_', str(itime))) # Displacement component name
    name = 'T3'
    disp_comp_z = ''.join((name, '_isubcase_', str(isubcase), 
                            '_itime_', str(itime))) # Displacement component name
    
    # Get displacement Components
    nnodes = grid.GetNumberOfPoints() # Number of nodes
    coor = np.zeros([nnodes,3]) # Initialize coordinate array
    for i in range(nnodes):
        coor[i,:] = grid.GetPoints().GetPoint(i) # Get point coordinates
        
    T1 = VN.vtk_to_numpy(grid.GetPointData().GetScalars(disp_comp_x)) # Get T1
    T2 = VN.vtk_to_numpy(grid.GetPointData().GetScalars(disp_comp_y)) # Get T2
    T3 = VN.vtk_to_numpy(grid.GetPointData().GetScalars(disp_comp_z)) # Get T3
    T_comp = np.vstack((T1, T2, T3))
    
    req_d = 50 # Desired Deformation Scale Factor
    factor = req_d/max(norm(T_comp[:,:],axis=0)) # As a % of max deformation
    
    coor_def = np.zeros(np.shape(coor)) # Deformed coordinates array
    for k in range(np.shape(coor)[1]):
        coor_def[:,k] = coor[:,k] + T_comp[k,:]*factor # Compute deformed coordinates

    targetPoints = vtk.vtkPoints() # vtk deformed points object
        #mag = np.zeros(nnodes)
    for i in range(nnodes):
        point = coor_def[i,:]
        id = targetPoints.InsertNextPoint(point) # Create deformed points vtk array
        
    if method == point_cond: # Get nodal data
        #============================= MAP SCALAR DATA ===============================#
        targetpointdata = vtk.vtkFloatArray() # point data object
        targetpointdata.SetName("Field_magnitude")
        targetpointdata.SetNumberOfComponents(1)
        targetpointdata.SetReferenceCount(2)

        for i in range(nnodes):
            mag = norm(T_comp[:,i])
            id = targetpointdata.InsertNextValue(mag) # Add deformation mag to data object
            
        target = vtk.vtkUnstructuredGrid() # Define new grid
        target.SetPoints(targetPoints) # Add deformed points to new grid
        target.SetCells(grid.GetCellType(1),grid.GetCells()) # grid.GetCellType(1) is common for all cells

        target.GetPointData().SetScalars(targetpointdata) # Add desired scalers to deformed grid

        #============================= MAP SCALAR DATA ===============================#

        scalar_range_p = targetpointdata.GetValueRange() # Set colormap range to mapped data range        
        
        mapper = vtk.vtkDataSetMapper()
        if vtk.VTK_MAJOR_VERSION <= 5:
            mapper.SetInput(target)
        else:
            mapper.SetInputData(target) # Provide deformed grid
        
        mapper.SelectColorArray(0) # Select Array to map
        mapper.SetColorModeToDefault()
        mapper.SetScalarModeToUsePointData() # Use point data as mapping method
        
        mapper.SetColorModeToMapScalars()
        mapper.SetScalarRange(scalar_range_p) # scalar_range
        mapper.Update()
        
        lookuptable = mapper.GetLookupTable() # Retrieve color lookuptable
        lookuptable.SetHueRange((2./3),0) # Reverse color Bar
    
        max_disp = max(norm(T_comp[:,:],axis=0))
    
        return max_disp, lookuptable, mapper, target

    #-----------------------------------------------------------------------------#
    #READ STRESS RESULTS
    elif method == cell_cond: # Cell data method

        vmstr = VN.vtk_to_numpy(grid.GetCellData().GetArray(isubcase-1)) # Get von mises stress
        
        targetcelldata = vtk.vtkFloatArray() # Cell data object
        targetcelldata.SetName("von_mises")
        targetcelldata.SetNumberOfComponents(1)
        targetcelldata.SetReferenceCount(2)
        
        ncells = grid.GetNumberOfCells()
        for i in range(ncells):
            vm = vmstr[i]
            id = targetcelldata.InsertNextValue(vm) # Update vonmises values to data object

        target = vtk.vtkUnstructuredGrid() # Define new grid
        target.SetPoints(targetPoints) # Add deformed points to new grid
        target.SetCells(grid.GetCellType(1),grid.GetCells()) # grid.GetCellType(1) is common for all cells

        data = target.GetCellData()
        data.AddArray(targetcelldata) # Include cell data in target grid

        scalar_range_c = targetcelldata.GetValueRange() # Set colormap range to mapped data range

        mapper = vtk.vtkDataSetMapper()
        if vtk.VTK_MAJOR_VERSION <= 5:
            mapper.SetInput(target)
        else:
            mapper.SetInputData(target) # Provide deformed grid
        
        mapper.SelectColorArray(0) # Select Array to map
        mapper.SetColorModeToDefault()
        mapper.SetScalarModeToUseCellFieldData() # Use cell field data as mapping method

        mapper.SetColorModeToMapScalars()
        mapper.Update()
        mapper.SetScalarRange(scalar_range_c) # scalar_range
        mapper.Update()
        
        lookuptable = mapper.GetLookupTable()
        lookuptable.SetHueRange((2./3),0) # Reverse color Bar
        
        max_stress = max(vmstr)
        return max_stress, lookuptable, mapper, target

#=============================================================================#
def plotresults(isubcase,filepath,method,output_path):
    #============================= USER-DEFINITION ===============================#
    az = 0
    el = 90  
    resolution = (1200, 1200)
    pngname = ''.join(("plot", "_subcase", str(isubcase), '_', method, ".png")) 
    textname = "target.vtk"
    
    if method == 'POINT':
        sc_title = 'Displacement (mm)'
    elif method == 'CELL':
        sc_title = 'S11 (MPa)'
    
    test_path = os.path.dirname(os.path.realpath(__file__)) # Working directory of file
    folder = 'Nastran_output'
    folder_path = os.path.join(test_path, folder)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    
    out_filename_png = os.path.join(output_path, pngname) # output png file name
    out_filename_vtk = os.path.join(output_path, textname) # output vtk file name
    #============================ create deformed grid ===========================#
    reader = read_unstructured_grid(filepath, 0) # Retrieve vtk grid
    grid = reader.GetOutput()

    # retrieve mapper and target grid for desired FEA result type
    [max_v, lookuptable, mapper, target] = evaluate_mapper_on_defogrid(grid, isubcase, method)
    #================================ ADD ACTORS =================================#

    actor = vtk.vtkActor()
    actor.SetMapper(mapper) 
    actor.GetProperty().SetLineWidth(15)
    
    camera = set_view(az,el) # Set camera angles
    scalar_bar = scalarbar_act(lookuptable, sc_title) 
    renderer_window = renderer_win(target, actor, camera, resolution, scalar_bar)
    interactor = interactor_win(renderer_window)
    make_marker(interactor)
    scalarbar_wid(scalar_bar, interactor)

    #-------------------------------------------#
    interactor.Initialize()
    write_outs(renderer_window, target, out_filename_png, out_filename_vtk)
    
    # close_window(interactor)
    # del renderer_window, interactor
    
    return max_v

def plotmesh(filepath,output_path):
    #============================= USER-DEFINITION ===============================#
    az = 0
    el = 90
    resolution = (1200, 1200)
    pngname = ''.join(("plot", "_mesh", ".png")) 
    textname = "target.vtk"

    out_filename_png = os.path.join(output_path, pngname) # output png file name
    out_filename_vtk = os.path.join(output_path, textname) # output vtk file name
    #============================ create deformed grid ===========================#
    reader = read_unstructured_grid(filepath, 0) # Retrieve vtk grid
    grid = reader.GetOutput()

    # retrieve mapper and target grid for desired FEA result type
    #================================ ADD ACTORS =================================#

    camera = set_view(az,el) # Set camera angles
    
    renderer = vtk.vtkRenderer()
    renderer.AddActor(yield_background_grid_actor(grid))
    renderer.SetBackground(1, 1, 1)
    renderer.SetActiveCamera(camera);
    renderer.ResetCamera()
    
    renderer_window = vtk.vtkRenderWindow()
    renderer_window.AddRenderer(renderer)
    renderer_window.SetSize(resolution) # Set window size
    renderer_window.Render()
    
    interactor = interactor_win(renderer_window)
    #-------------------------------------------#
    interactor.Initialize()
    write_outs(renderer_window, grid, out_filename_png, out_filename_vtk)
    
    # close_window(interactor)
    # del renderer_window, interactor

    pass

def postprocess_vtk(isubcase,filepath,method):
    reader = read_unstructured_grid(filepath, 0) # Retrieve vtk grid
    grid = reader.GetOutput()
    
    point_cond = 'POINT'
    cell_cond = 'CELL'
    itime = 0
    
    if method == point_cond: # Get nodal data
        # Join relevant strings
        name = 'T1'
        disp_comp_x = ''.join((name, '_isubcase_', str(isubcase), 
                               '_itime_', str(itime))) # Displacement component name
        name = 'T2'
        disp_comp_y = ''.join((name, '_isubcase_', str(isubcase), 
                               '_itime_', str(itime))) # Displacement component name
        name = 'T3'
        disp_comp_z = ''.join((name, '_isubcase_', str(isubcase), 
                               '_itime_', str(itime))) # Displacement component name
        # Get displacement Components
        
        T1 = VN.vtk_to_numpy(grid.GetPointData().GetScalars(disp_comp_x)) # Get T1
        T2 = VN.vtk_to_numpy(grid.GetPointData().GetScalars(disp_comp_y)) # Get T2
        T3 = VN.vtk_to_numpy(grid.GetPointData().GetScalars(disp_comp_z)) # Get T3
        T_comp = np.vstack((T1, T2, T3))

        # Get nodal coordinates
        nnodes = grid.GetNumberOfPoints() # Number of nodes
        coor = np.zeros([nnodes,3]) # Initialize coordinate array
        for i in range(nnodes):
            coor[i,:] = grid.GetPoints().GetPoint(i) # Get point coordinates
            if grid.GetPoints().GetPoint(i) == tuple(3*[0.0]):
                i_center = i

        #---------------------------------------------------------------------#
        max_out = max(norm(T_comp[:,:],axis=0)) # Maximum displacement result
        center_U = norm(T_comp[:,:],axis=0)[i_center]
        return center_U
    #=========================== MAP SCALAR DATA =============================#
    elif method == cell_cond: # Cell data method
        #READ STRESS RESULTS
        vmstr = VN.vtk_to_numpy(grid.GetCellData().GetArray(isubcase-1)) # Get von mises stress
        #---------------------------------------------------------------------#
        max_out = max(abs(vmstr))  # Maximum von mises stress result
        
        return max_out

if __name__ == '__main__':

    base_path = os.getcwd() # Working directory
    folder = 'Nastran_output'
    output_path = os.path.join(base_path, 'examples', 'CAD', folder)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    base_sim_name = 'trs_sim1-bearing_load'

    path = os.getcwd()
    vtk_file = os.path.join(path,'examples','CAD',base_sim_name+'.vtk')

    output = postprocess_vtk(1,vtk_file,'POINT')

    print(output)

    plotmesh(vtk_file,output_path)

    plotresults(1,vtk_file,'CELL',output_path)
    plotresults(1,vtk_file,'POINT',output_path)

    # base_sim_name = 'trs_sim1-bearing_load_fixed_flange'

    # path = os.getcwd()
    # vtk_file = os.path.join(path,'examples','CAD',base_sim_name+'.vtk')

    # output = postprocess_vtk(1,vtk_file,'POINT')

    # print(output)

    # plotmesh(vtk_file,output_path)

    # plotresults(1,vtk_file,'CELL',output_path)
    # plotresults(1,vtk_file,'POINT',output_path)