# NX 1867
# Journal created by Khalil on Fri Feb 25 10:47:01 2022 W. Europe Standard Time
#
import math
import NXOpen
import NXOpen.CAE
import NXOpen.Preferences
import os, sys

import time

# from pyNastran.op2.export_to_vtk import export_to_vtk, export_to_vtk_filename

def list_expressions(workPart,part_design):

    objects = []

    for name, value in part_design.items():
        if not value['skip']:

            if value['unit'] == 'no_unit':
                unit = NXOpen.Unit.Null
            else:
                unit = workPart.UnitCollection.FindObject(value['unit'])

            expression = workPart.Expressions.FindObject(name)
            workPart.Expressions.EditWithUnits(expression, unit, str(value['value']))

            objects += [expression]

    return objects

def update_expressions(theSession,part_name,part_design):

	theSession.Preferences.Modeling.UpdatePending = False
	workPart = theSession.Parts.FindObject(part_name)
	objects = list_expressions(workPart,part_design)

	# UPDATE CHANGES
	[status1, partLoadStatus1] = theSession.Parts.SetDisplay(workPart, False, True)

	markId1 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Make Up to Date")
	theSession.UpdateManager.MakeUpToDate(objects, markId1)
	partLoadStatus1.Dispose()

def post_processing(base_sim_name):

    from cpylog import get_logger
    from .export_to_vtk import export_to_vtk_filename
    logger = get_logger(level='error')

    path = os.getcwd()
    op2_filename = os.path.join(path,'CAD',base_sim_name+'.op2')
    bdf_filename = os.path.join(path,'CAD',base_sim_name+'.dat')
    vtk_filename = os.path.join(path,'CAD',base_sim_name+'.vtk')

    export_to_vtk_filename(bdf_filename,op2_filename,vtk_filename,log=logger)

    return vtk_filename

def run_nx_simulation(vane_length,vane_height,lean_angle,n_struts,
    r_hub=346.5,r_shroud=536.5,yield_strength=460,youngs_modulus=156.3e3,
    poissons_ratio=0.33,material_density=8.19e-06,bearing_x=600,bearing_y=500,
    mesh_size=10,pid=None):

    old_stdout = sys.stdout
    devnull = open(os.devnull, 'w')
    sys.stdout = devnull

    theSession  = NXOpen.Session.GetSession()
    # ----------------------------------------------
    #   Menu: File->Open...
    # ----------------------------------------------
    path = os.getcwd()
    part_name = 'trs'
    part_path = os.path.join(path,'CAD',part_name+'.prt')

    theSession.Parts.OpenBaseDisplay(part_path)
    workPart = theSession.Parts.Work
    displayPart = theSession.Parts.Display
    theSession.ApplicationSwitchImmediate("UG_APP_MODELING")

	# ----------------------------------------------
	# Update journal part (controls TRS design)
    geometry_params = {
        'r_hub' 		        : {'value' : r_hub	        , 'skip' : False, 'unit' : 'MilliMeter'	},
        'r_shroud' 			    : {'value' : r_shroud		, 'skip' : False, 'unit' : 'MilliMeter'	},
        'lean_angle'            : {'value' : lean_angle	    , 'skip' : False, 'unit' : 'Degrees'	},
        'n_struts' 				: {'value' : n_struts		, 'skip' : False, 'unit' : 'no_unit'	},
    }

    update_expressions(theSession,part_name,geometry_params)

    # ----------------------------------------------
    #   Menu: File->Open...
    # ----------------------------------------------
    part_name = 'trs_fem1'
    part_path = os.path.join(path,'CAD',part_name+'.fem')

    theSession.Parts.OpenBaseDisplay(part_path)
    workFemPart = theSession.Parts.BaseWork
    displayFemPart = theSession.Parts.BaseDisplay
    theSession.ApplicationSwitchImmediate("UG_APP_SFEM")

    mesh_params = {
        'yield_strength'        : {'value' : yield_strength	    , 'skip' : False, 'unit' : 'StressNewtonPerSquareMilliMeter'},
        'youngs_modulus'        : {'value' : youngs_modulus	    , 'skip' : False, 'unit' : 'StressNewtonPerSquareMilliMeter'},
        'poissons_ratio'        : {'value' : poissons_ratio	    , 'skip' : False, 'unit' : 'no_unit'},
        'material_density'      : {'value' : material_density   , 'skip' : False, 'unit' : 'KilogramPerCubicMilliMeter'},
        'vane_height'           : {'value' : vane_height	    , 'skip' : False, 'unit' : 'MilliMeter'},
        'vane_width'            : {'value' : vane_length        , 'skip' : False, 'unit' : 'MilliMeter'},
    }

    update_expressions(theSession,part_name,mesh_params)

    # ----------------------------------------------
    #   Menu: File->Utilities->Activate Meshing
    # ----------------------------------------------
    femPart1 = theSession.Parts.FindObject("trs_fem1")
    theSession.Parts.SetWork(femPart1)

    workFemPart = theSession.Parts.BaseWork
    origin = NXOpen.Point3d(0.0, 0.0, 0.0)
    # ----------------------------------------------
    #   Menu: Insert->Mesh->1D Mesh...
    # ----------------------------------------------

    fEModel1 = workFemPart.FindObject("FEModel")
    meshManager1 = fEModel1.Find("MeshManager")
    mesh1dBuilder1 = meshManager1.CreateMesh1dBuilder(NXOpen.CAE.Mesh1d.Null)

    meshCollector1 = meshManager1.FindObject("MeshCollector[trs_beams]")
    mesh1dBuilder1.ElementType.DestinationCollector.ElementContainer = meshCollector1

    mesh1dBuilder1.FlipDirectionOption = False

    mesh1dBuilder1.ElementType.ElementTypeName = "CBEAM"

    mesh1dBuilder1.ElementType.DestinationCollector.AutomaticMode = False

    mesh1dBuilder1.PropertyTable.SetIntegerPropertyValue("mesh density option", 0)

    unit2 = workFemPart.UnitCollection.FindObject("MilliMeter")
    mesh1dBuilder1.PropertyTable.SetBaseScalarWithDataPropertyValue("mesh density", str(mesh_size), unit2)

    mesh1dBuilder1.PropertyTable.SetBaseScalarWithDataPropertyValue("mesh density number", "10", NXOpen.Unit.Null)

    mesh1dBuilder1.PropertyTable.SetBooleanPropertyValue("merge nodes bool", True)

    boundingVolumeSelectionRecipe1 = workFemPart.SelectionRecipes.FindObject("trs")

    mesh1dBuilder1.SelectionList.Add(boundingVolumeSelectionRecipe1, workFemPart.ModelingViews.WorkView, origin)

    mesh1dBuilder1.ElementType.ElementDimension = NXOpen.CAE.ElementTypeBuilder.ElementType.Beam

    mesh1dBuilder1.ElementType.ElementTypeName = "CBEAM"

    destinationCollectorBuilder1 = mesh1dBuilder1.ElementType.DestinationCollector

    destinationCollectorBuilder1.ElementContainer = meshCollector1

    destinationCollectorBuilder1.AutomaticMode = False

    mesh1dBuilder1.PropertyTable.SetBooleanPropertyValue("mesh coat bool", False)

    mesh1dBuilder1.PropertyTable.SetIntegerPropertyValue("mesh density option", 0)

    mesh1dBuilder1.PropertyTable.SetBaseScalarWithDataPropertyValue("mesh density", str(mesh_size), unit2)

    mesh1dBuilder1.PropertyTable.SetBaseScalarWithDataPropertyValue("mesh density number", "10", NXOpen.Unit.Null)

    mesh1dBuilder1.PropertyTable.SetBooleanPropertyValue("merge nodes bool", True)

    mesh1dBuilder1.PropertyTable.SetBaseScalarWithDataPropertyValue("node coincidence tolerance", "0.0001", unit2)

    mesh1dBuilder1.PropertyTable.SetIntegerPropertyValue("mesh edit allowed", 0)

    mesh1dBuilder1.PropertyTable.SetBooleanPropertyValue("use mid nodes bool", False)

    mesh1dBuilder1.PropertyTable.SetIntegerPropertyValue("mesh time stamp", 0)

    mesh1dBuilder1.PropertyTable.SetIntegerPropertyValue("2d phantom mesh", 0)

    mesh1dBuilder1.CommitMesh()

    mesh1dBuilder1.Destroy()

    # ----------------------------------------------
    #   Menu: Insert->Mesh->1D Connection...
    # ----------------------------------------------
    cAEConnectionBuilder1 = fEModel1.CaeConnections.CreateConnectionBuilder(NXOpen.CAE.CAEConnection.Null)

    cAEConnectionBuilder1.ElementType.DestinationCollector.ElementContainer = meshCollector1

    cAEConnectionBuilder1.ElementTypeRbe3.DestinationCollector.ElementContainer = meshCollector1

    cAEConnectionBuilder1.ElementType.ElementDimension = NXOpen.CAE.ElementTypeBuilder.ElementType.Connection

    cAEConnectionBuilder1.ElementTypeRbe3.ElementDimension = NXOpen.CAE.ElementTypeBuilder.ElementType.Spider

    meshCollector2 = meshManager1.FindObject("MeshCollector[rigid]")
    cAEConnectionBuilder1.ElementTypeRbe3.DestinationCollector.ElementContainer = meshCollector2

    cAEConnectionBuilder1.MidNode = True

    cAEConnectionBuilder1.Type = NXOpen.CAE.CAEConnectionBuilder.ConnectionTypeEnum.PointToEdge

    cAEConnectionBuilder1.ElementType.ElementTypeName = "RBE2"

    cAEConnectionBuilder1.ElementType.DestinationCollector.ElementContainer = meshCollector2

    cAEConnectionBuilder1.ElementType.DestinationCollector.AutomaticMode = False

    cAEConnectionBuilder1.ElementTypeRbe3.ElementTypeName = "RBE2"

    cAEConnectionBuilder1.Label = 241

    cAEConnectionBuilder1.ElementType.ElementTypeName = "RBE2"

    boundingVolumeSelectionRecipe2 = workFemPart.SelectionRecipes.FindObject("center_point")
    cAEConnectionBuilder1.SourceSelection.Add(boundingVolumeSelectionRecipe2, workFemPart.ModelingViews.WorkView, origin)

    boundingVolumeSelectionRecipe3 = workFemPart.SelectionRecipes.FindObject("inner_casing")
    cAEConnectionBuilder1.TargetSelection.Add(boundingVolumeSelectionRecipe3, workFemPart.ModelingViews.WorkView, origin)

    cAEConnectionBuilder1.ElementType.ElementTypeName = "RBE2"

    cAEConnectionBuilder1.NodeFaceProximity = 0.0

    cAEConnectionBuilder1.Commit()

    cAEConnectionBuilder1.Destroy()

    # ----------------------------------------------
    #   Menu: File->Open...
    # ----------------------------------------------
    part_name = 'trs_sim1'
    part_path = os.path.join(path,'CAD',part_name+'.sim')

    theSession.Parts.OpenBaseDisplay(part_path)
    workSimPart = theSession.Parts.BaseWork
    displaySimPart = theSession.Parts.BaseDisplay

    theSession.Post.UpdateUserGroupsFromSimPart(workSimPart)

    sim_params = {
        'bearing_x'             : {'value' : bearing_x*1e4, 'skip' : False, 'unit' : 'Newton'},
        'bearing_y'             : {'value' : bearing_y*1e4, 'skip' : False, 'unit' : 'Newton'},
    }

    update_expressions(theSession,part_name,sim_params)

    # ----------------------------------------------
    #   Menu: File->Utilities->Activate Simulation
    # ----------------------------------------------
    theSession.Parts.SetWork(displaySimPart)

    theSimSolveManager = NXOpen.CAE.SimSolveManager.GetSimSolveManager(theSession)

    suffix = ''
    if pid is not None:
        suffix += '_' + str(pid)

    base_sim_name_1 = 'bearing_load'+suffix
    base_sim_name_2 = 'bearing_load_fixed_flange'+suffix
    base_sim_name_1 = base_sim_name_1.lower()
    base_sim_name_2 = base_sim_name_2.lower()
    results_file_1 = 'trs_sim1-'+base_sim_name_1
    results_file_2 = 'trs_sim1-'+base_sim_name_2
    results_path_1 = os.path.join(path,'CAD',results_file_1+'.op2')
    results_path_2 = os.path.join(path,'CAD',results_file_2+'.op2')

    ######################################################
    #                Run both load cases                 #
    ######################################################

    # psolutions1 = [NXOpen.CAE.SimSolution.Null] * 2
    # simSimulation1 = workSimPart.FindObject("Simulation")
    
    # simSolution1 = simSimulation1.FindObject("Solution[bearing_load]")
    # simSolution1.Rename(base_sim_name_1, False)
    # theSession.Post.UnloadResultFile(results_path_1)
    # psolutions1[0] = simSolution1

    # simSolution2 = simSimulation1.FindObject("Solution[bearing_load_fixed_flange]")
    # simSolutionStep1 = simSolution2.FindObject("SolutionStep[Subcase - Static Loads 1]")
    # simSolution2.ActiveStep = simSolutionStep1
    # simSolution2.Rename(base_sim_name_2, False)
    # theSession.Post.UnloadResultFile(results_path_2)
    # psolutions1[1] = simSolution2

    ######################################################
    #            Run only the fixed flange case          #
    ######################################################

    psolutions1 = [NXOpen.CAE.SimSolution.Null] * 1
    simSimulation1 = workSimPart.FindObject("Simulation")

    simSolution2 = simSimulation1.FindObject("Solution[bearing_load_fixed_flange]")
    simSolutionStep1 = simSolution2.FindObject("SolutionStep[Subcase - Static Loads 1]")
    simSolution2.ActiveStep = simSolutionStep1
    simSolution2.Rename(base_sim_name_2, False)
    theSession.Post.UnloadResultFile(results_path_2)
    psolutions1[0] = simSolution2

    ######################################################
    #                    Execute jobs                    #
    ######################################################
    run_successful = False

    try:
        numsolutionssolved1, numsolutionsfailed1, numsolutionsskipped1 = theSimSolveManager.SolveChainOfSolutions(psolutions1, NXOpen.CAE.SimSolution.SolveOption.Solve, NXOpen.CAE.SimSolution.SetupCheckOption.CompleteCheckAndOutputErrors, NXOpen.CAE.SimSolution.SolveMode.Foreground)
        time.sleep(1.0)
        nx_exception = False
    except NXOpen.NXException:
        nx_exception = True

    if os.path.exists(results_path_2) and not nx_exception:
        run_successful = True

    # ----------------------------------------------
    #   Menu: File->Close->All Parts
    # ----------------------------------------------
    theSession.Parts.CloseAll(NXOpen.BasePart.CloseModified.CloseModified, None)

    workSimPart = NXOpen.BasePart.Null
    displaySimPart = NXOpen.BasePart.Null
    theSession.ApplicationSwitchImmediate("UG_APP_NOPART")
    
    # ----------------------------------------------
    #   POST-PROCESSING
    # ----------------------------------------------
    if run_successful:
        # vtk_filename_1 = post_processing(results_file_1)
        vtk_filename_2 = post_processing(results_file_2)
    else:
        # vtk_filename_1 = None
        vtk_filename_2 = None

    return None, vtk_filename_2, run_successful

if __name__ == '__main__':

    n_struts            = 5
    vane_length         = 120
    vane_height         = 20
    lean_angle          = 20
    bearing_x           = 50
    bearing_y           = 50

    n_struts            = 9.0
    vane_length         = 23.384918013185924
    vane_height         = 7.11823399887898
    lean_angle          = 27.61749330865683
    bearing_x           = 176.54130348897905
    bearing_y           = 176.89984971500544

    vane_length = 80.23418632896156
    vane_height = 13.129238383885887
    lean_angle = 33.012277168282715
    n_struts = 19.0
    bearing_x = 115.32711545052709
    bearing_y = 113.00224669379733

    _, vtk_filename, success = run_nx_simulation(vane_length=vane_length,vane_height=vane_height,
        lean_angle=lean_angle,n_struts=n_struts,r_hub=346.5,r_shroud=536.5,yield_strength=460,
        youngs_modulus=156.3e3,poissons_ratio=0.33,material_density=8.19e-06,bearing_x=bearing_x,
        bearing_y=bearing_y,mesh_size=10,pid=None)

    base_path = os.getcwd() # Working directory
    folder = 'Nastran_output'
    output_path = os.path.join(base_path, 'CAD', folder)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    if success:
        from export_to_vtk import export_to_vtk_filename
        from plot_results import plotmesh, plotresults,postprocess_vtk

        plotmesh(vtk_filename,output_path)
        plotresults(1,vtk_filename,'POINT',output_path)
        plotresults(1,vtk_filename,'CELL',output_path)

        stress = postprocess_vtk(1,vtk_filename,'CELL')
        print(stress)