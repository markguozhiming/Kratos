from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
#from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
#from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.CompressiblePotentialFlowApplication import *
import potential_flow_solver
import os
import loads_output
import shutil
#import compute_lift_process_new
import subprocess
from PyPDF2 import PdfFileReader, PdfFileMerger
from math import *
######################################################################################
######################################################################################
######################################################################################
Number_Of_Refinements = TBD
Number_Of_AOAS = TBD

Initial_AOA = TBD
AOA_Increment = TBD

Initial_Airfoil_MeshSize = TBD
Airfoil_Refinement_Factor = TBD

Initial_FarField_MeshSize = TBD
FarField_Refinement_Factor = TBD

work_dir = '/home/inigo/simulations/naca0012/07_salome/05_MeshRefinement/'
input_mdpa_path = work_dir + 'mdpas/'
output_gid_path = '/media/inigo/10740FB2740F9A1C/Outputs/03_MeshRefinement/'

cl_results_file_name = work_dir + 'plots/cl/data/cl/cl_results.dat'
with open(cl_results_file_name,'w') as cl_file:
    cl_file.flush()
cl_results_directory_name = work_dir + 'plots/cl/data/cl'

cd_results_file_name = work_dir + 'plots/cd/data/cd/cd_results.dat'
with open(cd_results_file_name, 'w') as cd_file:
    cd_file.flush()
cd_results_directory_name = work_dir + 'plots/cd/data/cd'

aoa_results_file_name = work_dir + 'plots/aoa/cl_aoa.dat'
cl_aoa_file = open(aoa_results_file_name,'w')
cl_aoa_file.flush()

loads_output.write_header_all_cases(work_dir)

merger_global = PdfFileMerger()

case = 0
AOA = Initial_AOA
for j in range(Number_Of_AOAS):
    Airfoil_MeshSize = Initial_Airfoil_MeshSize
    FarField_MeshSize = Initial_FarField_MeshSize

    merger = PdfFileMerger()

    mesh_refinement_file_name = work_dir + 'plots/results/mesh_refinement_AOA_' + str(AOA)
    cl_data_directory_name = 'data/cl_AOA_' + str(AOA)
    cd_data_directory_name = 'data/cd_AOA_' + str(AOA)
    loads_output.write_header(work_dir)

    cp_data_directory_start = work_dir + 'plots/cp/data/AOA_' + str(AOA)
    os.mkdir(cp_data_directory_start)

    cl_aoa_file = open(aoa_results_file_name,'a')
    cl_aoa_file.write('{0:15f}'.format(AOA))
    cl_aoa_file.flush()

    os.mkdir(output_gid_path + 'AOA_' + str(AOA))

    for i in range(Number_Of_Refinements):
        print("\n\tCase ", case, "\n")
        loads_output.write_case(case, AOA, FarField_MeshSize, Airfoil_MeshSize, work_dir)

        cp_results_file_name = work_dir + 'plots/cp/data/0_original/cp_results.dat'
        cp_file = open(cp_results_file_name,'w')
        cp_file.flush()
        cp_coordinates_file_name = work_dir + 'plots/cp/data/0_original/coordinates.dat'
        cp_coordinates = open(cp_coordinates_file_name,'w')
        cp_coordinates.flush()
        cp_results_directory_name = work_dir + 'plots/cp/data/0_original'
        cp_data_directory_name = cp_data_directory_start + '/Case_' + str(case) + '_AOA_' + str(
                AOA) + '_Far_Field_Mesh_Size_' + str(FarField_MeshSize) + '_Airfoil_Mesh_Size_' + str(Airfoil_MeshSize)

        ## Parse the ProjectParameters
        #parameter_file = open("ProjectParameters_compressibility.json",'r')
        parameter_file = open("ProjectParameters.json",'r')
        ProjectParameters = Parameters( parameter_file.read())

        ## Get echo level and parallel type
        verbosity = ProjectParameters["problem_data"]["echo_level"].GetInt()
        parallel_type = ProjectParameters["problem_data"]["parallel_type"].GetString()

        ## Fluid model part definition
        main_model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
        main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

        alpharad = AOA*pi/180.0
        chord_normal_X = sin(alpharad)
        chord_normal_Y = cos(alpharad)
        
        main_model_part.ProcessInfo.SetValue(Y1, chord_normal_X)
        main_model_part.ProcessInfo.SetValue(Y2, chord_normal_Y)

        ###TODO replace this "model" for real one once available
        Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : main_model_part}

        #Set Mesh input_filename
        #mdpa_file_name = "Meshes/naca0012Mesh" +str(case)
        mdpa_file_name = input_mdpa_path + 'naca0012_Case_' + str(case) + '_AOA_' + str(
                AOA) + '_Far_Field_Mesh_Size_' + str(FarField_MeshSize) + '_Airfoil_Mesh_Size_' + str(Airfoil_MeshSize)

        ProjectParameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(mdpa_file_name)

        ## Solver construction    
        solver = potential_flow_solver.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

        solver.AddVariables()

        ## Read the model - note that SetBufferSize is done here
        solver.ImportModelPart()

        ## Add AddDofs
        solver.AddDofs()

        
        ## Set output name
        #problem_name = output_gid_path + ProjectParameters["problem_data"]["problem_name"].GetString()+ "Mesh" + str(case)
        problem_name = output_gid_path + 'AOA_' + str(AOA) + '/' + ProjectParameters["problem_data"]["problem_name"].GetString()+ '_Case_' + str(case) + '_AOA_' + str(
                AOA) + '_Far_Field_Mesh_Size_' + str(FarField_MeshSize) + '_Airfoil_Mesh_Size_' + str(Airfoil_MeshSize)

        ## Initialize GiD  I/O
        if (parallel_type == "OpenMP"):
            from gid_output_process import GiDOutputProcess
            gid_output = GiDOutputProcess(solver.GetComputingModelPart(),
                                        problem_name ,
                                        ProjectParameters["output_configuration"])

        gid_output.ExecuteInitialize()

        ##TODO: replace MODEL for the Kratos one ASAP
        ## Get the list of the skin submodel parts in the object Model
        for i in range(ProjectParameters["solver_settings"]["skin_parts"].size()):
            skin_part_name = ProjectParameters["solver_settings"]["skin_parts"][i].GetString()
            Model.update({skin_part_name: main_model_part.GetSubModelPart(skin_part_name)})

        ## Get the list of the no-skin submodel parts in the object Model (results processes and no-skin conditions)
        for i in range(ProjectParameters["solver_settings"]["no_skin_parts"].size()):
            no_skin_part_name = ProjectParameters["solver_settings"]["no_skin_parts"][i].GetString()
            Model.update({no_skin_part_name: main_model_part.GetSubModelPart(no_skin_part_name)})

        ## Print model_part and properties
        if(verbosity > 1):
            print("")
            print(main_model_part)
            for properties in main_model_part.Properties:
                print(properties)

        ## Processes construction
        import process_factory
        # "list_of_processes" contains all the processes already constructed (boundary conditions, initial conditions and gravity)
        # Note 1: gravity is firstly constructed. Outlet process might need its information.
        # Note 2: conditions are constructed before BCs. Otherwise, they may overwrite the BCs information.
        #ProjectParameters["boundary_conditions_process_list"]["model_import_settings"]["input_filename"].SetString(mdpa_file_name)
        #ProjectParameters["boundary_conditions_process_list"][2]["Parameters"]["mesh_refinement_file_name"].SetString(mesh_refinement_file_name)
        list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["boundary_conditions_process_list"] )

        if(verbosity > 1):
            for process in list_of_processes:
                print(process)

        ## Processes initialization
        for process in list_of_processes:
            process.ExecuteInitialize()

        ## Solver initialization
        solver.Initialize()

        #TODO: think if there is a better way to do this
        fluid_model_part = solver.GetComputingModelPart()

        ## Stepping and time settings
        # Dt = ProjectParameters["problem_data"]["time_step"].GetDouble()
        start_time = ProjectParameters["problem_data"]["start_time"].GetDouble()
        end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()

        time = start_time
        step = 0
        out = 0.0

        gid_output.ExecuteBeforeSolutionLoop()

        for process in list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        ## Writing the full ProjectParameters file before solving
        if ((parallel_type == "OpenMP") or (KratosMPI.mpi.rank == 0)) and (verbosity > 0):
            f = open("ProjectParametersOutput.json", 'w')
            f.write(ProjectParameters.PrettyPrintJsonString())
            f.close()

        Dt = 0.01
        step += 1
        time = time + Dt
        main_model_part.CloneTimeStep(time)

        if (parallel_type == "OpenMP") or (KratosMPI.mpi.rank == 0):
            print("")
            print("STEP = ", step)
            print("TIME = ", time)

        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()

        gid_output.ExecuteInitializeSolutionStep()

        solver.Solve()

        for process in list_of_processes:
            process.ExecuteFinalizeSolutionStep()

        gid_output.ExecuteFinalizeSolutionStep()

        #TODO: decide if it shall be done only when output is processed or not
        for process in list_of_processes:
            process.ExecuteBeforeOutputStep()

        if gid_output.IsOutputStep():
            gid_output.PrintOutput()

        for process in list_of_processes:
            process.ExecuteAfterOutputStep()

        out = out + Dt

        gid_output.ExecuteFinalize()

        loads_output.write_cp_figures(cp_data_directory_name, AOA, case, Airfoil_MeshSize, FarField_MeshSize, work_dir)

        shutil.copytree(cp_results_directory_name, cp_data_directory_name)

        latex = subprocess.Popen(['pdflatex', '-interaction=batchmode',work_dir + 'plots/cp/cp.tex'])
        latex.communicate()
        #path = "/home/inigo/simulations/naca0012/07_salome/01_MeshRefinement/"

        cp_file_name = work_dir + 'plots/cp/plots/cp_Case_' + str(case) + '_AOA_' + str(
                AOA) + '_Far_Field_Mesh_Size_' + str(FarField_MeshSize) + '_Airfoil_Mesh_Size_' + str(Airfoil_MeshSize) + '.pdf'
        shutil.copyfile('cp.pdf',cp_file_name)
        merger.append(PdfFileReader(cp_file_name), 'case_' + str(case))
        merger_global.append(PdfFileReader(cp_file_name), 'case_' + str(case))

        Airfoil_MeshSize /= Airfoil_Refinement_Factor
        FarField_MeshSize /= FarField_Refinement_Factor
        
        case +=1
    
    for process in list_of_processes:
        process.ExecuteFinalize()
    
    os.rename(work_dir + "mesh_refinement_loads.dat", mesh_refinement_file_name)
    
    loads_output.write_figures_cl(cl_data_directory_name, AOA, work_dir)
    loads_output.write_figures_cd(cd_data_directory_name, AOA, work_dir)
    
    shutil.copytree(cl_results_directory_name, work_dir + 'plots/cl/' + cl_data_directory_name)
    os.remove(cl_results_file_name)

    shutil.copytree(cd_results_directory_name, work_dir + 'plots/cd/' + cd_data_directory_name)
    os.remove(cd_results_file_name)


    cp_final_file_name = work_dir + 'plots/cp/cp_AOA_' + str(AOA) + '.pdf'
    merger.write(cp_final_file_name)
    AOA += AOA_Increment

cp_final_global_file_name = work_dir + 'plots/cp/cp_all.pdf'
merger_global.write(cp_final_global_file_name)
