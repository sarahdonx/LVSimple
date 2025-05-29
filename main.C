// Copyright (c) 2011-2012, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// Config files
//#include <IBAMR_config.h>
#include <ibamr/config.h>
//#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petsc.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFECentroidPostProcessor.h>
#include <ibamr/IBFEMethod.h>

#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>
#include <ibtk/IBTKInit.h>

// Headers for application components.
#include <ActiveContraction.h>
#include <BoundaryConditions.h>
#include <MechanicsModel.h>
#include <ModelInitialization.h>
#include <tetVolumeCalculation.h>

//for linux system call
#include <stdlib.h>
#include <stdio.h>
#include <chrono>
#include <ctime>

// Forward declaration of the main driver routine that actually runs the
// simulation.
int
main_driver(
    int argc,
    char* argv[]);

// Main routine that initializes the simulation libraries and calls the main
// driver routine.
int
main(
    int argc,
    char* argv[])
{
    main_driver(argc, argv);
    return 0;
}// main

// Main driver routine that actually runs the simulation.
int
main_driver(
    int argc,
    char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    //LibMeshInit init(argc, argv);
    
    //SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    //SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    //SAMRAIManager::startup();

    // Increase maximum patch data component indices
    //SAMRAIManager::setMaxNumberPatchDataEntries(2500);

 { // cleanup dynamically allocated objects prior to shutdown


    // Parse command line options, set some standard options from the input
    // file, initialize the restart database (if this is a restarted run), and
    // enable file logging.
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "FELV.log");
    Pointer<Database> input_db = app_initializer->getInputDatabase();

    // Get various standard options set in the input file.
    const bool dump_viz_data = app_initializer->dumpVizData();
    const int viz_dump_interval = app_initializer->getVizDumpInterval();
    //const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();
    const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
    const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
    const string exodus_filename = app_initializer->getExodusIIFilename();

    const bool dump_restart_data = app_initializer->dumpRestartData();
    const int restart_dump_interval = app_initializer->getRestartDumpInterval();
    const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();
    const string restart_read_dirname = app_initializer->getRestartReadDirectory();
    const int restart_restore_num = app_initializer->getRestartRestoreNumber();

    const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
    const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
    const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
    if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
    {
        Utilities::recursiveMkdir(postproc_data_dump_dirname);
    }

    const bool dump_timer_data = app_initializer->dumpTimerData();
    const int timer_dump_interval = app_initializer->getTimerDumpInterval();

    // Configure the passive stress model.
    MechanicsModel::normalize_stress = input_db->getBool("NORMALIZE_STRESS");

    // Configure the active tension model.
    MechanicsModel::enable_active_tension = input_db->getBool("ENABLE_ACTIVE_TENSION");
    MechanicsModel::T_scale = input_db->getDouble("T_SCALE");

    // Configure the pressure loading.
    BoundaryConditions::P_load = input_db->getDouble("P_LOAD");
    BoundaryConditions::P_load_es = input_db->getDouble("P_LOAD_ES");
    BoundaryConditions::t_load = input_db->getDouble("T_LOAD");

    // Configure the penalty parameters.
    BoundaryConditions::kappa = input_db->getDouble("KAPPA");
    MechanicsModel::beta_s = input_db->getDouble("BETA_S");

    // Load and initialize the FE mesh.
    Mesh mesh(init.comm(), NDIM);
    mesh.allow_renumbering(false);
    ModelInitialization::initialize_mesh(mesh, input_db);


    // Create major IBAMR solver objects.
    Pointer<IBFEMethod> ib_method_ops = new IBFEMethod("IBFEMethod",
                     app_initializer->getComponentDatabase("IBFEMethod"),
                     &mesh,
                     app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                    /*register_for_restart*/  true, restart_read_dirname, restart_restore_num);
    Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                     "INSStaggeredHierarchyIntegrator",
                     app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
    Pointer<IBHierarchyIntegrator> time_integrator = new IBExplicitHierarchyIntegrator(
                     "IBExplicitHierarchyIntegrator",
                     app_initializer->getComponentDatabase("IBExplicitHierarchyIntegrator"),
                     ib_method_ops,
                     navier_stokes_integrator);

    // Create major SAMRAI algorithm and data objects.
    Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
                      "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>(
                      "PatchHierarchy", grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM> > error_detector =
                      new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                      time_integrator,
                      app_initializer->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM> > load_balancer = new LoadBalancer<NDIM>(
                       "LoadBalancer",
                       app_initializer->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm = new GriddingAlgorithm<NDIM>(
                        "GriddingAlgorithm",
                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                        error_detector, box_generator, load_balancer);

    // Build equation systems that store various auxiliary variables, such as
    // the material axes and the values of spatially distributed active
    // contraction model state variables.

    // Configure the IBFE solver.
    pout << "Configuring the solver...\n";

    IBFEMethod::PK1StressFcnData PK1_dev_stress_data;
    //need to prepare the system data for using fibre direction
    std::vector<int> vars_fibre(3);
    vars_fibre[0] = 0;
    vars_fibre[1] = 1;
    vars_fibre[2] = 2;
    std::vector<int> vars_sheet(3);
    vars_sheet[0] = 0;
    vars_sheet[1] = 1;
    vars_sheet[2] = 2;
    std::vector<int> vars_acT(1);
    vars_acT[0] = 0;
    std::vector<SystemData> fsnAct_data(3);
    fsnAct_data[0] = SystemData("fiber direction", vars_fibre);
    fsnAct_data[1] = SystemData("sheet direction", vars_sheet);
    fsnAct_data[2] = SystemData("active tension", vars_acT);
    // the SystemData requires the same system name, and vars_fibre includes the index of variables

    PK1_dev_stress_data.fcn = MechanicsModel::PK1_dev_stress_function;
    //MechanicsModel::get_PK1_dev_stress_function_systems(PK1_dev_stress_data.systems);
    PK1_dev_stress_data.quad_type = QGAUSS;
    PK1_dev_stress_data.quad_order = THIRD;  // full-order integration
    PK1_dev_stress_data.system_data = fsnAct_data;
    ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data);

    IBFEMethod::PK1StressFcnData PK1_dil_stress_data;
    PK1_dil_stress_data.fcn = MechanicsModel::PK1_dil_stress_function;
	  //MechanicsModel::get_PK1_dev_stress_function_systems(PK1_dil_stress_data.systems);
    PK1_dil_stress_data.quad_type = QGAUSS;
    PK1_dil_stress_data.quad_order = FIRST;  // reduced-order integration
    PK1_dil_stress_data.system_data = fsnAct_data;
    ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data);

    std::vector<int> vars_dis(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) vars_dis[d] = d;
    vector<SystemData> sys_data_dis(1, SystemData(IBFEMethod::COORD_MAPPING_SYSTEM_NAME, vars_dis));
    IBFEMethod::LagSurfacePressureFcnData surface_pressure_data(
              BoundaryConditions::loading_force_function, sys_data_dis, NULL);
    ib_method_ops->registerLagSurfacePressureFunction(surface_pressure_data);
    IBFEMethod::LagSurfaceForceFcnData surface_fcn_data(
              BoundaryConditions::tether_force_function,sys_data_dis,NULL);
    ib_method_ops->registerLagSurfaceForceFunction(surface_fcn_data);

    ib_method_ops->initializeFEEquationSystems();
    EquationSystems* equation_systems =
                         ib_method_ops->getFEDataManager()->getEquationSystems();
    ModelInitialization::initialize_equation_systems(equation_systems);


    // Set up a post-processor to reconstruct the various quantities of interest.
    Pointer<IBFEPostProcessor> ib_post_processor =
    new IBFECentroidPostProcessor("IBFEPostProcessor", ib_method_ops->getFEDataManager());

    //1) Deformation gradient tensor FF = dX/ds.
    ib_post_processor->registerTensorVariable(
        "FF", MONOMIAL, CONSTANT, IBFEPostProcessor::FF_fcn);

    // 2) Deformed fiber and sheet axes.
    vector<SystemData> f_system_data(1);
    f_system_data[0] = fsnAct_data[0];
    ib_post_processor->registerVectorVariable(
       "f", MONOMIAL, CONSTANT, IBFEPostProcessor::deformed_material_axis_fcn,f_system_data);
    vector<SystemData> s_system_data(1);
    s_system_data[0] = fsnAct_data[1];
    ib_post_processor->registerVectorVariable(
           "s", MONOMIAL, CONSTANT, IBFEPostProcessor::deformed_material_axis_fcn,s_system_data);

    // 3) Fiber and sheet stretches.
    ib_post_processor->registerScalarVariable(
               "lambda_f", MONOMIAL, CONSTANT, IBFEPostProcessor::material_axis_stretch_fcn,f_system_data);
    ib_post_processor->registerScalarVariable(
               "lambda_s", MONOMIAL, CONSTANT, IBFEPostProcessor::material_axis_stretch_fcn,s_system_data);

    // 4) Cauchy stress sigma = (1/J) PP FF^T.
    ib_post_processor->registerTensorVariable(
                   "sigma_dev", MONOMIAL, CONSTANT,
                   IBFEPostProcessor::cauchy_stress_from_PK1_stress_fcn,
                    PK1_dev_stress_data.system_data, &PK1_dev_stress_data);
    ib_post_processor->registerTensorVariable(
                    "sigma_dil", MONOMIAL, CONSTANT,
                    IBFEPostProcessor::cauchy_stress_from_PK1_stress_fcn,
                    PK1_dil_stress_data.system_data, &PK1_dil_stress_data);


    //5) Eulerian pressure p_f.
    Pointer<hier::Variable<NDIM> > p_var = navier_stokes_integrator->getPressureVariable();
    Pointer<VariableContext> p_current_ctx = navier_stokes_integrator->getCurrentContext();
    HierarchyGhostCellInterpolation::InterpolationTransactionComponent p_ghostfill(
        /*data_idx*/ -1,
        "LINEAR_REFINE",
        /*use_cf_bdry_interpolation*/ false,
        "CONSERVATIVE_COARSEN",
        "LINEAR");
    FEDataManager::InterpSpec p_interp_spec("PIECEWISE_LINEAR",
                                                            QGAUSS,
                                                            FIFTH,
                                                            /*use_adaptive_quadrature*/ false,
                                                            /*point_density*/ 2.0,
                                                            /*use_consistent_mass_matrix*/ true,
                                                          /*use_nodal_quadrature*/ false);

    ib_post_processor->registerInterpolatedScalarEulerianVariable("p_f", LAGRANGE, FIRST, p_var, p_current_ctx,
                                                              p_ghostfill, p_interp_spec);

		// Create Eulerian boundary condition specification objects when needed.
    const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
    std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
    const bool periodic_boundaries = periodic_shift.min() > 0;
    if (!periodic_boundaries)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            ostringstream bc_coefs_name_stream;
            bc_coefs_name_stream << "u_bc_coefs_" << d;
            const string bc_coefs_name = bc_coefs_name_stream.str();
            ostringstream bc_coefs_db_name_stream;
            bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
            const string bc_coefs_db_name = bc_coefs_db_name_stream.str();
            u_bc_coefs[d] = new muParserRobinBcCoefs(bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
        }
        navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
    }

    // Set up visualization plot file writers.
    Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
    if (uses_visit)
    {
        time_integrator->registerVisItDataWriter(visit_data_writer);
    }
    UniquePtr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : NULL);

    // Initialize FE data.
    ib_method_ops->initializeFEData();
    ib_post_processor->initializeFEData();

    // Setup the material axes.
    ModelInitialization::initialize_material_axes(mesh, equation_systems, input_db);

    // Initialize hierarchy configuration and data on all patches.
    time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

    // Deallocate initialization objects.
    app_initializer.setNull();

    // Print the input database contents to the log file.
    plog << "Input database:\n";
    input_db->printClassData(plog);

    // Write out initial visualization data.
    int iteration_num = time_integrator->getIntegratorStep();
    double loop_time = time_integrator->getIntegratorTime();
    if (dump_viz_data)
    {
        pout << "\nWriting visualization files...\n\n";
        if (uses_visit)
        {
            // pout << "\n set up the plot data for visit\n";
            time_integrator->setupPlotData();
            // pout << "\n finish setting up the plot data for visit\n";
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            // pout << "\n finish writing out the initial visit file \n";
        }
        if (uses_exodus)
        {
            ib_post_processor->postProcessData(loop_time);
            exodus_io->write_timestep(exodus_filename, *equation_systems, iteration_num/viz_dump_interval+1, loop_time);
        }
    }
    pout << "\nWriting visualization files done...\n\n";

    //time points for simualtion
    BoundaryConditions::t_end_diastole = input_db->getDouble("TIME_END_DIASTOLE");




    // read the endo_list for LV cavity calculation;
    //BoundaryConditions::readingPoints(mesh);
    //BoundaryConditions::updatePointsPosition(equation_systems);
    BoundaryConditions::readingPointsGeneral(mesh,
                                    BoundaryConditions::LV_endo_points_list,
                                    BoundaryConditions::LV_endo_points,
                                    BoundaryConditions::LV_NoOfEndoNode,
                                    input_db->getString("ENDO_POINTS_LIST"));
    BoundaryConditions::updatePointsPositionGeneral(equation_systems,
                                BoundaryConditions::LV_endo_points_list,
                                BoundaryConditions::LV_endo_points,
                                BoundaryConditions::LV_NoOfEndoNode);

    if (0 == mesh.processor_id())
    {

		    BoundaryConditions::LV_volume = tetVolumeCalculationByPoints(BoundaryConditions::LV_endo_points,
                                        BoundaryConditions::LV_NoOfEndoNode);

		    std::ofstream pv_data_file;
		    pv_data_file.open("pressure_volume_data.dat", std::ofstream::out);
		    pv_data_file << "#time \t" <<"pressure \t"<<"LV volume\t"<<std::endl;
		    pv_data_file <<loop_time<<"    "<<BoundaryConditions::P_current_loading<<"    "<<BoundaryConditions::LV_volume<<std::endl;
		    pv_data_file.close();
	}
  	MPI_Bcast(&BoundaryConditions::LV_volume, 1, MPI_DOUBLE, 0, init.comm().get());

    printf("processor %d, ini volume is %f\n", mesh.processor_id(), BoundaryConditions::LV_volume);
    MPI_Barrier(init.comm().get());

    // Main time step loop.
    pout << "Entering main time step loop...\n";
    const double loop_time_end = time_integrator->getEndTime();
    double dt = 0.0;


    auto perf_start_time = std::chrono::system_clock::now();
    auto perf_end_time   = std::chrono::system_clock::now();
    std::chrono::duration<double> perf_elapsed_seconds = perf_end_time - perf_start_time;

    while (!MathUtilities<double>::equalEps(loop_time,loop_time_end) && time_integrator->stepsRemaining())
    {
        iteration_num = time_integrator->getIntegratorStep();

        pout <<                                                       endl;
        pout << "++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        pout << "At beginning of timestep # " << iteration_num  << endl;
        pout << "Simulation time is " << loop_time                 << endl;


        pout << "\nloading pressure beginning: " << BoundaryConditions::loading_pressure(loop_time) << " mmHg" << endl << endl;

        dt = time_integrator->getMaximumTimeStepSize();
        if (loop_time>BoundaryConditions::t_end_diastole)
	{
		double dtN = input_db->getDouble("DT_SYSTOLE");
		if (dtN<dt)
		{
			dt = dtN;
		}
	}
        BoundaryConditions::boundary_info = mesh.boundary_info.get();

        MechanicsModel::I1_dev_max = -1.0e8;
        MechanicsModel::I1_dev_min = +1.0e8;
        MechanicsModel::I1_dil_max = -1.0e8;
        MechanicsModel::I1_dil_min = +1.0e8;

        MechanicsModel::J_dev_max = -1.0e8;
        MechanicsModel::J_dev_min = +1.0e8;
        MechanicsModel::J_dil_max = -1.0e8;
        MechanicsModel::J_dil_min = +1.0e8;

        time_integrator->advanceHierarchy(dt);

        pout << "I1_dev_max = " << MechanicsModel::I1_dev_max << "\n"
             << "I1_dev_min = " << MechanicsModel::I1_dev_min << "\n"
             << "I1_dil_max = " << MechanicsModel::I1_dil_max << "\n"
             << "I1_dil_min = " << MechanicsModel::I1_dil_min << "\n";

        pout << "J_dev_max = " << MechanicsModel::J_dev_max << "\n"
             << "J_dev_min = " << MechanicsModel::J_dev_min << "\n"
             << "J_dil_max = " << MechanicsModel::J_dil_max << "\n"
             << "J_dil_min = " << MechanicsModel::J_dil_min << "\n";

        pout << "Current Pressure Loading =  "<< BoundaryConditions::P_current_loading<<"\n";
        perf_end_time   = std::chrono::system_clock::now();
        perf_elapsed_seconds = perf_end_time - perf_start_time;
        pout << "Computational time used = "<< perf_elapsed_seconds.count() << " seconds \n";

        if (MechanicsModel::enable_active_tension)
        {
            ActiveContraction::update_active_tension_model_state_variables(equation_systems, loop_time, dt);
        }

		//calculate the LV cavity volume for pressure update
        BoundaryConditions::updatePointsPositionGeneral(equation_systems,
                                BoundaryConditions::LV_endo_points_list,
                                BoundaryConditions::LV_endo_points,
                                BoundaryConditions::LV_NoOfEndoNode);
        if (0 == mesh.processor_id())
        {
		   BoundaryConditions::LV_volume = tetVolumeCalculationByPoints(BoundaryConditions::LV_endo_points, BoundaryConditions::LV_NoOfEndoNode);
		   std::ofstream pv_data_file;
		   pv_data_file.open("pressure_volume_data.dat", std::ofstream::out|std::ofstream::app);
		   pv_data_file <<loop_time<<"    "<<BoundaryConditions::P_current_loading<<"    "<<BoundaryConditions::LV_volume<<std::endl<<std::flush;
		   pv_data_file.close();
	     }
	    MPI_Bcast(&BoundaryConditions::LV_volume, 1, MPI_DOUBLE, 0, init.comm().get());
        if (init.comm().rank() == init.comm().size()-1)
            printf("processor %d, LV cavity volume is %f\n", mesh.processor_id(), BoundaryConditions::LV_volume);
        MPI_Barrier(init.comm().get());

        loop_time += dt;
        pout <<                                                       endl;
        pout << "At end       of timestep # " << iteration_num + 1 << endl;
        pout << "Simulation time is " << loop_time                 << endl;
        pout << "LV cavity volume is "<<BoundaryConditions::LV_volume<< endl;
        pout << "++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        pout <<                                                       endl;

        // At specified intervals, write visualization and restart files,
        // print out timer data, and store hierarchy data for post
        // processing.
        iteration_num += 1;
        const bool last_step = !time_integrator->stepsRemaining();

	bool ED_step = false;
	if (loop_time <= BoundaryConditions::t_load && loop_time +dt >  BoundaryConditions::t_load)
	{
	  ED_step = true;
	}

        if (dump_viz_data && (iteration_num%viz_dump_interval == 0 || last_step || ED_step))
        {
            pout << "\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (uses_exodus)
            {
                ib_post_processor->postProcessData(loop_time);
                exodus_io->write_timestep(exodus_filename, *equation_systems, iteration_num/viz_dump_interval+1, loop_time);
            }
        }
        if (dump_restart_data && (iteration_num%restart_dump_interval == 0 || last_step|| ED_step))
        {
            pout << "\nWriting restart files...\n\n";
            RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            ib_method_ops->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
        }
        if (dump_timer_data && (iteration_num%timer_dump_interval == 0 || last_step|| ED_step))
        {
            pout << "\nWriting timer data...\n\n";
            TimerManager::getManager()->print(plog);
        }

        //will try to figure out how to calculate the LV cavity volume
        //const Parallel::Communicator& comm_temp = (*equation_systems).comm();
        //MPI_Barrier( equation_systems->comm().get() );
        //MPI_Barrier(init.comm().get());
        //printf("From processor %d\n", equation_systems->processor_id());
        //printf("From processor %d\n", (*equation_systems).comm().rank());

    }// loop

    // Clean up dynamically allocated objects.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        delete u_bc_coefs[d];
    }
}// cleanup dynamically allocated objects prior to shutdown

 //SAMRAIManager::shutdown();
 //PetscFinalize(); //seems there is no need to use PetscFinalize() since using LibMeshInit to manage
    return 0;
}// main_driver
