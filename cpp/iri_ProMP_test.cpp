ri_ProDMPtest.cpp
 *
 *  Created on: Aug 29, 2013
 *      Author: robot
 */

#include <iostream>
#include <string>

#include <barrett/units.h>
#include <barrett/systems.h>
#include <barrett/products/product_manager.h>
#include <Eigen/Core>
#include "iri_wam.hpp"
#include "ac_promp.h"

using namespace barrett;
using detail::waitForEnter;

BARRETT_UNITS_TYPEDEFS(7);

int main(int argc, char** argv)
{
    ProductManager pm;
    pm.waitForWam(true);
    pm.wakeAllPucks();

    /* Home Made getWam */
    if ( !pm.foundWam7() )
        throw std::logic_error("ProductManager::getWam7(): No WAM7 was found on the bus.");

    /* The iri_WAM */
    barrett::systems::iri::WAM<7> iri_wam(pm.getExecutionManager(),
            pm.getWamPucks(),
            pm.getSafetyModule(),
            pm.getConfig().lookup(pm.getWamDefaultConfigPath()));

    /* Initialization & Activation*/
    pm.startExecutionManager();
    pm.getSafetyModule()->waitForMode(SafetyModule::ACTIVE, true, 0.05);
    std::cout << "Ready to gravity ..." << std::endl;
    std::cout << "Press ENTER!" << std::endl;
    waitForEnter();

    /* Ready to work */
    iri_wam.gravityCompensate();

    // ------------------------------------------------------------------ 
    // LOCAL VARIABLES + pre-processing

    int b_explore(false);
    bool b_isCartesian(true);
    bool b_isCondition(false);
    std::string traj_filename("testgerard");
    std::vector<float> pose; // [0-1] 
    pose.resize(7);
    // Original position was (0.6368, 0.1252, -0.2237)
    pose[0] = 0.5168;
    pose[1] = 0.1252;
    pose[2] = -0.2137;

    float point_of_change(0.5f); // [0-1] 
    float cov_ratio(0.3f);
    std::vector<Eigen::VectorXd> eigen_trajectory_points;
    std::vector<double> local_times;

    ac_math::compute_promp(b_explore, b_isCartesian, b_isCondition, traj_filename, pose, point_of_change, cov_ratio, eigen_trajectory_points, local_times);

    // conversion from vector of 
    std::vector<jp_type> trajectory_points;
    std::cout << "Size trajectory_points: " << eigen_trajectory_points.size() << std::endl;
    for(int ii=0; ii<eigen_trajectory_points.size(); ii++) {
        trajectory_points.push_back(eigen_trajectory_points[ii]);
    }

    // ------------------------------------------------------------------ 
    // LOG ALL THE DATA
    // bin file for log
    const char BIN_FILE[] =   "./tmp/sequence_state.bin";

    // flag to control the file save
    bool save_mode = true;
    typedef boost::tuple<jp_type, jv_type, jt_type> state_type;
    systems::PeriodicDataLogger<state_type> dl(pm.getExecutionManager(),
            new log::RealTimeWriter<state_type>(BIN_FILE, pm.getExecutionManager()->getPeriod()));

    // Grouper for Log
    systems::TupleGrouper<jp_type, jv_type, jt_type> mi_TG("State-Action Grouper");
    /* connect(my_ramp.output,mi_TG.template getInput<0> ());
       connect(my_traj.posOutput,mi_TG.template getInput<1> ());
       connect(my_traj.velOutput,mi_TG.template getInput<2> ());
       connect(my_traj.accelOutput,mi_TG.template getInput<3> ()); */
    connect(iri_wam.jpOutput, mi_TG.getInput<0> ());
    connect(iri_wam.jvOutput, mi_TG.getInput<1> ());
    connect(iri_wam.supervisoryController.output, mi_TG.getInput<2> ());
    /* Ready for LOG */
    connect(mi_TG.output, dl.input);

    // TRAJECTORIES EXECUTIONS
    // ------------------------------------------------------------------ 

    std::cout << "Ready to move ..." << std::endl;
    std::cout << "Press ENTER!" << std::endl;
    waitForEnter();

    //iri_wam.moveTo(trajectory_points[0], 0.75, 0.75);
    //iri_wam.moveToTrajectory(trajectory_points);
    iri_wam.moveToTrajectory(trajectory_points, false, local_times);
    std::cout << "Press ENTER!" << std::endl;
    waitForEnter();
    //iri_wam.moveToTrajectory(trajectory_points, false, local_times);
    iri_wam.moveToTrajectoryCompliant(trajectory_points, false, local_times);
    //iri_wam.moveToTrajectoryCompliant(trajectory_points);

    // Wait for the user to press Shift-idle
    std::cout << "End of testing" << std::endl;
    std::cout << "Press ENTER to go Home!" << std::endl;
    waitForEnter();
    iri_wam.moveHome();

    // ------------------------------------------------------------------ 
    // LOG ALL THE DATA
    // Saving and Closing log Files
    std::string output_fileName = "promp_test_state_sequence.txt";
    std::cout << "Saving... " << output_fileName << " " << std::endl;
    dl.closeLog();
    log::Reader<state_type> lr(BIN_FILE);
    lr.exportCSV(output_fileName.c_str());

    // ------------------------------------------------------------------ 
    // WAIT UNTIL IDLE
    pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);
    return 0;
}
