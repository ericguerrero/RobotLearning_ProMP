/*iri_learnProMP.cpp
 * Learn (by demonstration) a parametrized movement with Probabilistic Motor Primitives. This computes
 * the weights, centers and widths of the gaussians.
 *
 * Author : Adrià Colomé
 * Date   : March 2015
 *
 *
 *   ToDo : 
 *          
 */


#include <iostream>
#include <fstream>
#include <string>

#include <barrett/units.h>
#include <barrett/log.h>
#include <barrett/systems.h>
#include <barrett/products/product_manager.h>

#include <barrett/standard_main_function.h>
#include "dp_ja_type_system.hpp"
#include "iri_trajectory_recorder.hpp"
#include "iri_wamik.h"
#include "ac_promp.h"
#include "ac_math.h"
#include <Eigen/Core>
#include <Eigen/SVD>

using namespace barrett;
using detail::waitForEnter;


template<size_t DOF>
int wam_main(int argc, char** argv, ProductManager& pm, systems::Wam<DOF>& wam) {
    BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);
    // Check for trajectory names
    if(argv[1]==NULL )
    {
        std::cout << "Missing file name  [Trajectory]" << std::endl;
        return 1;
    }
    std::string trajectory_filename(argv[1]);

    // NUMBER OF ACTUATED DMP
    int rdof=7;
    // bin file for log
    const char BIN_FILE[] =   "./tmp/sequence_state.bin";
    // flag to control the file save
    bool save_mode = true;
    //Compensate Gravity!
    wam.gravityCompensate();
    systems::Ramp my_ramp(pm.getExecutionManager());
    recorder<DOF> traj_rec(pm.getExecutionManager(),"trajectoryrecorder");
    connect(wam.jpOutput,traj_rec.jp_Input);
    connect(wam.jvOutput,traj_rec.jv_Input);
    traj_rec.SetDataRatio(50);


    std::cout<<"\n Enter the number of trajectories that will generate the ProMP=";
    //waitForEnter();
    int Ndemos;
    int foo;
    //std::cin>>Ndemos;
    std::string line;
    std::getline(std::cin, line);
    std::cout<<"\n press enter\n ";
    waitForEnter();
    //istringstream ( line ) >> Ndemos;
    Ndemos=atoi(line.c_str());
    usleep(1000000);

    Eigen::MatrixXd JP,JP2,JV,JV2,Y0,JPout,JVout,Tout,T0;
    Eigen::VectorXd w;
    std::vector< Eigen::MatrixXd > Demos_JP;
    std::vector< Eigen::MatrixXd > Demos_JV;
    int timeratio=traj_rec.GetDataRatio();


    int solved,Nred;
    Nred=1;
    /*
       std::cout<<"to load now"<<flush;
    //Y0=LoadDataMatrix("rawW1",718,7);
    Y0=LoadDataMatrix("test0",89,7);
    Demos_JP.push_back(Y0);
    Eigen::MatrixXd Yk,Yk_s;
    std::cout<<"loaded 0"<<flush;

    Yk=LoadDataMatrix("test1",89,7);
    std::cout<<"loaded 1";
    solved=ac_math::IterativeTimeWarping(Y0,Yk,w,JP2,timeratio);
    ac_math::SaveMatrixAsFile(JP2,"shifted",1);	
    Demos_JP.push_back(JP2);

    Yk=LoadDataMatrix("test2",89,7);
    std::cout<<"loaded 1";
    solved=ac_math::IterativeTimeWarping(Y0,Yk,w,JP2,timeratio);
    ac_math::SaveMatrixAsFile(JP2,"shifted",2);	
    Demos_JP.push_back(JP2);
    */


    // GET DATA FROM DEMONSTRATIONS
    for (int i=0;i<Ndemos;i++){
        std::cout<<"Press enter to START recording a desired trajectory (first move the robot to the initial position) \n";
        std::getline(std::cin, line);

        std::cout<<"started"<<std::endl;
        traj_rec.start();
        my_ramp.start();
        // teach motion to robot	
        std::cout<<"Press enter to STOP recording the desired trajectory"<<std::endl;
        std::getline(std::cin, line);
        std::cout<<"stopped"<<std::endl<<std::flush;

        traj_rec.stop();
        my_ramp.stop();
        JP=traj_rec.GetPosMatrix();
        JV=traj_rec.GetVelMatrix();

        ac_math::ExtractDataExperiment(JP,JV,timeratio,JPout,JVout,Tout);
        if(i==0){
            T0=Tout;
            Y0=JPout;
            Demos_JP.push_back(JPout);
            Demos_JV.push_back(JVout);
            ac_math::SaveMatrixAsFile(JPout,"test",i);
        }else{
            ac_math::SaveMatrixAsFile(JPout,"test",i);	
            //			std::cout<<"\n JPout=\n "<<JPout;
            solved=ac_math::IterativeTimeWarping(Y0,JPout,w,JP2,timeratio);
            //			if(solved==0){
            //				Nred=1;
            //				std::cout<<"Didn't work, recomputing"<<flush;
            //				solved=ac_math::IterativeTimeWarping(Y0,JPout,w,JP2,Nred,timeratio);
            //			}	
            ac_math::SaveMatrixAsFile(JP2,"shifted",i);		
            Demos_JP.push_back(JP2);

        }

    }

    std::cout<<"\n Time warping completed. Trajectories stored.\n ";

    Eigen::MatrixXd Jpaux;
    // ProMP parameters
    int cartesiancontrol,d,Nf;
    Nf=10;

    std::cout<<"\n Do you want to store trajectory as a cartesian trajectory (press 1 and enter) or joint trajectory (press 0 and enter)"<<std::endl;
    std::cin>>cartesiancontrol;
    Eigen::VectorXd Cs(Nf);
    double h_s;
    // define ProMP centers and widths, we set final_time=1 and we'll rescale it later.
    for(int s=0;s<Nf;s++){
        Cs(s)=((double) s)/((double) Nf);
    }
    h_s=1.0/((double) Nf)/4.0;
    Eigen::MatrixXd weights_aux,weights,DATAcart;

    if(cartesiancontrol==0){
        std::cout<<"\n You chose joint trajectory";
        d=7;		
        weights.resize(Nf*d,Ndemos);
        weights.fill(0.0);
    }else{
        std::cout<<"\n You chose cartesian trajectory";
        // for the cartesian version
        d=6;
        weights.resize(Nf*d,Ndemos);
        weights.fill(0.0);

    }
    WAMIK WAM_kinematics;
    Eigen::VectorXd q_aux(7),u_aux,u_aux_old;
    u_aux.resize(3);
    u_aux.fill(0.0);
    u_aux_old.resize(3);
    u_aux_old.fill(0.0);
    Eigen::MatrixXd T_aux,R_aux;
    R_aux.resize(3,3);
    //usleep(5000000);
    Eigen::VectorXd Y0mean;
    Y0mean.resize(7);
    Y0mean.fill(0.0);
    for(int k=0;k<Ndemos;k++){
        // compute weights of k-th trajectory
        Jpaux=Demos_JP[k];
        std::cout<<"\n CHECK 1, k="<<k<<"\n"<<std::flush;
        for(int joi=0;joi<7;joi++){		
            Y0mean(joi)=Y0mean(joi)+Jpaux(0,joi)/Ndemos;
        }
        if (cartesiancontrol==0){ // control in joints, we only consider joint space trajectory
            weights_aux=ac_math::FitParamsFromTraj(Cs,h_s,Jpaux);
            weights.block(0,k,Nf*d,1)=weights_aux.block(0,0,Nf*d,1); // every column is one experiment.				
        }else{ //control in cartesian space, we need their cartesian expression.
            std::cout<<"\n CHECK 2"<<std::flush;
            DATAcart.resize(Jpaux.rows(),d);
            for(int i=0;i<Jpaux.rows();i++){
                // convert Jpaux.row(i) to cartesian trajectory (Homogeneous transformation) with forward kinematics
                for(int j=0;j<7;j++){
                    q_aux(j)=Jpaux(i,j);							
                }
                std::cout<<"\n CHECK 3"<<std::flush;
                T_aux=WAM_kinematics.forwardkinematics(q_aux);
                std::cout<<"\n CHECK 4"<<std::flush;
                // convert homogeneous transformation to state vector
                R_aux.block(0,0,3,3)=T_aux.block(0,0,3,3);
                u_aux_old=u_aux;
                u_aux=ac_math::Quat2TriVec(ac_math::Rot2Quat(R_aux));
                std::cout<<"\n CHECK 5"<<std::flush;
                double v1,v2;
                std::cout<<"\n u_aux="<<u_aux<<std::flush;
                std::cout<<"\n u_aux2="<<u_aux_old<<std::flush;
                v1=(std::abs(u_aux_old(0)-u_aux(0))+std::abs(u_aux_old(1)-u_aux(1))+std::abs(u_aux_old(2)-u_aux(2)));
                std::cout<<"\n CHECK 5a"<<std::flush;
                v2=(std::abs(u_aux_old(0)+u_aux(0))+std::abs(u_aux_old(1)+u_aux(1))+std::abs(u_aux_old(2)+u_aux(2)));
                std::cout<<"\n CHECK 5b"<<std::flush;
                if (i>1 && v1>v2+0.01){
                    u_aux=-u_aux;
                    u_aux_old=u_aux;
                    std::cout<<"\n Changed sign in k,i="<<k<<", "<<i<<std::endl;
                }
                std::cout<<"\n CHECK 6"<<std::flush;
                //					u_aux=Rot2Quat(R_aux);	
                // store in DATAcart	
                DATAcart(i,0)=T_aux(0,3);// x coord
                DATAcart(i,1)=T_aux(1,3);// y coord
                DATAcart(i,2)=T_aux(2,3);// z coord
                DATAcart(i,3)=u_aux(0);// rx coord
                DATAcart(i,4)=u_aux(1);// ry coord
                DATAcart(i,5)=u_aux(2);// rz coord				
                //					DATAcart(i,6)=u_aux(3);// rz coord		
                std::cout<<"\n CHECK 7"<<std::flush;
            }
            std::cout<<"\n CHECK 8"<<std::flush;
            weights_aux=ac_math::FitParamsFromTraj(Cs,h_s,DATAcart);
            weights.block(0,k,Nf*d,1)=weights_aux.block(0,0,Nf*d,1); 		
            ac_math::SaveMatrixAsFile(DATAcart,"Cartesian",k);		
            std::cout<<"\n CHECK 9"<<std::flush;
        }
    }
    std::cout<<"\n Weights=\n"<<weights;

    // fit parameters with all the weights
    Eigen::MatrixXd Sw,mw;
    ac_math::FitNormalDistribution(weights,Sw,mw);


    // STORE PROMP DATA AS FILES
    std::string SfileName(trajectory_filename),mfileName(trajectory_filename);
    SfileName +="Sw";
    mfileName +="mw";
    ac_math::SaveMatrixAsFile(Sw,SfileName.c_str(),-1.0);
    ac_math::SaveMatrixAsFile(mw,mfileName.c_str(),-1.0);
    // save other things
    ofstream OFile;
    std::string OfileName(trajectory_filename);
    OfileName += "O.txt";
    OFile.open(OfileName.c_str());
    // write C
    OFile<<Cs(0);
    for (int s=1;s<Nf;s++){
        OFile<<","<<Cs(s);
    }
    OFile<<"\n";
    // write Nf
    OFile<<Nf;
    OFile<<"\n";
    // write time
    OFile<<((double)Y0.rows())*0.002*((double)timeratio);
    OFile<<"\n";
    // write width
    OFile<<h_s;
    OFile<<"\n";
    // write dimension
    OFile<<d;
    OFile<<"\n";
    // write initial Joint position
    OFile<<Y0mean(0);
    for (int jo=1;jo<7;jo++){
        OFile<<","<<Y0mean(jo);
    }
    OFile<<"\n";

    // close file
    OFile.close();

    // Saving and Closing Files
    std::cout<<"Press enter to move Home"<<std::endl;
    waitForEnter();
    wam.moveHome(true,0.75);
    wam.idle();
    pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);


    return 0;
}
