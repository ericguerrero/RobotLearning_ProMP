/*
 * iri_learnProMP.cpp
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
// local compile g++ iri_learnProMP_local.cpp -I/usr/local/include/iridrivers/ -I/usr/include/eigen2/ -std=c++11 -L/usr/local/lib/iridrivers/ -lac_math -liri_wamik -o iri_learnProMpLocal -Wl,-rpath,/usr/local/lib/iridrivers/

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

//#include "dp_ja_type_system.hpp"
#include "ac_promp.h"
#include "ac_math.h"
#include <Eigen/Core>
#include <Eigen/SVD>


int main(int argc, char** argv) {
    // Check for trajectory names
    if(argv[1] == NULL) {
        std::cout << "Missing file name  [Trajectory]" << std::endl;
        return 1;
    }
    else if (argc < 5) {
        std::cerr << "Wrong number of parameters!!" << std::endl;
        std::cerr << "Usage: ./iri_learnProMpLocal <trajectory_params_name> <cartesian_trajectories_path> <number_of_demos> <number_of_gaussians>" << std::endl;
        return -1;
    }
    std::string trajectory_filename(argv[1]);
    std::string cartesianpath(argv[2]);
    int Ndemos=atoi(argv[3]);
    int Nf = atoi(argv[4]);
    std::cout << "Trajectory name: " << trajectory_filename << std::endl;
    std::cout << "Casrtesian trajectory files path: " << cartesianpath << std::endl;
    std::cout << "Number of demos: " << Ndemos << std::endl;
    std::cout << "Number of gaussians: " << Nf << std::endl;
    //Compensate Gravity!

    /*std::cout << "Enter the number of trajectories that will generate the ProMP = ";
    //waitForEnter();
    int Ndemos;
    int foo;
    //std::cin>>Ndemos;
    std::string line;
    std::getline(std::cin, line);*/

    //Ndemos=atoi(line.c_str());

    Eigen::MatrixXd JP, JP2, JV, JV2, Y0, JPout, JVout, Tout, T0;
    Eigen::VectorXd w;
    //std::vector<Eigen::MatrixXd> Demos_JP;
    //std::vector<Eigen::MatrixXd> Demos_JV;

    int solved, Nred;
    Nred = 1;

    /*// GET DATA FROM DEMONSTRATIONS
    for (int i=0; i<Ndemos; i++){
        std::cout << "Press enter to START recording a desired trajectory (first move the robot to the initial position)" << std::endl;
        std::getline(std::cin, line);
        std::cout << "started" << std::endl;

        traj_rec.start();
        my_ramp.start();
        // teach motion to robot	
        std::cout << "Press enter to STOP recording the desired trajectory" << std::endl;
        std::getline(std::cin, line);

        std::cout << "stopped" << std::endl << std::endl;
        traj_rec.stop();
        my_ramp.stop();

        JP = traj_rec.GetPosMatrix();
        JV = traj_rec.GetVelMatrix();
        ac_math::ExtractDataExperiment(JP, JV, timeratio, JPout, JVout, Tout);

        if(i==0){
            T0 = Tout;
            Y0 = JPout;
            Demos_JP.push_back(JPout);
            Demos_JV.push_back(JVout);
            ac_math::SaveMatrixAsFile(JPout,"test",i);
        }else{
            ac_math::SaveMatrixAsFile(JPout,"test",i);	
            //	std::cout<<"\n JPout=\n "<<JPout;
            solved=ac_math::IterativeTimeWarping(Y0,JPout,w,JP2,timeratio);
            //	if(solved==0){
            //	  Nred=1;
            //	  std::cout<<"Didn't work, recomputing"<<endl;
            //	  solved=ac_math::IterativeTimeWarping(Y0,JPout,w,JP2,Nred,timeratio);
            //	}	
            ac_math::SaveMatrixAsFile(JP2,"shifted",i);		
            Demos_JP.push_back(JP2);
        }
    }

    std::cout << "Time warping completed. Trajectories stored." << std::endl;*/

    ////Eigen::MatrixXd Jpaux;
    // ProMP parameters

    // define ProMP centers and widths, we set final_time=1 and we'll rescale it later.
    Eigen::VectorXd Cs(Nf);
    for(int s=0; s<Nf; s++) {
        Cs(s)=((double) s)/((double) Nf);
    }

    double h_s;
    h_s = 1.0/(4.0f*double(Nf));
    
    Eigen::MatrixXd weights_aux, weights, DATAcart;

    int cartesiancontrol, d;
    //std::cout << "Do you want to store trajectory as a cartesian trajectory (press 1 and enter) or joint trajectory (press 0 and enter)" << std::endl;
    //std::cin >> cartesiancontrol;
    cartesiancontrol = 1; // FIXME note this as the test just was using cartesian
    if(cartesiancontrol==0) {
        std::cout<<"\n You chose joint trajectory";
        d=7;		
        weights.resize(Nf*d, Ndemos);
        weights.fill(0.0);
    }
    else {
        std::cout<<"\n You chose cartesian trajectory";
        // for the cartesian version
        d=6;
        weights.resize(Nf*d,Ndemos);
        weights.fill(0.0);
    }
    std::cout << std::endl;

    //WAMIK WAM_kinematics;
    /*Eigen::VectorXd q_aux(7),u_aux,u_aux_old;
    u_aux.resize(3);
    u_aux.fill(0.0);
    u_aux_old.resize(3);
    u_aux_old.fill(0.0);
    Eigen::MatrixXd T_aux, R_aux;
    R_aux.resize(3,3);*/
    //usleep(5000000);
    /*Eigen::VectorXd Y0mean;
    Y0mean.resize(7);
    Y0mean.fill(0.0);*/
    for(int k=0; k<Ndemos; k++){
        // compute weights of k-th trajectory
        /*Jpaux=Demos_JP[k];
        std::cout << "CHECK 1, k= " << k << std::endl;

        for(int joi=0;joi<7;joi++){	
            Y0mean(joi)=Y0mean(joi)+Jpaux(0,joi)/Ndemos;
        }*/

        /*if (cartesiancontrol==0){ // control in joints, we only consider joint space trajectory
            weights_aux=ac_math::FitParamsFromTraj(Cs,h_s,Jpaux);
            weights.block(0,k,Nf*d,1) = weights_aux.block(0,0,Nf*d,1); // every column is one experiment.				
        }else{ //control in cartesian space, we need their cartesian expression.
            /*std::cout << "CHECK 2" << std::endl;
            DATAcart.resize(Jpaux.rows(),d);
            for(int i=0; i<Jpaux.rows(); i++){
                // convert Jpaux.row(i) to cartesian trajectory (Homogeneous transformation) with forward kinematics
                for(int j=0; j<7; j++){
                    q_aux(j)=Jpaux(i,j);							
                }
                std::cout << "CHECK 3" << std::endl;
                T_aux=WAM_kinematics.forwardkinematics(q_aux);
                std::cout << "CHECK 4" << std::endl;
                // convert homogeneous transformation to state vector
                R_aux.block(0,0,3,3)=T_aux.block(0,0,3,3);
                u_aux_old=u_aux;
                u_aux=ac_math::Quat2TriVec(ac_math::Rot2Quat(R_aux));
                std::cout << "CHECK 5" << std::endl;
                double v1,v2;
                std::cout << "u_aux= " << u_aux << std::endl;
                std::cout << "u_aux2= " << u_aux_old << std::endl;
                v1=(std::abs(u_aux_old(0)-u_aux(0))+std::abs(u_aux_old(1)-u_aux(1))+std::abs(u_aux_old(2)-u_aux(2)));
                std::cout << "CHECK 5a" << std::endl;
                v2=(std::abs(u_aux_old(0)+u_aux(0))+std::abs(u_aux_old(1)+u_aux(1))+std::abs(u_aux_old(2)+u_aux(2)));
                std::cout << "CHECK 5b" << std::endl;
                if (i>1 && v1>v2+0.01){
                    u_aux=-u_aux;
                    u_aux_old=u_aux;
                    std::cout << "Changed sign in k,i= " << k << ", " << i << std::endl;
                }
                std::cout << "CHECK 6" << std::endl;
                //u_aux=ac_math::Rot2Quat(R_aux);	
                // store in DATAcart	
                DATAcart(i,0) = T_aux(0,3);// x coord
                DATAcart(i,1) = T_aux(1,3);// y coord
                DATAcart(i,2) = T_aux(2,3);// z coord
                DATAcart(i,3) = u_aux(0);// rx coord
                DATAcart(i,4) = u_aux(1);// ry coord
                DATAcart(i,5) = u_aux(2);// rz coord				
                //					DATAcart(i,6)=u_aux(3);// rz coord		
                std::cout << "CHECK 7" << std::endl;
            }
            std::cout << "CHECK 8" << std::endl;*/

            // GCC: Count rows and columns and load data
            std::string cartesianfname = cartesianpath + "/Cartesian" + std::to_string(k);
            int rows = 0;
            std::ifstream in(cartesianfname + ".txt");
            if (!in.is_open()) {
                std::cerr << "Error: Cartesian traj file \""  << cartesianfname << ".txt\" was not found." << std::endl;
                return -1;
            }
            std::string unused;
            std::getline(in, unused);
            int cols = std::count(unused.begin(), unused.end(), ',') + 1;
            ++rows;
            while (std::getline(in, unused)) ++rows;
            
            std::cout << cartesianfname << " -- rows: " << rows << " cols: " << cols << " -- " << unused << std::endl;
            DATAcart = ac_math::LoadDataMatrix(cartesianfname, rows, cols);
            // GCC
            weights_aux = ac_math::FitParamsFromTraj(Cs,h_s,DATAcart);
            weights.block(0,k,Nf*d,1) = weights_aux.block(0,0,Nf*d,1); 		
            //ac_math::SaveMatrixAsFile(DATAcart,"Cartesian",k);		
            
        //}
    }
    std::cout << "Weights = " << weights << std::endl;

    // Fit parameters with all the weights
    Eigen::MatrixXd Sw,mw;
    ac_math::FitNormalDistribution(weights,Sw,mw);

    // STORE PROMP DATA AS FILES
    std::string SfileName(trajectory_filename), mfileName(trajectory_filename);
    SfileName += "Sw";
    mfileName += "mw";
    ac_math::SaveMatrixAsFile(Sw, SfileName.c_str(), -1.0);
    ac_math::SaveMatrixAsFile(mw, mfileName.c_str(), -1.0);

    // Save other things
    /*std::ofstream OFile;
    std::string OfileName(trajectory_filename);
    OfileName += "O.txt";
    OFile.open(OfileName.c_str());
    // write C
    OFile << Cs(0);
    for (int s=1;s<Nf;s++){
        OFile << "," << Cs(s);
    }
    OFile << "\n";
    // write Nf
    OFile << Nf;
    OFile << "\n";
    // write time
    OFile << ((double)Y0.rows())*0.002*((double)timeratio);
    OFile << "\n";
    // write width
    OFile << h_s;
    OFile << "\n";
    // write dimension
    OFile << d;
    OFile << "\n";
    // write initial Joint position
    OFile << Y0mean(0);
    for (int jo=1;jo<7;jo++){
        OFile << "," << Y0mean(jo);
    }
    OFile << "\n";
    // close file
    OFile.close();*/

    return 0;
}
