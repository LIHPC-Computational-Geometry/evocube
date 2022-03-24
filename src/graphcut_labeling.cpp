
#include "graphcut_labeling.h"
#include <GCoptimization.h>
#include <chrono>

#include "logging.h"
#include "flagging_utils.h"

Eigen::VectorXi graphcutFlagging(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
        const Eigen::MatrixXd& N, const Eigen::MatrixXi& TT,
        const Eigen::VectorXi& locked_flags, const Eigen::VectorXi& forbidden_flags,
        int compact_coeff, int fidelity_coeff){

    // See example.cpp in gco-v3.0
    // function called "GridGraph_DArraySArraySpatVarying"

    // Note: it would probably be worth looking for a more recent library

    std::chrono::steady_clock::time_point before_precomp = std::chrono::steady_clock::now();

    Eigen::MatrixXd axes = axesMatrix();

    Eigen::VectorXi result(F.rows());
    int num_elem = F.rows();
    int num_labels = 6;
	// first set up the array for data costs
	int *data = new int[num_elem*num_labels];
	for ( int i = 0; i < num_elem; i++ ){
        Eigen::RowVector3d nt = N.row(i);
        for (int l = 0; l < num_labels; l++ ){
            double dot = (nt.dot(axes.row(l)) - 1.0)/0.2;
            double cost = 1 - std::exp(-(1./2.)*std::pow(dot,2));
            data[i*num_labels+l] = (int) (fidelity_coeff*100*cost);
        }
    }

    //Use data costs to enforce locked flags
    if (locked_flags.rows() > 0 && locked_flags.maxCoeff() > -1){
        if (locked_flags.minCoeff() == locked_flags.maxCoeff()){
            coloredPrint("Warning in Graph-cut flagging: locked on a single color.", "yellow");
        }
        //std::cout << "Locked percentage: " << 100*((locked_flags.array() != -1).count())/F.rows() << std::endl;
        for (int i = 0; i < num_elem; i++){
            if (locked_flags(i) == -1) continue;
            for (int l = 0; l < num_labels; l++ ){
                if (l!=locked_flags(i)){
                    data[i*num_labels+l] = 10e4;
                }
                else {
                    data[i*num_labels+l] = 0;
                }
            }
        }
    }

    //Use data costs to enforce forbidden flags
    if (forbidden_flags.rows() > 0 && forbidden_flags.maxCoeff() > -1){
        if (forbidden_flags.minCoeff() == forbidden_flags.maxCoeff()){
            coloredPrint("Warning in Graph-cut flagging: single color forbidden.", "yellow");
        }

        for (int i = 0; i < num_elem; i++){
            if (forbidden_flags(i) == -1) continue;
            for (int l = 0; l < num_labels; l++ ){
                if (l==forbidden_flags(i)){
                    data[i*num_labels+l] = 10e4;
                }
            }
        }
    }

	// next set up the array for smooth costs
    // 1 if neighbors, 0 otherwise. Will be 
    // multiplied by neighbor coefficient provided in
    // setNeighbors call.

    bool prevent_opposite_neighbors = false;

	int *smooth = new int[num_labels*num_labels];
	for ( int l1 = 0; l1 < num_labels; l1++ ){
		for (int l2 = 0; l2 < num_labels; l2++ ){ 
			if (l1==l2) smooth[l1+l2*num_labels] = 0;
            else if (axes.row(l1) == - axes.row(l2)) {
                if (prevent_opposite_neighbors)
                    smooth[l1+l2*num_labels] = 10000;
                else smooth[l1+l2*num_labels] = 1;
            }
            else smooth[l1+l2*num_labels] = 1;
        }
    }


	try{
		GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph(num_elem,num_labels);
		gc->setDataCost(data);
		gc->setSmoothCost(smooth);
		
        for (int i=0; i<TT.rows(); i++){
            Eigen::RowVector3d np = N.row(i);
            for (int j=0; j<TT.cols(); j++){
                if (TT(i,j) == -1) continue;
                Eigen::RowVector3d nq = N.row(TT(i,j));
                double dot = (np.dot(nq)-1)/0.25;
                double cost = std::exp(-(1./2.)*std::pow(dot,2));
                gc->setNeighbors(i, TT(i,j), (int) (compact_coeff*100*cost));
            }
        }


        std::chrono::steady_clock::time_point after_precomp = std::chrono::steady_clock::now();   
        //precompute_time = std::chrono::duration_cast<std::chrono::milliseconds>(after_precomp - before_precomp).count();

		//printf("\nBefore optimization energy is %d",gc->compute_energy());
		gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		//printf("\nAfter optimization energy is %d\n",gc->compute_energy());

		for ( int  i = 0; i < num_elem; i++ )
			result[i] = gc->whatLabel(i);

		delete gc;
	}
	catch (GCException e){
		e.Report();
	}

	delete [] smooth;
	delete [] data;

    return result;    
}

Eigen::VectorXi graphcutFlagging(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                                 const Eigen::MatrixXd& N, const Eigen::MatrixXi& TT,
                                 int compact_coeff, int fidelity_coeff){
    //int precompute_time;
    Eigen::VectorXi locked_flags, forbidden_flags;

    return graphcutFlagging(V, F, N, TT, locked_flags, forbidden_flags, 
                            compact_coeff,  fidelity_coeff);
}

std::vector<int> graphcutTurningPoints(const std::vector<int>& bnd, const Eigen::MatrixXd& V,
                                       const Eigen::RowVector3d& desired_dir){
    Eigen::VectorXi result(bnd.size() - 1); //first vertex can't be a turning point
    int num_elem = bnd.size() - 1;
    int num_labels = 2;

	// unary costs
	int *data = new int[num_elem*num_labels];
    Eigen::RowVector3d dir = desired_dir;

	for (int i = 0; i < num_elem; i++ ){
        Eigen::RowVector3d edge = (V.row(bnd[i+1]) - V.row(bnd[i])).normalized();
        for (int l = 0; l < num_labels; l++){
            double dot = (dir.dot(edge.row(l)) - 1.0)/0.9;
            double cost = 1 - std::exp(-(1./2.)*std::pow(dot,2));
            data[i*num_labels+l] = (int) (100*cost);
        }
    }

	// binary costs
	int *smooth = new int[num_labels*num_labels];
	for ( int l1 = 0; l1 < num_labels; l1++ ){
		for (int l2 = 0; l2 < num_labels; l2++ ){ 
			if (l1==l2) smooth[l1+l2*num_labels] = 0;
            else smooth[l1+l2*num_labels] = 1;
        }
    }

    double falloff_binary = 1.0; // smaller value -> more turning points

	try{
		GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph(num_elem,num_labels);
		gc->setDataCost(data);
		gc->setSmoothCost(smooth);
		
        for (int i=0; i<bnd.size() - 2; i++){
            Eigen::RowVector3d edge1 = (V.row(bnd[i+1]) - V.row(bnd[i])).normalized();
            Eigen::RowVector3d edge2 = (V.row(bnd[i+2]) - V.row(bnd[i+1])).normalized();
            edge1 = edge1.normalized();
            edge2 = edge2.normalized();
            double dot = (edge1.dot(edge2)-1)/falloff_binary;
            double cost = std::exp(-(1./2.0)*std::pow(dot,2));
            gc->setNeighbors(i, i+1, (int) (100*cost));
        }

		gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		
		for ( int  i = 0; i < num_elem; i++ )
			result[i] = gc->whatLabel(i);

		delete gc;
	}
	catch (GCException e){
		e.Report();
	}

	delete [] smooth;
	delete [] data;

    std::vector<int> turning_points;
    for (int i=1; i<result.rows(); i++){
        if (result(i) != result(i-1)){
            turning_points.push_back(i);
        }
    }

    // handle loops
    int last_id = bnd.size() - 1;
    if (bnd[0] == bnd[last_id] && result(0) != result(last_id - 1)){
        turning_points.push_back(0);
    }

    return turning_points;    
}