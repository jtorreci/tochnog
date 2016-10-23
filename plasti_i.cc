/*
	This routine enables to implement kinematic hardening, bounding
	surface plasticity models into Tochnog. 
	Source code is fully compatible with single element program Triax.
	
	Uses forward euler scheme with error controll and automatic substepping.
	Initial intersection with the yield surface found using nephton-raphson
		iterations for scalar multiplier alpha.
	
	options_element_dof were added to store dofs in integration points, rather 
	than only in nodes. 

	These changes into standard distribution of Tochnog are necessary:
	-include this file in makefile
	-add functions matrix_inverse_general, make_dev, matrix4_ab, matrix_a_contr_b
		matrix_ab4 and cinv into math.cc and declare these functions in tochnog.h
	-include matrix.h in tochnog.h
	-declare plasti_incr in tochnog.h
 	-cut-off backward--euler iterations in stress.cc 
	-add new models into database.cc
	-add new elastic models into stress.cc
	It is possible to use constitutive models source code of single element
	program Triax after minor modification
	
	For options_element_dof these files had to be modified:
	-tochnog.h, database.cc, create.cc, top.cc, materi.cc, initia.cc, elem.cc
	Modifications are mostly labeled with comment '// added for options_element_dof'
	search for 'element_dof' as well. As comparison with the professional version
	discovered some bug in node_dof from free version, whereas the 
	options_element_dof implemented here gives the same results as the
	professional version, default was set to options_element_dof -yes.
	
	For nonlocal calculation try to grep -i for:
		options_nonlocal_softvar
		materi_plasti_softvar_local, materi_plasti_softvar_nonlocal,
		element_nonlocal_weight, element_nonlocal_ipoint, nonlocal_element_info
		find_local_softvar, nonlocal_first_set, find_nonlocal_weights,
	hope it is enough

	To print materi_plasti_incremental_substeps into output change input.cc, 
	database.cc, tochnog.h, initia.cc, general.cc, materi.cc, plasti_incr.cc. Search
	for the word 'substeps'
	
	groundflow_addtopressure introduced in order to model hydrostatic pore pressure distribution 
	above phreatic level. Added into tochnog.h, database.cc and groundfl.cc. Search for 'addtopressure'.

	Messages from SuperLU solver were suppressed in get_perm_c.c and dgssv.c. Labelled by
	comment 'commented out!' :-). (Just in my copy) 
	
		David Masin masin@natur.cuni.cz
		
    This program is free softvare; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free softvare Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.


    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free softvare Foundation 
    59 Temple Place, Suite 330, Boston, MA, 02111-1307, USA
*/

/**********************header***************************************************************/

#include "tochnog.h"

const double PI=3.14159265358979323846264338327950288;
const double EULER = 2.7182818284590452354;	/* e */
const double onethird=0.333333333333333333333333333333333333333333333333333333333333333333333333333 /*:-)*/;
const double LIMITF=1.e-5;
const int MAX_INTERSEC_ITER=40;

struct Variable {
	double  sig[MDIM*MDIM], new_sig[MDIM*MDIM],	 	
	geo_sig[MDIM*MDIM], d_epp[MDIM*MDIM], epp[MDIM*MDIM], d_ept[MDIM*MDIM];
	double delta_sig_trial[MDIM*MDIM];
	double df_trial[MDIM*MDIM];
	double dgdsig[MDIM*MDIM];
	double kron_delta[MDIM][MDIM];
	double lambda_mltp [MDIM*MDIM];
	double matrix_Cep[MDIM][MDIM][MDIM][MDIM];  
	double matrix_C[MDIM][MDIM][MDIM][MDIM];  
	long int element, usematrix;
	double f, g, softv_nonl, softv_l, ddum[1];
	bool plasti, nonloc;
	vector<double> elasti_data;
	vector<double> plasti_data;        					
	vector<double> hisv;
	vector<double> other_f;		

	Variable() {
		inicialize();
	}
	void inicialize();
	void recount(); 
};

struct F_euler_errorc_data {
	double  rotated_old_sig[MDIM*MDIM], inc_ept[MDIM*MDIM], init_dgdsig[MDIM*MDIM];	
	vector<double> old_hisv;	//input
	double  double_sig[MDIM*MDIM], double_epp[MDIM*MDIM], double_dgdsig[MDIM*MDIM];
	vector<double> double_hisv;	//output after one step
	double  single_sig[MDIM*MDIM], single_epp[MDIM*MDIM], single_dgdsig[MDIM*MDIM];
	vector<double> single_hisv;	//output after half step
	double total_inc_epp[MDIM*MDIM], df_init[9]; 
	double init_lambda[MDIM*MDIM], single_lambda[MDIM*MDIM], double_lambda[MDIM*MDIM];
	bool largeerror;
	double num_substeps;
	double max_substeps;
	double min_substeps;
	double sigerror;
	F_euler_errorc_data() {
		inicialize();
	}
	void inicialize();
};

enum Num_pl_model {
	camclay_incremental,
	tskh_ai_plasti,
	tskh_std_plasti,
};

enum Integration_scheme {
	f_euler_errorcont,
	forw_eul_subst,
	forward_eul,
};

struct Plasti_rule {         							
	double df[9];
	Plasti_rule() {};
	virtual void give_f(double sig[MDIM*MDIM], double epp[9], double d_epp[9],
		 double &f, double &g, bool recount_hisv) {};
	virtual void give_matrix(bool only_f, bool elpl_matrix) 
		{cout<<"Plasti doesn't work virtually ";};					
	virtual void make_matrix_Cep (double df[], double dg[], double H, bool elpl_matrix);
	virtual void inic_hisv(double hisv[], bool only_substep) {};
	virtual void erase_rec_hist() {};
};

struct Count {			
	Plasti_rule* plasti;
	Count(Plasti_rule* b) : plasti(b) {
	};
	virtual void count_elasti(bool trial);		
	virtual void count_plasti();	
	virtual bool controll_plasti();
	virtual void calc_results();
};

struct Forward_euler : Count {			
	Forward_euler(Plasti_rule* bn) : Count(bn) {};
};

struct Forw_eul_substep : Count {			
	Forw_eul_substep(Plasti_rule* bn) : Count(bn) {};
	virtual void calc_results();
	virtual void finish_substeps();
	virtual void apply_errorc();
};

struct F_euler_errorc : Forw_eul_substep {			
	F_euler_errorc(Plasti_rule* bn) : Forw_eul_substep(bn) {};
	void apply_errorc();
	void calc_substeps();
};

enum Tskh_mode {
	tskh_std,
	tskh_genai,	
};

struct Camclay_incremental : Plasti_rule {					
	Camclay_incremental () : Plasti_rule() {}; 
	void give_f(double sig[MDIM*MDIM], double epp[9], double d_epp[9],
		 double &f, double &g, bool recount_hisv);
	void give_matrix(bool only_f, bool elpl_matrix);
	void inic_hisv(double hisv[], bool only_substep);
};

struct Tskh_general_plasti : Plasti_rule {	
	Tskh_mode mode;
	double m, kappa, lambda, T, S, psi, N, m_flow, r_m, r_mfl, m_cmp, m_flow_K0,
		rolodeb, rolodefl, drolodenbdlodenb, drolodefldlodenb;
	double sig_b[9], sig_a[9], I_0, sig_0[9], v;
	Tskh_general_plasti (Tskh_mode m) : Plasti_rule(), mode(m){}; 
	void give_f(double sig[MDIM*MDIM], double epp[9], double d_epp[9],
		 double &f, double &g, bool recount_hisv);
	void give_matrix(bool only_f, bool elpl_matrix);
	void inic_hisv(double hisv[], bool only_substep);
	void erase_rec_hist();
	void calc_rolode(double lode, double beta, double &rolode, double &drolodedlode, 
		bool alsoderivative);
};

/**********************end header***********************************************/

vector<Plasti_rule*> vector_plasti;
vector<Count*> countv;
Integration_scheme ischem;

Camclay_incremental camclay_incr;
Tskh_general_plasti tskh_ai_pl(tskh_genai);
Tskh_general_plasti tskh_std_pl(tskh_std);

//accesible for everyone
Variable vari;	
F_euler_errorc_data errorc;	

void Variable :: inicialize() {
	array_set(sig, 0, MDIM*MDIM);
	array_set(new_sig, 0, MDIM*MDIM);
	array_set(d_epp, 0, MDIM*MDIM);
	array_set(d_ept, 0, MDIM*MDIM);
	array_set(epp	, 0, MDIM*MDIM);
	array_set(delta_sig_trial, 0, MDIM*MDIM);
	array_set(df_trial, 0, MDIM*MDIM);
	array_set(dgdsig, 0, MDIM*MDIM);
	array_set(lambda_mltp, 0, MDIM*MDIM);
	array_set(&matrix_C[0][0][0][0], 0, MDIM*MDIM*MDIM*MDIM);
	array_set(&matrix_Cep[0][0][0][0], 0, MDIM*MDIM*MDIM*MDIM);
	f=-10;
	g=-1.e10;
	softv_nonl=softv_l=ddum[0]=0;
	plasti=nonloc=false;
	element=0;
	usematrix=-NO;
	for(int i=0; i<MDIM; i++) {
		for(int j=0; j<MDIM; j++) {
			if(i==j) kron_delta[i][j]=1;
			else kron_delta[i][j]=0;
		}
	}
	recount();

	vector_plasti.push_back(&camclay_incr);
	vector_plasti.push_back(&tskh_ai_pl);
	vector_plasti.push_back(&tskh_std_pl);
}

void F_euler_errorc_data :: inicialize() {
	array_set(rotated_old_sig, 0, MDIM*MDIM);
	array_set(inc_ept, 0, MDIM*MDIM);
	array_set(double_sig, 0, MDIM*MDIM);
	array_set(double_epp, 0, MDIM*MDIM);
	array_set(single_sig, 0, MDIM*MDIM);
	array_set(single_epp, 0, MDIM*MDIM);
	array_set(init_dgdsig, 0, MDIM*MDIM);
	array_set(single_dgdsig, 0, MDIM*MDIM);
	array_set(double_dgdsig, 0, MDIM*MDIM);
	array_set(total_inc_epp, 0, MDIM*MDIM);
	array_set(df_init, 0, MDIM*MDIM);
	array_set(init_lambda, 0, MDIM*MDIM);
	array_set(single_lambda, 0, MDIM*MDIM);
	array_set(double_lambda, 0, MDIM*MDIM);
	largeerror=true;
	num_substeps=1;
	max_substeps=5000000;
	min_substeps=2;
	sigerror=0.00001;
}

void Variable::recount() {
	array_multiply(sig, geo_sig, -1, MDIM*MDIM);
}


void plasti_incr( double rotated_old_sig[], double inc_ept[], double old_epp[], 
	    double new_ept[], double C[3][3][3][3], double old_hisv[], long int plasti_type, 
	    double plasti_dt[], long lenght_pl, long int gr, long int element, 
	    double softvar_nonl, double &softvar_l,
            double new_sigv[], double inc_epp[], double new_hisv[], double &new_f, double &new_substeps,
	    double Cep_cons[3][3][3][3]) {

	/*******************************************************
	This function is an iterface between single element program
	Triax (D.Masin@city.ac.uk) and Tochnog. 
	*******************************************************/
	Num_pl_model num_plasti=tskh_std_plasti;
	double elasti_dt[DATA_ITEM_SIZE], ddum[1]={0};
	for(int i=0; i<DATA_ITEM_SIZE; i++) elasti_dt[i]=0;
	long lenght_el=0;
	long int ldum=0, idum[1]={0};

	if(plasti_type==GROUP_MATERI_PLASTI_TSKH) {
		num_plasti=tskh_std_plasti;
	  	if(!get_group_data( GROUP_MATERI_ELASTI_TSKH, gr, element, ddum, 
		      elasti_dt, lenght_el, GET_IF_EXISTS )) {
		      cout<<"3-SKH model must be used with the 3-SKH elasticity"<<endl;
		      exit(1);
		}
	}
	else if(plasti_type==GROUP_MATERI_PLASTI_AITSKH) {
		num_plasti=tskh_ai_plasti;
	  	if(!get_group_data( GROUP_MATERI_ELASTI_TSKH, gr, element, ddum, 
		      elasti_dt, lenght_el, GET_IF_EXISTS )) {
		      cout<<"AI3-SKH model must be used with the 3-SKH elasticity"<<endl;
		      exit(1);
		}
	}
	else if(plasti_type==GROUP_MATERI_PLASTI_CAMCLAY_INCREMENTAL) {
		num_plasti=camclay_incremental;
	}
	else {
		cout<<"New model was not added into plasti_incr"<<endl;
		exit(1);
	}
	F_euler_errorc f_euler_errorc(vector_plasti[num_plasti]);
	Forw_eul_substep forw_eul_substep(vector_plasti[num_plasti]);
	Forward_euler forward_euler(vector_plasti[num_plasti]);
	countv.push_back(&f_euler_errorc);
	countv.push_back(&forw_eul_substep);
	countv.push_back(&forward_euler);

	//ischem=forward_eul;		// not in use actually - without substeps
	ischem=f_euler_errorcont;	// default with error 0.00001
	
	if(get_group_data( GROUP_MATERI_PLASTI_INCREMENTAL_FESUBSTEPS, gr, element, ddum, 
		      &errorc.num_substeps, ldum, GET_IF_EXISTS )) {
		ischem=forw_eul_subst;		      
	}
	else if(get_group_data( GROUP_MATERI_PLASTI_INCREMENTAL_FEERROR, gr, element, ddum, 
		      &errorc.sigerror, ldum, GET_IF_EXISTS )) {
		ischem=f_euler_errorcont;
	}
	if(find_local_softvar) {
		ischem=forw_eul_subst;		      
		errorc.num_substeps=1;
	}
	
	get_group_data( GROUP_MATERI_PLASTI_INCREMENTAL_MAXSUBSTEPS, gr, element, ddum, 
		      &errorc.max_substeps, ldum, GET_IF_EXISTS );
	get_group_data( GROUP_MATERI_PLASTI_INCREMENTAL_MINSUBSTEPS, gr, element, ddum, 
		      &errorc.min_substeps, ldum, GET_IF_EXISTS );
	if (db_active_index( GROUP_MATERI_PLASTI_INCREMENTAL_USEMATRIX, gr, VERSION_NORMAL) ) 
	      db( GROUP_MATERI_PLASTI_INCREMENTAL_USEMATRIX, gr, &vari.usematrix, ddum, ldum, 
	        VERSION_NORMAL, GET );
	if (scalar_dabs(options_nonlocal_softvar)>TINY) {
		//vari.usematrix=-NO;
		if(!find_local_softvar) vari.nonloc=true;
		else vari.nonloc=false;
	}

		//Feeding vari with current values
	for(int i=0; i<lenght_el; i++) vari.elasti_data.push_back(elasti_dt[i]);
	for(int i=0; i<lenght_pl; i++) vari.plasti_data.push_back(plasti_dt[i]);
	array_move(rotated_old_sig, vari.sig, MDIM*MDIM);
	array_move(rotated_old_sig, vari.new_sig, MDIM*MDIM);
	array_move(inc_ept, vari.d_ept, MDIM*MDIM);
	array_move(&C[0][0][0][0], &vari.matrix_C[0][0][0][0], MDIM*MDIM*MDIM*MDIM);
	countv[ischem]->plasti->inic_hisv(old_hisv, false);
	vari.element=element;
	vari.softv_nonl=softvar_nonl;
	vari.softv_l=softvar_l;
	vari.recount();

		//Erase the influence of recent history
        if ( db_active_index( GROUP_MATERI_PLASTI_INCREMENTAL_ERASERECENTHISTORY, 
             gr, VERSION_NORMAL ) ) {
	     	long int eraserecenthistory=-NO;
		double time[1]={1};
		db( GROUP_MATERI_PLASTI_INCREMENTAL_ERASERECENTHISTORY, gr,
        		&eraserecenthistory, ddum, ldum, VERSION_NORMAL, GET );
	        db( TIME_CURRENT, 0, idum, time, ldum, VERSION_NORMAL, GET );
		if(time[0]==0 && eraserecenthistory==-YES) countv[ischem]->plasti->erase_rec_hist();
        }

		//calculation

	countv[ischem]->calc_results();
	
	long int printf;
	if (db_active_index( GROUP_MATERI_PLASTI_INCREMENTAL_PRINTF, gr, VERSION_NORMAL) ) {
	      db( GROUP_MATERI_PLASTI_INCREMENTAL_PRINTF, gr, &printf, ddum, ldum, 
	        VERSION_NORMAL, GET );
	      if ( printf==element && !find_local_softvar) 
	      	cout<<vari.f<<" "<<vari.other_f[0]<<" "<<vari.other_f[1]<<endl;
	}

		//Feeding output
        array_move(vari.new_sig, new_sigv, MDIM*MDIM);
        array_move(vari.d_epp, inc_epp, MDIM*MDIM);
	
	double d_sig[MDIM*MDIM];
	array_set(d_sig, 0, MDIM*MDIM);
	array_subtract(vari.new_sig, rotated_old_sig, d_sig, MDIM*MDIM);

	for(int i=0; i<int(vari.hisv.size()); i++)
		new_hisv[i]=vari.hisv[i];	 
	new_f=vari.f;
	
	if(new_f>0) new_substeps=errorc.num_substeps;
	else new_substeps=1;
	
	array_move(&vari.matrix_Cep[0][0][0][0], &Cep_cons[0][0][0][0], MDIM*MDIM*MDIM*MDIM);

		//cleaning vari for the next step
	vari.elasti_data.clear();
	vari.plasti_data.clear();        					
	vari.hisv.clear();
	vari.other_f.clear();		
	countv.clear();
}

void Count::calc_results() {	
	if(controll_plasti()==true) count_plasti();
	else count_elasti(false);

	//Calculates new history variables
	if(vari.plasti) plasti->give_f(vari.sig, vari.epp, vari.d_epp, vari.ddum[0], vari.ddum[0], true);	
}

bool Count::controll_plasti() {
	double dfdsig;	
	plasti->give_f(vari.sig, vari.epp, vari.d_epp, vari.f, vari.g, false);
	count_elasti(true);		//calculates trial stress rate
	if (vari.f>=0) {
		plasti->give_matrix(false, false);	//calculates df/dsig
		dfdsig=array_inproduct(vari.df_trial, vari.delta_sig_trial, MDIM*MDIM);	
		if( dfdsig > 0) vari.plasti=true;
		else vari.plasti=false;
	}
	else vari.plasti = false;
	return(vari.plasti);
}

void Count::count_elasti(bool trial) {	
	double delta_epe[MDIM*MDIM], delta_sig[MDIM*MDIM];
	array_set(delta_epe, 0, MDIM*MDIM);
	array_set(delta_sig, 0, MDIM*MDIM);
	array_move(vari.d_ept, delta_epe, MDIM*MDIM);
	matrix_a4b(vari.matrix_C, delta_epe, delta_sig);

	if(trial) array_move(delta_sig, vari.delta_sig_trial, MDIM*MDIM);
	else array_add(vari.sig, delta_sig, vari.new_sig, MDIM*MDIM);

	array_set(vari.d_epp, 0, MDIM*MDIM);
}	

void Count::count_plasti() {			
	array_move(vari.sig, vari.new_sig, MDIM*MDIM);
	double delta_ept[MDIM*MDIM];
	array_set(delta_ept, 0, MDIM*MDIM);

	double delta_sig[MDIM*MDIM];
	array_set(delta_sig, 0, MDIM*MDIM);
	double delta_epe[MDIM*MDIM];
	array_set(delta_epe, 0, MDIM*MDIM);

	double delta_epp[MDIM*MDIM];
	array_set(delta_epp, 0, MDIM*MDIM);

	array_move(vari.d_ept, delta_ept, MDIM*MDIM);

	double used_lambda=array_inproduct(vari.lambda_mltp, vari.d_ept, MDIM*MDIM);
	array_multiply(vari.dgdsig, delta_epp, used_lambda, MDIM*MDIM);
	array_subtract(delta_ept, delta_epp, delta_epe, MDIM*MDIM);
	matrix_a4b(vari.matrix_C, delta_epe, delta_sig);

	array_add(vari.sig, delta_sig, vari.new_sig, MDIM*MDIM);
	array_move(delta_epp, vari.d_epp, 9);
}

/******Forward Euler, const given number of substeps*************/

void Forw_eul_substep::calc_results() {	
	array_set(errorc.total_inc_epp, 0, MDIM*MDIM);

	controll_plasti();

	if(!vari.plasti) {
		double testf_after = -10;
		count_elasti(false);		//is stress outside after elastic step?
		plasti->give_f(vari.new_sig, vari.epp, vari.d_epp, testf_after, vari.ddum[0], false);	

			//if yes -> calculates initial intersection with the yield surface
		if(testf_after>0 && vari.f<0) {
			double delta_sig_e[MDIM*MDIM];
			double sig_init[MDIM*MDIM];
			array_subtract(vari.new_sig, vari.sig, delta_sig_e, MDIM*MDIM);
			array_move(vari.sig, sig_init, MDIM*MDIM);
			double alpha=-vari.f/(testf_after-vari.f);
			double dalpha=0, dftrdsige=0;
			double num_iter=1;
			double dsige_times_alpha[MDIM*MDIM];
			double df_transposed[MDIM*MDIM];
			array_set(dsige_times_alpha, 0, MDIM*MDIM);
			array_set(df_transposed, 0, MDIM*MDIM);
			
				//newton-raphson iterations for scalar multiplier alpha			
			while(1) {
				array_multiply(delta_sig_e, dsige_times_alpha, alpha, MDIM*MDIM);
				array_add(sig_init, dsige_times_alpha, vari.sig, MDIM*MDIM);
				vari.recount();
				plasti->give_f(vari.sig, vari.epp, vari.d_epp, vari.f, vari.ddum[0], false);	
				plasti->give_matrix(true, false);	//calculates df/dsig_i

				if(num_iter>MAX_INTERSEC_ITER || (vari.f>0 && vari.f<LIMITF)) break;
				
				for(int i=0; i<MDIM; i++) {	//transpose df/dsig_i
					for(int j=0; j<MDIM; j++) {
						df_transposed[3*i+j]=vari.df_trial[3*j+i];
					}
				}
				dftrdsige=array_inproduct(df_transposed, delta_sig_e, MDIM*MDIM);
				if(!dftrdsige==0) dalpha=-vari.f/dftrdsige;
				else dalpha=0;
				num_iter++;
				alpha+=dalpha;
			}
				//new d_ept for el-pl time integration, vari.sig is now on yield surface 
			array_multiply(vari.d_ept, vari.d_ept, (1-alpha), MDIM*MDIM);	
			vari.plasti=true;
		}
	}
	if(vari.plasti) {
		apply_errorc();
		array_set(errorc.total_inc_epp, 0, MDIM*MDIM);		
		array_multiply(vari.d_ept, vari.d_ept, (1/errorc.num_substeps), MDIM*MDIM);

		for(int i=0; i<errorc.num_substeps; i++) 
				finish_substeps();

		array_move(vari.new_sig, vari.sig, MDIM*MDIM);
	}

	if(vari.plasti && vari.usematrix==-YES) plasti->give_matrix(false, true);	
	else array_move(&vari.matrix_C[0][0][0][0], &vari.matrix_Cep[0][0][0][0], MDIM*MDIM*MDIM*MDIM);
	array_move(errorc.total_inc_epp, vari.d_epp, MDIM*MDIM);
}

void Forw_eul_substep::finish_substeps() {	
	vari.recount();
	plasti->give_f(vari.sig, vari.epp, vari.d_epp, vari.f, vari.g, false);	//before step->calc f
	plasti->give_matrix(false, false);	
	count_plasti();
	plasti->give_f(vari.sig, vari.epp, vari.d_epp, vari.ddum[0], vari.ddum[0], true);	//after step-> calc hisv
	array_move(vari.new_sig, vari.sig, MDIM*MDIM);
	array_add(errorc.total_inc_epp, vari.d_epp, errorc.total_inc_epp, MDIM*MDIM);
}

void Forw_eul_substep::apply_errorc() {	
}

/******************Forward Euler, control of local error, automatic substepping*********************/

void F_euler_errorc::apply_errorc() {	
	errorc.largeerror=true;
	errorc.num_substeps=2;

		//initial values
	array_move(vari.sig, errorc.rotated_old_sig, MDIM*MDIM);
	array_move(vari.d_ept, errorc.inc_ept, MDIM*MDIM);
	for(int i=0; i<int(vari.hisv.size()); i++) errorc.old_hisv.push_back(vari.hisv[i]);
	for(int i=0; i<int(vari.hisv.size()); i++) errorc.single_hisv.push_back(vari.hisv[i]);
	for(int i=0; i<int(vari.hisv.size()); i++) errorc.double_hisv.push_back(vari.hisv[i]);
	plasti->give_matrix(false, false);	
	array_move(vari.dgdsig, errorc.init_dgdsig, MDIM*MDIM);
	array_move(vari.lambda_mltp, errorc.init_lambda, MDIM*MDIM);
       	array_move(plasti->df, errorc.df_init, MDIM*MDIM);
	
	count_plasti();
	plasti->give_f(vari.sig, vari.epp, vari.d_epp, vari.f, vari.ddum[0], true);	//new_hisv

		//saves results after one step
       	array_move(vari.new_sig, errorc.double_sig, MDIM*MDIM);
       	array_move(vari.d_epp, errorc.double_epp, MDIM*MDIM);
	for(int i=0; i<int(vari.hisv.size()); i++) errorc.double_hisv[i]=vari.hisv[i];

	while(errorc.largeerror) calc_substeps();

		//refeed with init. values
	array_move(errorc.rotated_old_sig, vari.sig, MDIM*MDIM);
	array_move(errorc.inc_ept, vari.d_ept, MDIM*MDIM);
	array_set(errorc.total_inc_epp, 0, MDIM*MDIM);
	for(int i=0; i<int(vari.hisv.size()); i++) vari.hisv[i]=errorc.old_hisv[i];
	plasti->inic_hisv(vari.ddum, true);

	errorc.old_hisv.clear();
	errorc.double_hisv.clear();
	errorc.single_hisv.clear();
}

void F_euler_errorc::calc_substeps() {	
	
		//refeed vari with init values
	array_move(errorc.rotated_old_sig, vari.sig, MDIM*MDIM);
	array_move(errorc.df_init, plasti->df, MDIM*MDIM);
	array_move(errorc.init_dgdsig, vari.dgdsig, MDIM*MDIM);
	array_move(errorc.init_lambda, vari.lambda_mltp, MDIM*MDIM);
	array_multiply(vari.d_ept, vari.d_ept, 0.5, MDIM*MDIM);	//half step size
	for(int i=0; i<int(vari.hisv.size()); i++) vari.hisv[i]=errorc.old_hisv[i];
	plasti->inic_hisv(vari.ddum, true);
	
		//calc first substep
	plasti->give_f(vari.sig, vari.epp, vari.d_epp, vari.f, vari.g, false);
	count_plasti();
	plasti->give_f(vari.sig, vari.epp, vari.d_epp, vari.ddum[0], vari.ddum[0], true);	
	array_move(vari.new_sig, vari.sig, MDIM*MDIM);

		//saves results after substep
        array_move(vari.new_sig, errorc.single_sig, MDIM*MDIM);
       	array_move(vari.d_epp, errorc.single_epp, MDIM*MDIM);
	for(int i=0; i<int(vari.hisv.size()); i++) errorc.single_hisv[i]=vari.hisv[i];
	vari.recount();
	plasti->give_matrix(false, false);	
	array_move(vari.dgdsig, errorc.single_dgdsig, MDIM*MDIM);
	array_move(vari.lambda_mltp, errorc.single_lambda, MDIM*MDIM);

		//calc second substep
	plasti->give_f(vari.sig, vari.epp, vari.d_epp, vari.f, vari.ddum[0], false);	
	count_plasti();
	plasti->give_f(vari.sig, vari.epp, vari.d_epp, vari.ddum[0], vari.ddum[0], true);	
	array_move(vari.new_sig, vari.sig, MDIM*MDIM);
	array_add(errorc.total_inc_epp, vari.d_epp, errorc.total_inc_epp, MDIM*MDIM);
	
		//calc local trucation error
	double errorsig[MDIM*MDIM];
	array_set(errorsig, 0, MDIM*MDIM);
	array_subtract(errorc.double_sig, vari.new_sig, errorsig, MDIM*MDIM);
	double error = sqrt(array_inproduct(errorsig, errorsig, MDIM*MDIM));
	
	if((error<errorc.sigerror|| errorc.num_substeps>errorc.max_substeps)&&
		(errorc.num_substeps>(errorc.min_substeps-0.1))) errorc.largeerror=false;
	else {
		errorc.largeerror=true;

			//half step is main step for next turn
        	array_move(errorc.single_sig, errorc.double_sig, MDIM*MDIM);
	       	array_move(errorc.single_epp, errorc.double_epp, MDIM*MDIM);
		for(int i=0; i<int(vari.hisv.size()); i++) errorc.double_hisv[i]=errorc.single_hisv[i];
        	array_move(errorc.single_dgdsig, errorc.double_dgdsig, MDIM*MDIM);
	       	array_move(errorc.single_lambda, errorc.double_lambda, MDIM*MDIM);

		array_set(errorc.total_inc_epp, 0, MDIM*MDIM);
		errorc.num_substeps*=2;
	}
}

/********************Elasto-plastic stiffness matrix***********************************/

void Plasti_rule :: make_matrix_Cep (double df[], double dg[], double H, 
		bool elpl_matrix) { 	

	double tmp1[MDIM*MDIM], tmp2[MDIM*MDIM], tmp3[MDIM][MDIM][MDIM][MDIM], chi=1;
	double Cp[MDIM][MDIM][MDIM][MDIM];

	matrix_a4b(vari.matrix_C, dg, tmp1);
	matrix_ab4(df, vari.matrix_C, tmp2);
	matrix4_ab(tmp1, tmp2, tmp3);
	matrix_a_contr_b(df, tmp1, chi);

	chi=chi+H;
	if(chi==0) {
		cout<<" Problem in calculation of elasto--plastic stiffness matrix -> chi=0 ";
		//exit_tn( -YES ); 
	}
	array_multiply(tmp2, vari.lambda_mltp, 1/chi, MDIM*MDIM);
	
	if(elpl_matrix) {
	  for (int i=0; i<3; i++ ) {
	    for (int j=0; j<3; j++ ) {
	      for (int k=0; k<3; k++ ) {
        	for (int l=0; l<3; l++ ) {
	          Cp[i][j][k][l] = tmp3[i][j][k][l]/chi;
        	}
	      }
	    }
	  }

	  for (int i=0; i<3; i++ ) {
	    for (int j=0; j<3; j++ ) {
	      for (int k=0; k<3; k++ ) {
        	for (int l=0; l<3; l++ ) {
	          vari.matrix_Cep[i][j][k][l] = vari.matrix_C[i][j][k][l] - Cp[i][j][k][l];
        	}
	      }
	    }
	  }
	}  
	array_move(dg, vari.dgdsig, MDIM*MDIM);
}

/***************************MODELS***********************************************/
/******** Cam-Clay -- incremental solution for comparison with iterative ********/

void Camclay_incremental::inic_hisv(double hisv[], bool only_substep) {
	array_set(df, 0, 9);
	vari.hisv.push_back(hisv[0]);
	vari.hisv.push_back(hisv[1]);
}

void Camclay_incremental :: give_f(double sig[MDIM*MDIM], double epp[9], double d_epp[9],
		 double &f, double &g, bool recount_hisv) {
		 
	double geo_sigd[9], d_eppd[9];
	array_set(geo_sigd, 0, 9);
	array_set(d_eppd, 0, 9);
	array_multiply(sig, geo_sigd, -1, MDIM*MDIM);
	array_multiply(d_epp, d_eppd, -1, MDIM*MDIM);
	double m = vari.plasti_data[0];
	double kappa = vari.plasti_data[1];
	double lambda = vari.plasti_data[2];
	double N = vari.plasti_data[3];
	if ( scalar_dabs(lambda-kappa)==0. ) {
		cout<<"Error in Cam--Clay: lambda==kappa";
		exit(1);
	}
	double dev = vari.d_ept[0]+vari.d_ept[4]+vari.d_ept[8];
	double p = (geo_sigd[0]+geo_sigd[4]+geo_sigd[8])/3;
	double e = vari.hisv[0];
	//double p0 = vari.hisv[1];
	double p0 = exp((N-kappa*log(p)-(1+e))/(lambda-kappa));
	double devp = (d_eppd[0]+d_eppd[4]+d_eppd[8]);
	double de = dev*(1.+e);
	double dp0 = devp * p0 * (1.+e)/(lambda-kappa);
	e += de;
	p0 += dp0;
	double q = sqrt( 
        	0.5*( scalar_square(geo_sigd[0]-geo_sigd[4])+
         	     scalar_square(geo_sigd[4]-geo_sigd[8])+
	              scalar_square(geo_sigd[0]-geo_sigd[8]) ) + 
        	3.*( scalar_square(geo_sigd[1]) +
	             scalar_square(geo_sigd[2]) +
        	     scalar_square(geo_sigd[5]) ) );
	f = q*q - m*m*(p*(p0-p));
	g = q*q - m*m*(p*(p0-p));	

	if(recount_hisv) {	
	        vari.hisv[0] += de;
        	if (vari.f>0) vari.hisv[1] = p0;
	}
}

void Camclay_incremental::give_matrix(bool only_f, bool elpl_matrix) {
	double m = vari.plasti_data[0];
	double kappa = vari.plasti_data[1];
	double lambda = vari.plasti_data[2];
	double N = vari.plasti_data[3];
	double e = vari.hisv[0];	
	double p = (vari.geo_sig[0]+vari.geo_sig[4]+vari.geo_sig[8])/3;
	//double p0 = vari.hisv[1];
	double p0 = exp((N-kappa*log(p)-(1+e))/(lambda-kappa));

	double df[MDIM*MDIM], dg[MDIM*MDIM], sigd[MDIM*MDIM], zero[MDIM*MDIM], eppd[MDIM*MDIM];
	array_set(zero, 0, MDIM*MDIM);
	array_set(dg, 0, MDIM*MDIM);
	array_set(df, 0, MDIM*MDIM);
	array_set(sigd, 0, MDIM*MDIM);
	array_set(eppd, 0, MDIM*MDIM);
	array_move(vari.sig, sigd, MDIM*MDIM);
	double f_right=0, f_left=0;

	double signorm=array_size(vari.sig, MDIM*MDIM)/100000;	
	if(signorm==0) signorm=0.0000001;
	for(int i=0; i<MDIM*MDIM; i++) {
		sigd[i] += signorm;
		give_f(sigd, zero, zero, f_right, f_right, false);
		sigd[i] -= 2*signorm;
		give_f(sigd, zero, zero, f_left, f_left, false);
		df[i]=(f_right-f_left)/(2*signorm);
		array_move(vari.sig, sigd, MDIM*MDIM);
	}
	array_move(df, dg, MDIM*MDIM);

	double dfdep[MDIM*MDIM];
	array_set(dfdep, 0, MDIM*MDIM);

	dfdep[0]=dfdep[4]=dfdep[8]=(1+e)*m*m*p*p0/(lambda-kappa);

	double H=0;
	for(int i=0; i<MDIM*MDIM; i++) {
		H -= dfdep[i]*dg[i];
	}
	array_move(df, vari.df_trial, MDIM*MDIM);
	if(!only_f) make_matrix_Cep(df, dg, H, elpl_matrix);
}

/***************************3-SKH********************************************************/

void Tskh_general_plasti::erase_rec_hist() {
	vari.hisv[2]=vari.hisv[8]=-vari.sig[0];
	vari.hisv[3]=vari.hisv[9]=-vari.sig[1];
	vari.hisv[4]=vari.hisv[10]=-vari.sig[2];
	vari.hisv[5]=vari.hisv[11]=-vari.sig[4];
	vari.hisv[6]=vari.hisv[12]=-vari.sig[5];
	vari.hisv[7]=vari.hisv[13]=-vari.sig[8];

	inic_hisv(vari.ddum, true);
}

void Tskh_general_plasti::inic_hisv(double hisv[], bool only_substep) {

	if(!only_substep) {
		m_cmp = vari.plasti_data[0];
		kappa = vari.plasti_data[1];
		lambda = vari.plasti_data[2];
		T = vari.plasti_data[3];
		S = vari.plasti_data[4];
		psi = vari.plasti_data[5];
		N = vari.plasti_data[6];
		m_flow = m = m_flow_K0 = m_cmp;

		r_m=r_mfl=3/(3+m_cmp);
		
		array_set(df, 0, 9);
		array_set(sig_b, 0, 9);
		array_set(sig_a, 0, 9);
		array_set(sig_0, 0, 9);
	
		for(int i=0; i<14; i++)	vari.hisv.push_back(hisv[i]);

		vari.other_f.push_back(-10);
		vari.other_f.push_back(-10);
	}

	v=vari.hisv[0]+1;  //e stored as a history variable
	I_0=3*vari.hisv[1];
	sig_0[0]=sig_0[4]=sig_0[8]=I_0/3;
	sig_a[0]=vari.hisv[2];
	sig_a[1]=sig_a[3]=vari.hisv[3];
	sig_a[2]=sig_a[6]=vari.hisv[4];
	sig_a[4]=vari.hisv[5];
	sig_a[5]=sig_a[7]=vari.hisv[6];
	sig_a[8]=vari.hisv[7];
	sig_b[0]=vari.hisv[8];
	sig_b[1]=sig_b[3]=vari.hisv[9];
	sig_b[2]=sig_b[6]=vari.hisv[10];
	sig_b[4]=vari.hisv[11];
	sig_b[5]=sig_b[7]=vari.hisv[12];
	sig_b[8]=vari.hisv[13];
}

void Tskh_general_plasti :: give_f(double sig[MDIM*MDIM], double epp[9], double d_epp[9],
		 double &f, double &g, bool recount_hisv) {
		 
	double geo_sig[9], d_geo_epp[9];
	array_set(geo_sig, 0, 9);
	array_set(d_geo_epp, 0, 9);
	array_multiply(sig, geo_sig, -1, MDIM*MDIM);
	array_multiply(d_epp, d_geo_epp, -1, MDIM*MDIM);
	double f_hist=vari.other_f[0]; 
	double f_bound=vari.other_f[1]; 
	double f_yield=vari.f;
	double dsig_b[9], dsig_a[9];
	if ( scalar_dabs(lambda-kappa)==0. ) {
		cout<<"Problem in the 3-SKH model: lambda must differ from kappa";
		exit(1);
	}

	double p = (geo_sig[0]+geo_sig[4]+geo_sig[8])/3;
	if(p<=0) p=TINY;
	double A = vari.elasti_data[0];
	double n= vari.elasti_data[1];
	double m_elasti = vari.elasti_data[2];
	double R0=2*I_0 / (3*p);
	double G = A*scalar_power(p, n)*scalar_power(R0, m_elasti);
	double K = p/kappa;
	double poisson = (3*K - 2*G)/(2*G + 6*K);
	if(poisson>=0.5) {
		cout<<"error in AI3-SKH, location 1, element "<<vari.element<<endl;
		//exit_tn( -YES );
	}
	
	//m_cmp = vari.plasti_data[0];

	double mfl_dnm=
		(((1-2*poisson)*(6-m_cmp)*lambda-m_cmp*kappa*(1+poisson))*((6-m_cmp)*(6-m_cmp)-9));
	if(mfl_dnm<=0) {
		cout<<"error in AI3-SKH, location 2, element "<<vari.element<<endl;
		//exit_tn( -YES );
	}
	m_flow_K0=3*(6-m_cmp)*sqrt((1-2*poisson)*(lambda-kappa)*m_cmp/mfl_dnm);
	if(m_flow_K0<=0)  {
		cout<<"error in AI3-SKH, location 3, element "<<vari.element<<endl;
		//exit_tn( -YES );
	}

	if(mode==tskh_std) {
		r_m=1;
		r_mfl=1;
		m_flow_K0=m_cmp;
	}
	
	double sig_sigb[9];
	double sig_siga[9];
	array_set(sig_sigb, 0, 9);
	array_set(sig_siga, 0, 9);
	array_subtract(geo_sig, sig_b, sig_sigb, 9);
	array_subtract(geo_sig, sig_a, sig_siga, 9);
	
	if(mode==tskh_genai) {
		r_mfl=m_cmp/m_flow_K0;
		// In 2D y axis is vertical
		// in 3D z axis is vertical
		if(ndim==2) {
			if(!((sig_sigb[4]>sig_sigb[8])&&(sig_sigb[4]>sig_sigb[0]))) {
				r_mfl=1;
				m_flow_K0=m_cmp;
			}
		}
		else if(ndim==3) {
			if(!((sig_sigb[8]>sig_sigb[4])&&(sig_sigb[8]>sig_sigb[0]))) {
				r_mfl=1;
				m_flow_K0=m_cmp;
			}
		}
		else {
			cout<<"Can not use AI3-SKH in 1D"<<endl;
			exit_tn( -YES );
		}
	}
	/*if(mode==tskh_genai) {	//MN3-SKH
		r_mfl=1;
		m_flow_K0=m_cmp;
	}*/

	if(recount_hisv && vari.plasti) {	

		/************* isotropic hardening ******************/

		double devp = (d_geo_epp[0]+d_geo_epp[4]+d_geo_epp[8]);
		double dI0 = devp * I_0/(lambda-kappa); 	

		if(I_0<=0) {
			cout<<"error in AI3-SKH, location 4, element "<<vari.element<<endl;
			//exit_tn( -YES );
		}

		array_set(dsig_b, 0, 9);
		array_set(dsig_a, 0, 9);
		for (int j=0; j<9; j++ ) {
			dsig_b[j]=dI0*sig_b[j]/I_0;
			dsig_a[j]=dI0*sig_a[j]/I_0;
		}
		for (int j=0; j<9; j++ ) {
			sig_b[j]+=dsig_b[j];
			sig_a[j]+=dsig_a[j];
		}
		I_0+=dI0;
		sig_0[0]=sig_0[4]=sig_0[8]=I_0/3;

		/****************** direction of translation ******************/

		double beta[9];		//direction for history surface
		array_set(beta,0,9);
		for (int j=0; j<9; j++ ) {
			beta[j]=(vari.geo_sig[j]-sig_a[j])/T+sig_0[j]-vari.geo_sig[j];
		}
		double gama[9];		//direction for yield surface
		array_set(gama,0,9);
		for (int j=0; j<9; j++ ) {
			gama[j]=(vari.geo_sig[j]-sig_b[j])/S+sig_a[j]-vari.geo_sig[j];
		}

		/************* kinematic hardening ******************/

		array_set(dsig_b, 0, 9);
		array_set(dsig_a, 0, 9);
		double d_geosig[9];
		array_set(d_geosig, 0, 9);
		array_subtract(vari.new_sig, vari.sig, d_geosig, MDIM*MDIM);
		array_multiply(d_geosig, d_geosig, -1, 9);
		double tmpup=0, tmpdown=0;

		if (f_yield>=0 && f_hist<0 && f_bound<0) {	
			tmpup=0;
			tmpdown=0;
			double dfdI0=-T*T*S*S*I_0/9;
			for (int j=0; j<9; j++ ) {
				tmpup +=-df[j]*d_geosig[j]+df[j]*dI0*sig_b[j]/I_0;
				tmpdown +=-df[j]*gama[j];	//-df <- back to geotechnical convention
			}
			if(tmpdown==0) {
				cout<<"error in AI3-SKH, location 5, element "<<vari.element<<endl;
				//exit_tn( -YES );
			}

			double Zs=(tmpup+dfdI0*dI0)/tmpdown;
			for (int j=0; j<9; j++ ) {
				dsig_b[j]=Zs*gama[j];
			}
			for (int j=0; j<9; j++ ) {
				sig_b[j]+=dsig_b[j];
			}
		}
		else if (f_hist>=0 && f_bound < 0 && f_yield>=0) {	
			tmpup=0;
			tmpdown=0;
			double dfhdI0=-T*T*I_0/9;
			double dfhdsig[9];
			array_set(dfhdsig, 0, 9);
			for (int j=0; j<9; j++ ) {
				dfhdsig[j]=df[j]/S;
			}
			for (int j=0; j<9; j++ ) {
				tmpup +=-dfhdsig[j]*d_geosig[j]+dfhdsig[j]*dI0*sig_a[j]/I_0;
				tmpdown +=-dfhdsig[j]*beta[j];	//-df <- back to geotechnical convention
			}
			if(tmpdown==0) {
				cout<<"error in AI3-SKH, location 6, element "<<vari.element<<endl;
				//exit_tn( -YES );
			}
			double Ws=(tmpup+dfhdI0*dI0)/tmpdown;
			for (int j=0; j<9; j++ ) {
				dsig_a[j]=Ws*beta[j];
			}
			for (int j=0; j<9; j++ ) {
				sig_a[j]+=dsig_a[j];
				sig_b[j]=geo_sig[j]-(geo_sig[j]-sig_a[j])*S;
			}
		}
		else if(f_bound >= 0) {
			//A position of the yield and hist. surf. on bounding s., touching at conj. points			
			for (int j=0; j<9; j++ ) {
				sig_b[j]=geo_sig[j]-(geo_sig[j]-sig_0[j])*T*S;
				sig_a[j]=geo_sig[j]-(geo_sig[j]-sig_0[j])*T;
			}
		}
		else if(f_yield>0) {
			cout<<endl<<"some problem"<<endl;
			//exit_tn( -YES );
		}

		/*********************************************/

		//double dev = vari.d_ept[0]+vari.d_ept[4]+vari.d_ept[8];
		//double dv = dev*v;
		//v += dv;
		
		vari.hisv[1]=I_0/3;
		vari.hisv[2]=sig_a[0];
		vari.hisv[3]=sig_a[1];
		vari.hisv[4]=sig_a[2];
		vari.hisv[5]=sig_a[4];
		vari.hisv[6]=sig_a[5];
		vari.hisv[7]=sig_a[8];
		vari.hisv[8]=sig_b[0];
		vari.hisv[9]=sig_b[1];
		vari.hisv[10]=sig_b[2];
		vari.hisv[11]=sig_b[4];
		vari.hisv[12]=sig_b[5];
		vari.hisv[13]=sig_b[8];
	}
	double lodeb=0, lodea=0, I=0, J=0, lode=0; 	
	double rolode=0, rolodea=0;
	double tmparray[3][3]={{0,0,0},{0,0,0},{0,0,0}};
	calc_IJlode(geo_sig, I, J, lode, false, tmparray, tmparray, tmparray);
	double Jnorm_a=0, Jnorm_b=0, Inorm_a=0, Inorm_b;
	calc_IJlode(sig_sigb, Inorm_b, Jnorm_b, lodeb, false, tmparray, tmparray, tmparray);
	calc_IJlode(sig_siga, Inorm_a, Jnorm_a, lodea, false, tmparray, tmparray, tmparray);
	calc_rolode(lode, r_m, rolode, vari.ddum[0], false); 
	calc_rolode(lodeb, r_m, rolodeb, drolodenbdlodenb, true); 
	calc_rolode(lodea, r_m, rolodea, vari.ddum[0], false); 
	calc_rolode(lodeb, r_mfl, rolodefl, drolodefldlodenb, true); 

	m=m_cmp*rolodeb;
	m_flow=m_flow_K0*rolodefl;
	double m_bound=m_cmp*rolode;
	double m_hist=m_cmp*rolodea;
	double m_yield=m;
	if((m_bound<=0)||(m_hist<=0)||(m_yield<=0)||(m_flow<=0)) {
		cout<<"error in AI3-SKH, location 7, element "<<vari.element<<endl;
		//exit_tn( -YES );
	}

	f_yield=(3*Jnorm_b*Jnorm_b/(m_yield*m_yield)+scalar_power(Inorm_b/3, 2)-
		T*T*S*S*I_0*I_0/9)/2;
	f_hist=(3*Jnorm_a*Jnorm_a/(m_hist*m_hist)+scalar_power(Inorm_a/3, 2)-
		T*T*I_0*I_0/9)/2;
	f_bound=(3*J*J/(m_bound*m_bound)+scalar_power((I-I_0)/3, 2)-
		I_0*I_0/9)/2;
	g=(3*Jnorm_b*Jnorm_b/(m_flow*m_flow)+scalar_power(Inorm_b/3, 2)-
		T*T*S*S*I_0*I_0/9)/2;
	f=f_yield;

	vari.other_f[0] = f_hist;
	vari.other_f[1] = f_bound;

	double powv=N-lambda*log(2*I_0/3)+kappa*log(2*I_0/I);
	v=scalar_power(EULER, powv);
	vari.hisv[0]=v-1;
}

void Tskh_general_plasti::give_matrix(bool only_f, bool elpl_matrix) {
	
	double dg[MDIM*MDIM];
	array_set(dg, 0, MDIM*MDIM);
	array_set(df, 0, MDIM*MDIM);
	double lodeb=0; 	
	double sig_sigb[9];
	array_set(sig_sigb, 0, 9);
	array_subtract(vari.geo_sig, sig_b, sig_sigb, 9);
	double Jnorm_b=0, Inorm_b;
	double dJnbdsig[3][3];
	double dInbdsig[3][3];
	double dlodenbdsig[3][3];
	array_set(*dJnbdsig, 0, 9);
	array_set(*dInbdsig, 0, 9);
	array_set(*dlodenbdsig, 0, 9);
	calc_IJlode(sig_sigb, Inorm_b, Jnorm_b, lodeb, true, dInbdsig, dJnbdsig, dlodenbdsig);
	double dfdJnb=0, dfdInb=0, dfdlodenb=0, dgdJnb=0, dgdInb=0, dgdlodenb=0;
	if((m<=0)||(m_flow<=0)||(rolodeb<=0)||(rolodefl<=0)) {
		cout<<"error in AI3-SKH, location 8, element "<<vari.element<<endl;
		//exit_tn( -YES );
	}
	dfdInb=dgdInb=Inorm_b/9;
	dfdJnb=3*Jnorm_b/(m*m);
	dgdJnb=3*Jnorm_b/(m_flow*m_flow);	
	
	double dfdrolodenb=-3*Jnorm_b*Jnorm_b/(m_cmp*m_cmp*scalar_power(rolodeb,3));
	double dgdrolodenb=-3*Jnorm_b*Jnorm_b/(m_flow_K0*m_flow_K0*scalar_power(rolodefl,3));
	dfdlodenb = dfdrolodenb * drolodenbdlodenb;
	dgdlodenb = dgdrolodenb * drolodefldlodenb;

	for (int i=0; i<3; i++ ) {
		for (int j=0; j<3; j++ ) {
	          	df[3*i+j]= dfdInb*dInbdsig[i][j]+dfdJnb*dJnbdsig[i][j]+
				dfdlodenb*dlodenbdsig[i][j];
			df[3*i+j]*=-1;	//Triax requires solid mechanics convention
	          	dg[3*i+j]= dgdInb*dInbdsig[i][j]+dgdJnb*dJnbdsig[i][j]+
				dgdlodenb*dlodenbdsig[i][j];
			dg[3*i+j]*=-1;	//Triax requires solid mechanics convention
		}
	}

	double dfdepp[9];
	array_set(dfdepp, 0, 9);
	double dfdI0=-T*T*S*S*I_0/9;
	double dI0devp=I_0/(lambda-kappa);
	double dsigbdepv[3][3];
	array_set(*dsigbdepv, 0, 9);	
	for (int i=0; i<3; i++ ) {
		for (int j=0; j<3; j++ ) {
			dsigbdepv[i][j]=sig_b[3*i+j]/(lambda-kappa);
		}
	}
	double tmpA[3][3];
	array_set(*tmpA, 0, 9);
	for (int i=0; i<3; i++ ) {
		for (int j=0; j<3; j++ ) {
			tmpA[i][j]=-dfdInb*dInbdsig[i][j]-dfdJnb*dJnbdsig[i][j]-dfdlodenb*dlodenbdsig[i][j];
		}
	}	
	double tmpB=array_inproduct(*tmpA, *dsigbdepv, 9)+dfdI0*dI0devp;
	for (int i=0; i<3; i++ ) {
		for (int j=0; j<3; j++ ) {
			dfdepp[3*i+j]=tmpB*vari.kron_delta[i][j];
			dfdepp[3*i+j] *= -1; 	//Triax requires solid mechanics convention
		}
	}	

	double b1max=2*I_0*(1-T)/3;
	if(b1max==0) b1max=TINY;
	double beta[9];
	array_set(beta,0,9);
	for (int j=0; j<9; j++ ) {
		beta[j]=(vari.geo_sig[j]-sig_b[j])/(T*S)+sig_0[j]-
			(vari.geo_sig[j]-sig_b[j])/S - sig_a[j];
	}
	double b1=array_inproduct(beta, df, 9)*3/(S*T*I_0);		
	if(vari.other_f[1]>0) b1=0;
	double H1=S*S*scalar_power(scalar_dabs(b1)/b1max, psi)*I_0*I_0*I_0/((lambda-kappa)*27);	
	double b2max=2*T*I_0*(1-S)/3;
	if(b2max==0) b2max=TINY;
	double gama[9];
	array_set(gama,0,9);
	for (int j=0; j<9; j++ ) {
		gama[j]=(vari.geo_sig[j]-sig_b[j])/S+sig_a[j]-vari.geo_sig[j];
	}
	double b2=array_inproduct(gama, df, 9)*3/(T*S*I_0);		
	if(vari.other_f[0]>0) b2=0;
	double H2=scalar_power(T*scalar_dabs(b2)/b2max, psi)*I_0*I_0*I_0/((lambda-kappa)*27);	
	
	double H0=0;
	for(int i=0; i<MDIM*MDIM; i++) H0-=dg[i]*dfdepp[i];
	double H = H0+H1+H2;

	array_move(df, vari.df_trial, MDIM*MDIM);
	if (!only_f) make_matrix_Cep(df, dg, H, elpl_matrix);
}

void Tskh_general_plasti :: calc_rolode(double lode, double beta, double &rolode, 
	double &drolodedlode, bool alsoderivative) {
	//eliptical lode dependence after william and warnkle
	lode*=-1;	//w&w dependence defined for stress negative in compression	
	double tmpA = 2*(1-beta*beta)*cos(PI/6-lode);
	double insqrt=(4*(1-beta*beta)*scalar_power(cos(PI/6-lode),2)+beta*(5*beta-4));
	if(insqrt<0) {
		cout<<"error in AI3-SKH, location 9, element "<<vari.element<<endl;
		//exit_tn( -YES );
	}
	double tmpB = (2*beta-1)*sqrt(insqrt);
	double tmpC = 4*(1-beta*beta)*scalar_power(cos(PI/6-lode),2)+(2*beta-1)*(2*beta-1);
	if(tmpC==0) {
		cout<<"error in AI3-SKH, location 10, element "<<vari.element<<endl;
		//exit_tn( -YES );
	}
	
	rolode=(tmpA+tmpB)/tmpC;
	drolodedlode=0;

	if(alsoderivative) {	
		double insqrt=4*(1-beta*beta)*scalar_power(cos(PI/6-lode),2)+beta*(5*beta-4);
		if(insqrt<0) {
			cout<<"error in AI3-SKH, location 11, element "<<vari.element<<endl;
			//exit_tn( -YES );
		}		
		double u = (2*beta-1)*sqrt(insqrt)+ 2*(1-beta*beta)*cos(PI/6-lode);
		double v = 4*(1-beta*beta)*scalar_power(cos(PI/6-lode),2)+(2*beta-1)*(2*beta-1);
		double udash_den=4*(1-beta*beta)*scalar_power(cos(PI/6-lode),2)+beta*(5*beta-4);
		if((udash_den<=0)||(v==0)) {
			cout<<"error in AI3-SKH, location 12, element "<<vari.element<<endl;
			//exit_tn( -YES );
		}		
		double udash = 2*(1-beta*beta)*sin(PI/6-lode)+(2*(2*beta-1)*(1-beta*beta)*sin(PI/3-2*lode))/
			sqrt(udash_den);
		double vdash = 4*(1-beta*beta)*sin(PI/3-2*lode);
		drolodedlode = (udash*v-vdash*u)/(v*v);
	}
	drolodedlode*=-1;	//bloody solid mechanics convention :-)
}
