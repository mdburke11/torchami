#include "ami_spec.hpp"

std::complex<double> AmiSpec::evaluate_sp_term(AmiBase::ami_parms &parms, AmiSpec::ami_sp_term &sp_term, NewAmiCalc::ext_vars &ev,   AmiBase::ami_vars &external, NewAmiCalc::k_vect_list_t &klist,   xi_t &xi_list){

std::complex<double> output(0,0);

AmiBase::ami_vars gprod_external=external;	

// TODO: is this sign swap necessary?	
for(int i=0; i< gprod_external.energy_.size(); i++){
 
gprod_external.energy_[i]=-xi_list[i];	
	
}		

std::complex<double> A_prod=eval_Aprod(sp_term.aprod_, xi_list, external.frequency_, klist, ev.MU_);

std::complex<double> gprod;
gprod=amibase.eval_gprod(parms, sp_term.ami_term_.g_list, gprod_external);
	
	
std::complex<double> fprod;
fprod=amibase.eval_fprod(parms, sp_term.ami_term_.p_list, gprod_external);	

std::complex<double> term_val(0,0);
std::complex<double> norm(0,0);

term_val=sp_term.ami_term_.sign*gprod*fprod;

std::complex<double> imag(0.,1.0);
norm=std::pow(-imag*M_PI/(2.0*xi_cutoff), sp_term.delta_count);

term_val=term_val*A_prod*norm;

output+= term_val;

return output;	
	
	
}





std::complex<double> AmiSpec::eval_Aprod(A_prod_t &Ap, xi_t &xi, AmiBase::frequency_t &freq, NewAmiCalc::k_vect_list_t &klist, std::complex<double> &mu){

	std::complex<double> output(1,0);

	// A(Sigma, X, E)  : Sigma: self-energy, X: is frequency, E: energy from k_vector
	for(int i=0; i< Ap.size(); i++){

		std::complex<double> this_X=get_X( Ap[i].x_, xi, Ap[i].x_alpha_, freq);

		NewAmiCalc::k_vector_t this_k=ami.construct_k(Ap[i].alpha_, klist);
		std::complex<double> this_E=eval_tb(1.,0., this_k, mu);

		std::complex<double> this_sigma=get_sigma(this_k, this_X);


		output=output*A_eval(this_sigma, this_X, this_E);



	}


	return output;

}

// todo: probably don't need to pass A if this function takes in X and E already
std::complex<double> AmiSpec::A_eval( std::complex<double> &sigma, std::complex<double> &X, std::complex<double> &E){

std::complex<double> output(0,0);

output=- sigma.imag()/M_PI/( std::pow(X-E-sigma.real(),2) +std::pow(sigma.imag(),2) );
return output;

}



std::complex<double> AmiSpec::eval_tb(double t, double tp, NewAmiCalc::k_vector_t &k, std::complex<double> &mu){

	std::complex<double> output(0,0);

	for(int i=0; i<k.size();i++){

	output+=-2.0*t*cos(k[i]);

	}

	// uncomment this block if want t' in dispersion
	// double term=-4.0*tp;
	// for(int i=0; i<k.size(); i++){

		// term=term*cos(k[i]);
	// }
	// output+=term;


	output -= mu;

	return output;

}


std::complex<double> AmiSpec::construct_energy(AmiBase::alpha_t &alpha, NewAmiCalc::k_vect_list_t &klist, std::complex<double> &mu){

std::complex<double> result=0;

NewAmiCalc::k_vector_t this_k=ami.construct_k(alpha, klist);
result=eval_tb(1.,0., this_k, mu);


return result;

}