
#include <boost/algorithm/string/trim.hpp>
/////////////////////////////////////////////////////////////////////////////////////
// Does what it says on the tin.
void save_output( phi_state &phiS, chi_state &chiS ){
	phiS.calc_dt_C();
	phiS.checkGauge();
	chiS.calcGradients();
	//make sure we are up to date with things
	chiS.fit_GCOC();
	phiS.calc_exp_Uen();
	phiS.calc_norms();
	//~ chiS.calc_nuc_phase();
	// print info to screen
	cout<<step<<"\t"<<t<<"\t"<<1.-chiS.nd_norm()<<endl;
	
	
	
	// begin writing things of interest to file 
	// NOTE: the magnitude squared of the elec coeffs
	// are printed to file in phiS.calc_norms()
	//~ for(int c = 0; c<n_states; c++){
		//~ write_vector1D(phiS.coeff_set[c], toString(file_n) + "C"+toString(c+1));
		
		//~ write_vector1D(phiS.T1_vec[c], toString(file_n) +  "T1_"+toString(c+1));
		//~ write_vector1D(phiS.T2_vec[c], toString(file_n) +  "T2_"+toString(c+1));
		//~ write_vector1D(phiS.T3_vec[c], toString(file_n) +  "T3_"+toString(c+1));
	//~ }		
	
	//~ write_vector1D(chiS.ft_grad, toString(file_n) +  "chi_ft_grad");
	//~ write_vector1D(chiS.ft_grad_sq, toString(file_n) +  "chi_ft_grad_sq");
	
	
	//~ write_vector1D(phiS.gauge_cond, toString(file_n) +  "gauge_cond");
	
	write_vector1D(phiS.Avec, toString(file_n) +  "Avec");
	write_vector1D(phiS.divAvec, toString(file_n) +  "divAvec");
	
	
	write_vector1D(phiS.TDPES, toString(file_n) +  "TDPES");
	//~ write_vector1D(phiS.dt_C_set[0], toString(file_n) +  "dt_C_1");
	//~ write_vector1D(phiS.dt_C_set[1], toString(file_n) +  "dt_C_2");
	
	//~ write_vector1D(phiS.UenC[0], toString(file_n) +  "Uen1");
	//~ write_vector1D(phiS.UenC[1], toString(file_n) +  "Uen2");
	
	vector<double> nd = chiS.density();
	write_vector1D( nd , toString(file_n) +  "nd");
	write_vector1D( chiS.GCOC , toString(file_n) +  "GCOC");
	//~ write_vector1D( chiS.H_chi , toString(file_n) +  "H_chi");
	write_vector1D(global_mask, toString(file_n) +  "global_mask");
	
	write_vector1D(phiS.exp_Uen, toString(file_n) +  "expUen");
	//~ write_vector1D(phiS.exp_Uen_trunc, toString(file_n) +  "expUen_trunc");
	write_vector1D(chiS.vec, toString(file_n) +  "chi");
	//~ write_vector1D(chiS.grad_R, toString(file_n) +  "grad_chi");
	//~ write_vector1D(chiS.grad_sq_R, toString(file_n) +  "grad_sq_chi");
	
	//~ write_vector1D(chiS.b, toString(file_n) +  "chi_b");
	
	//~ write_vector1D(phiS.stability_track, toString(file_n) +  "stability_track");
	
	// construct and write nuclear norm time series (norm dev as fn of time)
	nd_arr[0] = t;
	nd_arr[1] = abs(1.-chiS.nd_norm());
	nd_norm_timeseries.push_back(nd_arr);
	file_n++;	
}

/////////////////////////////////////////////////////////////////////////////////////
// Parses config file (in ad-hoc,fixed, timey-wimey kind of way)
void parse_config( config &c, string prev_str, string str ){
	stringstream iss( str );
	if(prev_str == "sim_name"){
		c.sim_name = str;
		
	} else if(prev_str == "basis"){
		int number;
		while ( iss >> number ){
			c.basis.push_back( number );
		}
	} else if(prev_str == "dt"){
		c.dt = stod(str);
		
	} else if(prev_str == "tmax"){
		c.tmax = stod(str);
		
	} else if(prev_str == "model"){
		c.model = str;
		
	} else if(prev_str == "sub_model"){
		if(str == "none"){
			str.clear();
		}
		c.sub_model = str;
		
	}  else if(prev_str == "step_gap"){
		c.step_gap = stoi(str);
		
	} else if(prev_str == "n_threads"){
		c.n_threads = stoi(str);
		
	} else if(prev_str == "have_phi"){
		c.have_phi = stoi(str);
		
	} else if(prev_str == "edge_tols"){
		double number;
		while ( iss >> number ){
			c.edge_tols.push_back( number );
		}
		
	} else if(prev_str == "init_state"){
		c.init_state = stoi( str );
		
	} else if(prev_str == "R0"){
		c.R0 = stod( str );
		
	} else if(prev_str == "k0"){
		c.k0 = stod( str );
		
	} else if(prev_str == "sigsq"){
		c.sigsq = stod( str );
		
	} else if(prev_str == "mask_width"){
		c.mask_width = stod( str );
	
	} else if(prev_str == "nuc_mass"){
		c.nuc_mass = stod( str );
	
	} else if(prev_str == "Uen_weights"){
		double number = 0.;
		while ( iss >> number ){
			c.weights.push_back( number );
		}
	
	} else if(prev_str == "A_re"){
		c.A_re = stoi(str);
	
	} else if(prev_str == "sub_imUen"){
		c.sub_imUen = stoi(str);
	
	} else if(prev_str == "mask_coeffs"){
		c.mask_coeffs = stoi(str);
	
	} else if(prev_str == "mask_uen"){
		c.mask_uen = stoi(str);
	
	} else if(prev_str == "mask_gcoc"){
		c.mask_GCOC = stoi(str);
	
	
	} else if(prev_str == "fixed_chi"){
		c.fixed_chi = stoi(str);
		
	} else if(prev_str == "mask_dt"){
		c.mask_dt = stoi(str);
	
	} else if(prev_str == "subGauge"){
		c.subGauge = stoi(str);
	
	} else if(prev_str == "free_chi"){
		c.free_chi = stoi(str);
		
	} else if(prev_str == "bdry_mask_edge"){
		c.bdry_mask_edge = stod(str);
	} 
	else if(prev_str == "dir_uen_tdpes"){
		c.dir_uen_tdpes = stod(str);
	} 
	else if(prev_str == "ctmqc_uen"){
		c.ctmqc_uen = stoi(str);
	}	
	else if(prev_str == "mask_A"){
		c.mask_A = stoi(str);
	}
	else if(prev_str == "re_tdpes"){
		c.re_tdpes = stoi(str);
	} 
	else if(prev_str == "phi_RK45"){
		c.phi_RK45 = stoi(str);
	}
	else if(prev_str == "phi_BE"){
		c.phi_BE = stoi(str);
	}
	else if(prev_str == "phi_CN"){
		c.phi_CN = stoi(str);
	}
	else if(prev_str == "ft_deriv"){
		c.ft_deriv = stoi(str);
	}
	
	
}
/////////////////////////////////////////////////////////////////////////////////////
// Reads and parses config file and also prints config file to screen
// Can also save a copy of the config to output folder
void get_config( config &in, string config_file ){
	print("#################################################");
	print("        Current configuration (config.dat)");
	print("#################################################");
	int line_num = 0;
	vector<string> config_text;
	ifstream file;
	string line, prev_line, head;
	file.open( config_file );
	
	if( file ){
		while(getline(file,line)){
			config_text.push_back(line);
			head = line.at(0);
			if(head == "#"){
				prev_line = line;
				prev_line.erase(0, 1); // gets rid of # bit - probably not needed given manual parsing above
				boost::algorithm::trim(prev_line); //trims spaces
				// do nothing, for now - might want to print to log later
			} else {
				parse_config( in, prev_line, line);
			}
		}
	} else {
		cout<<"Config file '" + in.model +"/"+ config_file + "' not found, exiting."<<endl;
		exit(0);
	}
	
	
	
	// Check if output folder exists with this name already
	// Some other crap to check if we want to overwrite what's in the folders or just remove 
	// Also does some stuff with the imgs folder creation (still need better/global script for this)
	struct stat info;
	outputFolder = in.model + "/output_" + in.sim_name +"/";
	string outputFolder_imgs = in.model + "/imgs_" + in.sim_name +"/";
	string out_command = "mkdir "+outputFolder; //+ " " + outputFolder_imgs ;
	string rm_command = "rm "+outputFolder+"*" ;//+ " " + outputFolder_imgs+"*" ;
	int systemRes;
	if( stat( outputFolder.c_str(), &info ) != 0 ){
		systemRes = system( out_command.c_str() ); // stores error val of this function - not used, just stops warning when compiled :)
	} else {
		cout<<"Output folder with that sim_name already exists! Overwrite?  (y/n)"<<endl;
		string ow;
		cin>>ow;
		if(ow == "y"){
			cout<<"Empty existing directories (output and imgs)? (y/n)"<<endl;
			cin>>ow;
			if(ow == "y"){
				systemRes = system( rm_command.c_str() );
			}
		} else {
			cout<<"Well lah-de-dah...exiting."<<endl;
			exit(0);
		}
	}

	//////////////////////////////////////////////////////////////
	// Apply configuration settings
	BO_basis = in.basis;
	n_states = BO_basis.size();
	have_phi = in.have_phi;
	model_loc = in.model+"/";
	submodel = in.sub_model;
	
	for(int c = 0; c<config_text.size(); c++){
		cout<<config_text[c]<<endl;
	}
	
	ofstream config_write;
	string config_out_file = outputFolder + in.sim_name + ".used_config";
	config_write.open(config_out_file.c_str(), ios::out);
	for(int i = 0; i<config_text.size(); i++){
		config_write<<config_text[i]<<endl;
	}
	config_write.close();
	
}

/////////////////////////////////////////////////////////////////////////////////////
//END OF FILE
/////////////////////////////////////////////////////////////////////////////////////
