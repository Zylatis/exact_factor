/////////////////////////////////////////////////////////////////////////////////////
//STORES MISC FUNCTIONS USED IN EXACT FACTOR CRANK NICHOLSON CODE
/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
// Forward declarations stuff
// NACVs used in this library but also those functions rely on definitions in *this* library
// so fixes circular dependencies
void NACVs( int i, int j, int size, vector<vector< dcomp > > &flat_BO_vecs );
void calc_NACVs();

/////////////////////////////////////////////////////////////////////////////////////
// Fn to convert to string
string toString(int param){
 std::stringstream ss(stringstream::in|stringstream::out);
 ss<<param;
 string t=ss.str();
 return t;
}

/////////////////////////////////////////////////////////////////////////////////////
// polynomimal fitting function and deriv
double poly_fn( double x, vector<double> &pars ){
		return 1.*pars[0] + 1.*pars[1]*x  + 1.*pars[2]*pow(x,2.) + 1.*pars[3]*pow(x,3.);
	}


/////////////////////////////////////////////////////////////////////////////////////
// forward declarations of file writing functions
void write_vector1D(  vector< double>  array, string file_name, string outF = outputFolder);
void write_vector1D(  vector< dcomp>  array, string file_name, string outF = outputFolder);
void write_vector1D(  vector< int>  array, string file_name, string outF = outputFolder);
//~ template <typename T> 
//~ void write_vector2D(  vector< vector< T > >  array, string file_name, string outF = outputFolder);

/////////////////////////////////////////////////////////////////////////////////////
// Print 1d vector
template <typename T> 
void print_col_vector(T  &col_vec){
	for(int i = 0; i<col_vec.size();i++){
		cout<<col_vec[i]<<endl;
	}
}


/////////////////////////////////////////////////////////////////////////////////////
vector<dcomp> mean_vec( vector<dcomp> &v_old, vector<dcomp> &v_new, double a){
	vector<dcomp> out;
	int s = v_old.size();
	out.resize( s );
	for(int R = 0; R<s; R++){
		out[R] = (1.-a)*v_old[R] + a*v_new[R];
	}
	return out;
}

/////////////////////////////////////////////////////////////////////////////////////
// Find position of max val in a vector
int get_max_pos( vector<double> &vec ){
	double cur_max = vec[0];
	int pos = 0;
	for(int c = 1; c< vec.size(); c++){
		if( vec[c] > cur_max ){
			pos = c;
			cur_max = vec[c];
		}
	}
	return pos;
}	
/////////////////////////////////////////////////////////////////////////////////////
// Find edges of density with fed in centre point to start - absolute tolerance

vector<int> find_edges( vector<double> &density, int centre, double tol){
	vector<int> out;
	out.resize(2);
	int s = density.size();
	for(int c = centre; c >= 0; c--){
		if( density[c] < tol){
			out[0] = c;
			break;
		}
	}
	
	for(int c = centre; c < s; c++){
		if( density[c] < tol){
			out[1] = c;
			break;
		}
	}
	return out;
}

/////////////////////////////////////////////////////////////////////////////////////
vector<int> find_mask_points( vector<double> &density, double tol ){
	vector<int> points;
	for( int R = 1; R < nR; R++ ){
		if( (density[R-1] < tol && density[R] > tol) || (density[R-1] > tol && density[R] < tol)){
			points.push_back(R);
		}
	}
	
	return points;
}

/////////////////////////////////////////////////////////////////////////////////////
// Convert 1D vector coord to 2D r-R coord
vector<int> get_coord( int c){
	vector<int> out;
	out.resize(2);
	int rpos = floor(c/nR);
	int Rpos = c-rpos*nR;
	out[0] = rpos;
	out[1] = Rpos;
	
	return out;
}

/////////////////////////////////////////////////////////////////////////////////////
//Read in BO wfn file
vector< vector<double> > FileRead( string filename){
    vector< vector<double> > vvv;
    int i = 0;
    string line;

    fstream in(filename.c_str());
    if(!in.is_open()){
		cout<<"Bugger: unable to open file "+filename+" - exiting."<<endl;
		exit(0);
	}
	while (std::getline(in, line))
    {
        double value;
        stringstream ss(line);

        vvv.push_back(vector<double>());

        while (ss >> value)
        {
            vvv[i].push_back(value);
        }
        ++i;
    }
    return vvv;
}

/////////////////////////////////////////////////////////////////////////////////////
// Use above code to read in 1D complex function where first col is reals, 2nd imag
vector<dcomp>  FileRead1D_dcomp( string filename){
	vector<vector< double> > raw = FileRead(filename);
	vector<dcomp> out;
	int s = raw.size();
	out.resize( s );
	for(int c = 0;c<s; c++){
		out[c] = raw[c][0]+ii*raw[c][1];
	}
	return out;
}

/////////////////////////////////////////////////////////////////////////////////////
// Calculate time in seconds
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

/////////////////////////////////////////////////////////////////////////////////////
// Drop duplicates from a vector
void get_unique( vector<double> &vec){
	set<double> s( vec.begin(), vec.end() );
	vec.assign( s.begin(), s.end() );
}

/////////////////////////////////////////////////////////////////////////////////////
// Extract column from multidimensional array
vector<double> get_col( vector<vector<double> > &array, int col ){
	vector<double>  temp;
	temp.resize( array.size() );
	for (int i=0; i < array.size(); i++){
		if(array[i].size()!=0){
			temp[i] = array[i][col];
		}
	}
	return temp;
}
/////////////////////////////////////////////////////////////////////////////////////
// calculate a full 1D norm
double calc_vec_norm( vector< dcomp > &vec, double spacing){
	double norm = 0.;
	for(int c = 0; c<vec.size(); c++){
		norm += abs(pow(vec[c],2.))*spacing;
	}
	return norm;
}

/////////////////////////////////////////////////////////////////////////////////////
// Lazy print fucntion
void print( string str ){
	cout<<str<<endl;
}

///////////////////////
vector<dcomp> calc_grad_vec( vector<dcomp> &vec, int dim, int l_boundary, int r_boundary ){
	
	vector<dcomp> grad;
	grad.resize(vec.size());
	for(int c = 0; c<vec.size(); c++){
		grad[c] = deriv(vec, c, 1, dim, l_boundary, r_boundary);
	}
	return grad;
}

///////////////////////
vector<dcomp> calc_grad_sq_vec( vector<dcomp> &vec, int dim, int l_boundary, int r_boundary ){
	
	vector<dcomp> grad;
	grad.resize(vec.size());
	for(int c = 0; c<vec.size(); c++){
		grad[c] = laplacian(vec, c, 1, dim, l_boundary, r_boundary);
	}
	
	return grad;
}

/////////////////////////////////////////////////////////////////////////////////////
// Run initial setup to allocate memory and other fun things
void do_setup(){
	print("#################################################");
	print("                 Running setup");
	print("#################################################");
	
	print("Import adiabatic PES: ");

    BOPES.resize( n_states ); 
	vector<vector<double> >  t1;
    for( int c = 0; c<n_states;c++){
		t1 = FileRead( model_loc + "BOPES/"+toString(BO_basis[c])+"_bopes"+submodel+".dat");
		vector<double> t2 = get_col(t1,1);
		BOPES[c].resize(t2.size());
		for(int x = 0; x<t2.size(); x++){	
			BOPES[c][x] = t2[x];
		}
	}
	
	R_points = get_col(t1,0);
	get_unique(R_points);
	nR = R_points.size();
	global_mask.resize(nR);
	
	R_min = R_points[0];
	R_max = R_points[nR-1];
	dR = R_points[1] - R_points[0];
	
	cout<<"nR = "<<nR<<endl;
	cout<<"dR = "<<dR<<endl;
	cout<<"n_states = "<<n_states<<endl;
	
	// initialize NACV tensor (or rather matrix of functions)
	NACV_1.resize( n_states);
	NACV_2.resize( n_states);
	grad_eps.resize( n_states );
	grad_sq_eps.resize( n_states );
	vector<dcomp> test_vec, out2;
	test_vec.resize(nR);
	out2.resize(nR);
	
	for(int i = 0; i<n_states; i++){
		NACV_1[i].resize(n_states);
		NACV_2[i].resize(n_states);
		
		grad_eps[i].resize( nR );
		grad_sq_eps[i].resize( nR );
		
		for(int j = 0; j<n_states; j++){
			NACV_1[i][j].resize(nR);
			NACV_2[i][j].resize(nR);
			fill(NACV_1[i][j].begin(), NACV_1[i][j].end(), 0.);
			fill(NACV_2[i][j].begin(), NACV_2[i][j].end(), 0.);
		}
		
		grad_eps[i] = calc_grad_vec(BOPES[i], 1, 0, nR-1);	
		grad_sq_eps[i] = calc_grad_sq_vec(BOPES[i], 1, 0, nR-1);			
	}
	
	
	if(have_phi){
		cout<<"Calc NACV's:"<<endl;
		calc_NACVs();
		
		empty_vec_nr.resize( nr*nR );
		unity_vec_nr.resize( nr*nR );
	
		fill(empty_vec_nr.begin(), empty_vec_nr.end(), 0.);
		fill(unity_vec_nr.begin(), unity_vec_nr.end(), 1.);

	} else if(!have_phi){
		cout<<"Importing NACVs:"<<endl;
		for(int i = 0; i<n_states; i++){
			for(int j = 0; j<n_states; j++){
				NACV_1[i][j] =  FileRead1D_dcomp( model_loc + "NACV/NACV1-"+toString( i + 1 )+toString( j + 1 )+submodel+".txt");
				NACV_2[i][j] =  FileRead1D_dcomp( model_loc + "NACV/NACV2-"+toString( i + 1 )+toString( j + 1 )+submodel+".txt");
				
				if(i == j ){
					fill(NACV_1[i][j].begin(),NACV_1[i][j].end(),0.);
				}
			}
		}
	}
	gradNACV_1.resize(n_states);
	for(int i = 0; i<n_states; i++){
		gradNACV_1[i].resize(n_states);
		for( int j = 0; j < n_states; j++ ){
			gradNACV_1[i].resize(nR);
			gradNACV_1[i][j] = calc_grad_vec( NACV_1[i][j], 1, 0, nR-1);
		}
	}

	nd_arr.resize( 2 );
	
	nd_edges_inner_smooth.resize(2);
	nd_edges_outer_smooth.resize(2);
	
	empty_vec_nR.resize( nR );	
	unity_vec_nR.resize( nR );
	fill(empty_vec_nR.begin(), empty_vec_nR.end(), 0.);
	fill(unity_vec_nR.begin(), unity_vec_nR.end(), 1.);
		
}

/////////////////////////////////////////////////////////////////////////////////////
// Functions to write 1D vectors as a 2D list with col 1 = Real, col 2 = Imag
void write_vector1D(  vector< dcomp>  array, string file_name, string outF){
	ofstream o;
	int size = array.size();
	string file = outF+file_name+".txt";
	o.open(file.c_str(), ios::out);
	for(int i = 0; i<size; i++){
		o<< std::fixed << std::setprecision( FILE_P ) << real(array[i])<<"\t"<<imag(array[i])<<endl;
	}
	o.close();
}

void write_vector1D(  vector< double>  array, string file_name, string outF){
	ofstream o;
	int size = array.size();
	string file = outF+file_name+".txt";
	o.open(file.c_str(), ios::out);
	for(int i = 0; i<size; i++){
		o<< std::fixed << std::setprecision( FILE_P ) << array[i]<<"\t"<<0.<<endl;
	}
	o.close();
}

void write_vector1D(  vector< int>  array, string file_name, string outF){
	ofstream o;
	int size = array.size();
	string file = outF+file_name+".txt";
	o.open(file.c_str(), ios::out);
	for(int i = 0; i<size; i++){
		o<< std::fixed << std::setprecision( FILE_P ) << array[i]<<"\t"<<0.<<endl;
	}
	o.close();
}

void write_vector2D( vector< vector< double >  >  array, string file_name, string outF = outputFolder){

	ofstream o;
	int size = array.size();
	string file = outF+file_name+".txt";
	o.open(file.c_str(), ios::out);
	for(int i = 0; i<size; i++){
		o<< std::fixed << std::setprecision( FILE_P ) << array[i][0]<<"\t"<<array[i][1]<<endl;
	}
	o.close();
}

void write_vector2D( vector< vector< int>  >  array, string file_name, string outF = outputFolder){

	ofstream o;
	int size = array.size();
	string file = outF + file_name+".txt";
	o.open(file.c_str(), ios::out);
	for(int i = 0; i<size; i++){
		o<< std::fixed << std::setprecision( FILE_P ) << array[i][0]<<"\t"<<array[i][1]<<endl;
	}
	o.close();
}


/////////////////////////////////////////////////////////////////////////////////////
//END OF FILE
/////////////////////////////////////////////////////////////////////////////////////
