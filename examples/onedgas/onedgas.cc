///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../../libmd.cc"
#include "../../tools/BaX/BaX.h"
#include "../../tools/springio/springio.cc"
#include <boost/program_options.hpp>
#include <random>
#include <cmath>
#include <string>
#include <cstdio>
#include "pcg/pcg_random.hpp"
#include <sys/stat.h>

namespace po = boost::program_options;

pcg64 RNG;
std::normal_distribution<double> NormalDistribution {0, 1};

ldf globalTime=0;


template <ui dim> void brownianforce(ui i, vector<ui>* particles, vector<ldf>* parameters, void* sys) {
    (void) particles;
    ldf sqrtDt = parameters->at(0);
    for (ui d=0; d<dim; d++) {
	ldf rf = NormalDistribution(RNG);
	SYS->particles[i].F[d] += sqrtDt*rf;
    }
}

template <ui dim> void weimueller(ui i, vector<ui>* particles, vector<ldf>* parameters, void* sys) {
    ldf v1 = parameters->at(0); ldf v2 = parameters->at(1);
    ldf p = parameters->at(2); ldf q = parameters->at(3); ldf T = parameters->at(4);
    ldf xi = SYS->particles[i].x[0];
    SYS->particles[i].F[0] += -v1*q*sin(q*xi) - v2*p*sin(p*xi-2*M_PI*globalTime/T);
}

template <ui dim> void weimuellerdir(ui i, vector<ui>* particles, vector<ldf>* parameters, void* sys) {
    ldf v1 = parameters->at(0); ldf v2 = parameters->at(1);
    ldf p = parameters->at(2); ldf q = parameters->at(3); ldf T = parameters->at(4);
    ldf xi = SYS->particles[i].x[0];
    ldf yi = SYS->particles[i].x[1];
    SYS->particles[i].F[0] += -v1*q*sin(q*xi) - v2*p*sin(p*xi-2*M_PI*yi/T);
    SYS->particles[i].F[1] += v2*2*M_PI/T*sin(p*xi-2*M_PI*yi/T);
}

template <ui dim> void linetension(ui i, vector<ui>* particles, vector<ldf>* parameters, void* sys) {
    ldf tau_x = parameters->at(0);
    ldf tau_y = parameters->at(1);
    SYS->particles[i].F[0] += tau_x;
    SYS->particles[i].F[1] += tau_y;
}


ldf urand()
{

    static unsigned long int rseed=666UL;
    return (rseed=16807UL*rseed%2147483647UL)/2147483648.0;
}

int main(int argc, char** argv)
{
    __libmd__info();


    // parse options
    
    po::options_description desc("Allowed options");
    po::variables_map vm; 
    try {
        
        desc.add_options()
            ("help", "produce help message")
	    ("png", "save pngs")
            ("kT", po::value<ldf>()->default_value(0.), "kbT")
            ("gamma", po::value<ldf>()->default_value(0.), "drag")
            ("Lx", po::value<ldf>()->default_value(1.), "box size x")
            ("Ly", po::value<ldf>()->default_value(10.), "box size y")
            ("tau", po::value<ldf>()->default_value(0.), "tension")

            ("v1", po::value<ldf>()->default_value(0.), "potential v1")
            ("v2", po::value<ldf>()->default_value(0.), "potential v2")
            ("repx", po::value<ui>()->default_value(10), "reps of potential along x")
            ("repy", po::value<ldf>()->default_value(5), "reps of potential along y")
            ("wmp", po::value<ldf>()->default_value(1.), "Wei-Mueller p")
            ("wmq", po::value<ldf>()->default_value(2.), "Wei-Mueller q")


	    
            ("N", po::value<ui>()->default_value(1), "# polymer chains")
            ("np", po::value<ui>()->default_value(1), "# monomers in chain")	    

	    ("h", po::value<ldf>()->default_value(0.001), "timestep")
	    ("frames", po::value<ui>()->default_value(1000), "# frames")
	    ("step", po::value<ui>()->default_value(1000), "# steps per frame")

	    ("yb", po::value<ldf>()->default_value(0.), "yukawa b")
	    ("yk", po::value<ldf>()->default_value(200.), "yukawa kappa")
	    ("linksize", po::value<ldf>()->default_value(-1.), "link size")
	    ("linkk", po::value<ldf>()->default_value(1000000.), "link spring constant")

        ;

        po::store(po::parse_command_line(argc, argv, desc), vm);
        //~ po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);
        po::notify(vm);    

        if (vm.count("help")) {
            cout << desc << "\n";
            exit(1);
        }
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        exit(1);
    }
    catch(...) {
        cerr << "options: exception of unknown type!\n";
    }

    
    bool png = false;
    if (vm.count("png")) png = true;

    ui N = vm["N"].as<ui>();

    ldf *x = new ldf[N]();

    ldf *zero = new ldf[N]();
    
    unsigned int W=500,H=100;

    bitmap bmp(W,H);
//    color pix[]={RED,GREEN};

    color pix[]={NAVY,DARKBLUE,MEDIUMBLUE,BLUE,DARKGREEN,GREEN,TEAL,DARKCYAN,DEEPSKYBLUE,DARKTURQUOISE,MEDIUMSPRINGGREEN,LIME,SPRINGGREEN,AQUA,CYAN,MIDNIGHTBLUE,DODGERBLUE,LIGHTSEAGREEN,FORESTGREEN,SEAGREEN,DARKSLATEGRAY,LIMEGREEN,MEDIUMSEAGREEN,TURQUOISE,ROYALBLUE,STEELBLUE,DARKSLATEBLUE,MEDIUMTURQUOISE,INDIGO,DARKOLIVEGREEN,CADETBLUE,CORNFLOWERBLUE,MEDIUMAQUAMARINE,DIMGRAY,SLATEBLUE,OLIVEDRAB,SLATEGRAY,LIGHTSLATEGRAY,MEDIUMSLATEBLUE,LAWNGREEN,CHARTREUSE,AQUAMARINE,MAROON,PURPLE,OLIVE,GRAY,SKYBLUE,LIGHTSKYBLUE,BLUEVIOLET,DARKRED,DARKMAGENTA,SADDLEBROWN,DARKSEAGREEN,LIGHTGREEN,MEDIUMPURPLE,DARKVIOLET,PALEGREEN,DARKORCHID,YELLOWGREEN,SIENNA,BROWN,DARKGRAY,LIGHTBLUE,GREENYELLOW,PALETURQUOISE,LIGHTSTEELBLUE,POWDERBLUE,FIREBRICK,DARKGOLDENROD,MEDIUMORCHID,ROSYBROWN,DARKKHAKI,SILVER,MEDIUMVIOLETRED,INDIANRED,PERU,CHOCOLATE,TAN,LIGHTGREY,PALEVIOLETRED,THISTLE,ORCHID,GOLDENROD,CRIMSON,GAINSBORO,PLUM,BURLYWOOD,LIGHTCYAN,LAVENDER,DARKSALMON,VIOLET,PALEGOLDENROD,LIGHTCORAL,KHAKI,ALICEBLUE,HONEYDEW,AZURE,SANDYBROWN,WHEAT,BEIGE,WHITESMOKE,MINTCREAM,GHOSTWHITE,SALMON,ANTIQUEWHITE,LINEN,LIGHTGOLDENRODYELLOW,OLDLACE,RED,FUCHSIA,MAGENTA,DEEPPINK,ORANGERED,TOMATO,HOTPINK,CORAL,DARKORANGE,LIGHTSALMON,ORANGE,LIGHTPINK,PINK,GOLD,PEACHPUFF,NAVAJOWHITE,MOCCASIN,BISQUE,MISTYROSE,BLANCHEDALMOND,PAPAYAWHIP,LAVENDERBLUSH,SEASHELL,CORNSILK,LEMONCHIFFON,FLORALWHITE,SNOW,YELLOW,LIGHTYELLOW,IVORY,WHITE};

    if (png) bmp.fillup(BLACK);
    md<1> sys(N);
    sys.simbox.L[0]=vm["Lx"].as<ldf>();
    sys.simbox.bcond[0]=BCOND::PERIODIC;

    vector<ldf> selfavoid={vm["yb"].as<ldf>(),vm["yk"].as<ldf>(),12};

    if (selfavoid[0] > 0) sys.add_typeinteraction(0,0,POT::POWERLAW,selfavoid);

    ldf posx = -sys.simbox.L[0]/2.+sys.simbox.L[0]/N/2.;
    for (ui j=0; j < N; j++) {
	x[j] = posx;
	posx += sys.simbox.L[0]/N;
    }

    sys.import_pos(x);
//    sys.indexdata.method=INDEX::BRUTE_FORCE;
    sys.set_rco(sys.simbox.L[0]*5./N);
    sys.set_ssz(sys.simbox.L[0]*10./N);
    

    // damping
    ldf gam = vm["gamma"].as<ldf>();
    ldf kT = vm["kT"].as<ldf>();
    sys.set_damping(gam);

    sys.integrator.h = vm["h"].as<ldf>();
    
    vector<vector<ui> > dummyptr;
    // brownian external force?
    vector<ldf> Dt = {sqrt(2*kT*gam/sys.integrator.h)};
    
    ui bf = sys.f.add((extforceptr<1>) brownianforce<1>);
    ui bftype = sys.add_forcetype(bf,dummyptr,Dt);
    sys.assign_all_forcetype(bftype);
    
    // wei-mueller potential
    ldf onestepfactor = 2*M_PI/(sys.simbox.L[0]/vm["repx"].as<ui>());

    vector<ldf> WMparam = {vm["v1"].as<ldf>(),vm["v2"].as<ldf>(),
			   vm["wmp"].as<ldf>()*onestepfactor,vm["wmq"].as<ldf>()*onestepfactor,
			   vm["repy"].as<ldf>()};
    ui wmf = sys.f.add((extforceptr<1>) weimueller<1>);
    ui wmftype = sys.add_forcetype(wmf,dummyptr,WMparam);
    sys.assign_all_forcetype(wmftype);

    // data folder
    char fname[100];
    sprintf(fname, "./onedgas_Lx%1.2f_N%d_kT%g_v1%g_v2%g_yb%g_p%1f_q%1f/",sys.simbox.L[0],
	    N,kT,WMparam[0],WMparam[1],selfavoid[0],vm["wmp"].as<ldf>(),vm["wmq"].as<ldf>());
    string datafolder = string(fname);
    string pngfolder = datafolder + "img/";
    string outfolder = datafolder + "pts/";
	
    mkdir(datafolder.c_str(), S_IRWXU);
    mkdir(pngfolder.c_str(), S_IRWXU);
    mkdir(outfolder.c_str(), S_IRWXU);

    // output all parameter values
    ofstream cf;
    cf.open(datafolder+"config_file");
    for (const auto& it:vm) {
	cf << it.first.c_str() << ":\t";
	auto& value = it.second.value();
	if (auto v = boost::any_cast<ui>(&value))
	    cf << *v << endl;
	else if (auto v = boost::any_cast<ldf>(&value))
	    cf << *v << endl;
    }
    cf.close();

    
    sys.import_vel(zero);
    sys.export_pos(x);

    if (png) {
	for(ui i=0;i<N;i++) {
	    bmp.solidkykel(2.0,W*x[i]/sys.simbox.L[0]+W/2.0,H/2., pix[i]);
	}
    
	bmp.save_png_seq(const_cast<char*>((pngfolder+"sim").c_str()));
    }
    char frameid[100];
    cout << endl;
    
    for(ui k=0;k<vm["frames"].as<ui>();k++)
    {
	for (ui step=0; step < vm["step"].as<ui>();step++) {
	    sys.timestep();
	    globalTime += sys.integrator.h;
	}

	if (png) {
	    sys.export_pos(x);
	    bmp.fillup(BLACK);

	    for(ui i=0;i<N;i++) {
		bmp.solidkykel(2.0,W*x[i]/sys.simbox.L[0]+W/2.0,
				   H/2., pix[i]);
	    }
	    bmp.save_png_seq(const_cast<char*>((pngfolder+"sim").c_str()));
	}
	sprintf(frameid,"%08d",k);

	write_points_x(outfolder+"sim"+string(frameid)+".pts", sys);
	cout << "\r" << k+1 << "/" << vm["frames"].as<ui>();
	cout.flush();
    }
    cout << endl;
    return EXIT_SUCCESS;
}
