// Strings library
#include <string>
using std::string;
using std::getline;

// Containers library
#include <vector>
using std::vector;

// Numerics library
#include <cmath>

// Input/output library
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <sstream>
using std::istringstream;

// Filesystem library
#include <filesystem>
using std::filesystem::exists;



void printlog(const string lfname = "pxrd.log");

void csvtotable(const string ifname, vector<vector<double>>& table);

struct peakdata {
    int number;     // peak number
    double theta2;  // 2theta (deg)
    double d2;      // d^-2 (angstrom^-2)
    int index;      // h^2 + k^2 + l^2
    double a;       // lattice constant
    double nr;      // Nelson-Riley function
};

void tabletopeaks(const vector<vector<double>>& table, vector<peakdata>& peaks, string& lattice_type, const double lambda = 1.5405);

void refine(const vector<peakdata>& peaks, double& lattice_const, double& sigma);

void addlog(const string ifname, const string ofname, const string lattice_type, const double lattice_const, const double sigma, const string lfname = "pxrd.log");



int main(int argc, char **argv){
    // Processing command line arguments
    if (argc == 1) {
        cerr << "Missing filename" << endl;
        return 1;
    }

    string arg1 = argv[1];
    if (arg1 == "-l") {
        printlog();
        return 0;
    }

    int argnum = 1;
    if (arg1 == "-o") {
        argnum = 2;
    }

    if ((argnum == 2) && (argc == 2)) {
        cerr << "Missing filename" << endl;
        return 1;
    }

    // Main process
    int errcode = 0;
    for (int i = argnum; i < argc; i++) {
        string ifname = argv[i];
        if (!exists(ifname)) {
            cerr << ifname << ": No such file" << endl;
            errcode = 2;
            continue;
        }

        string ofname = ifname;
        if (argnum == 1) {
            ofname.insert(ofname.rfind("."), "_out");
        }

        // Reading input file
        vector<vector<double>> table;
        csvtotable(ifname, table);

        // Determination of the index and lattice constant for each peak
        vector<peakdata> peaks;
        string lattice_type;
        tabletopeaks(table, peaks, lattice_type);

        // Writing output file
        cout << "Writing " << ofname << " ..." << endl;
        ofstream ofile(ofname);
        ofile << "number,2theta,d^-2,index,nr,a" << endl;
        for (int j = 0; j < peaks.size(); j++) {
            ofile << peaks.at(j).number << ',';
            ofile << peaks.at(j).theta2 << ',';
            ofile << peaks.at(j).d2 << ',';
            ofile << peaks.at(j).index << ',';
            ofile << peaks.at(j).nr << ',';
            ofile << peaks.at(j).a << endl;
        }
        ofile.close();

        // Refinement of the lattice constant
        double lattice_const, sigma;
        refine(peaks, lattice_const, sigma);
        addlog(ifname, ofname, lattice_type, lattice_const, sigma);
    }

    if (errcode == 0) {
        cout << "The command completed successfully!" << endl;
    }
    return errcode;
}



void printlog(const string lfname) {
    if (!exists(lfname)) {
        cout << "No log" << endl;
        return;
    }

    ifstream logfile(lfname);

    string line_s;
    getline(logfile, line_s);
    while (getline(logfile, line_s)) {
        istringstream line_iss(line_s);

        string field_s;
        getline(line_iss, field_s, ':');
        cout << field_s << "->";

        getline(line_iss, field_s, ':');
        cout << field_s << ":\n";

        getline(line_iss, field_s, ':');
        cout << "   " << field_s << ": ";

        getline(line_iss, field_s, ':');
        istringstream field_iss_1(field_s);
        double field_d;
        field_iss_1 >> field_d;
        cout << field_d << " ± ";

        getline(line_iss, field_s);
        istringstream field_iss_2(field_s);
        field_iss_2 >> field_d;
        cout << field_d << " Å\n" << endl;
    }
    logfile.close();
    return;
}



void csvtotable(const string ifname, vector<vector<double>>& table) {
    ifstream ifile(ifname);

    string line_s;
    while (getline(ifile, line_s)) {
        istringstream line_iss(line_s);
        vector<double> record;

        string field_s;
        while (getline(line_iss, field_s, ',')) {
            istringstream field_iss(field_s);
            double field_d = 0;
            field_iss >> field_d;
            record.push_back(field_d);
        }

        if (record.size()) {
            if (record.at(0)) {
                table.push_back(record);
            }
        }
    }

    ifile.close();
    return;
}



void tabletopeaks(const vector<vector<double>>& table, vector<peakdata>& peaks, string& lattice_type, const double lambda) {
    const double pi = acos(-1);

    for (int i = 0; i < table.size(); i++) {
        peakdata peak;
        peak.number = static_cast<int>(table.at(i).at(0));
        peak.theta2 = table.at(i).at(1);
        double angle = peak.theta2 * pi / 360;
        peak.d2 = pow((2*sin(angle)/lambda), 2);
        peak.nr = (1/sin(angle) + 1/angle) * pow(cos(angle), 2);
        peaks.push_back(peak);
    }

    int parameter;
    double ratio21 = peaks.at(1).d2 / peaks.at(0).d2;
    if (ratio21 > 2.33) {
        lattice_type = "diamond structure";
        parameter = 3;
    } else if (ratio21 < 1.67) {
        lattice_type = "FCC (face-centered cubic)";
        parameter = 3;
    } else if (peaks.size() >= 7) {
        double ratio71 = peaks.at(6).d2 / peaks.at(0).d2;
        if ( ratio71 < 7.5) {
            lattice_type = "BCC (body-centered cubic)";
            parameter = 2;
        } else {
            lattice_type = "SC (simple cubic)";
            parameter = 1;
        }
    } else {
        lattice_type = "uncertain";
        parameter = 1;
    }

    for (int i = 0; i < peaks.size(); i++) {
        double index_d = peaks.at(i).d2 * parameter / peaks.at(0).d2;
        peaks.at(i).index = round(index_d);
        peaks.at(i).a = sqrt(peaks.at(i).index / peaks.at(i).d2);
    }

    return;
}



void refine(const vector<peakdata>& peaks, double& lattice_const, double& sigma) {
    double A = 0;
    double A2 = 0;
    double NR = 0;
    double NR2 = 0;
    double NRA = 0;
    for (int i = 0; i < peaks.size(); i++) {
        A += peaks.at(i).a;
        A2 += pow(peaks.at(i).a, 2);
        NR += peaks.at(i).nr;
        NR2 += pow(peaks.at(i).nr, 2);
        NRA += (peaks.at(i).nr * peaks.at(i).a);
    }

    double N = peaks.size();
    double intercept = (NR2*A - NRA*NR) / (N*NR2 - pow(NR, 2));
    double slope = (N*NRA - NR*A) / (N*NR2 - pow(NR, 2));
    double sigma0 = sqrt(NR2 / (N*NR2 - pow(NR, 2)));
    double sigmai = sqrt((A2 - N*pow(intercept, 2) - NR2*pow(slope, 2) - 2*intercept*slope*NR) / (N-2));

    lattice_const = intercept;
    sigma = sigma0 * sigmai;
    return;
}



void addlog(const string ifname, const string ofname, const string lattice_type, const double lattice_const, const double sigma, const string lfname) {
    int flag_first = 0;
    if (!exists(lfname)) {
        flag_first = 1;
    }
    ofstream logfile(lfname, std::ios::app);
    if (flag_first) {
        logfile << "input_file : output_file : lattice_type : lattice_const : error" << endl;
    }
    logfile << ifname << " : ";
    logfile << ofname << " : ";
    logfile << lattice_type << " : ";
    logfile << lattice_const << " : ";
    logfile << sigma << endl;
    logfile.close();
    return;
}
