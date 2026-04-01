#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <cctype>
using namespace std;



string strip(const string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

vector<string> splitCSV(const string& line) {
    vector<string> tokens;
    string token;
    stringstream ss(line);
    while (getline(ss, token, ',')) {
        string t = strip(token);
        if (!t.empty()) tokens.push_back(t);
    }
    return tokens;
}


bool safeParseDouble(const string& s, double& out) {
    try {
        size_t idx = 0;
        out = stod(s, &idx);

        return idx > 0;
    }
    catch (...) {
        return false;
    }
}


bool findFirstNumberInLine(const string& line, double& val) {
    for (size_t i = 0; i < line.size(); ++i) {
        if ((line[i] >= '0' && line[i] <= '9') || line[i] == '-' || line[i] == '+' || line[i] == '.') {
            string tail = strip(line.substr(i));
            return safeParseDouble(tail, val);
        }
    }
    return false;
}



class GaussQuadrature {
public:
    int N;
    vector<double> xi;
    vector<double> w;
    GaussQuadrature(int n = 2) : N(n) {
        xi.clear(); w.clear();
        switch (n) {
        case 1:
            xi = { 0.0 }; w = { 2.0 };
            break;
        case 2: {
            double a = sqrt(1.0 / 3.0);
            xi = { -a, a }; w = { 1.0, 1.0 };
            break;
        }
        case 3: {
            double a = sqrt(3.0 / 5.0);
            xi = { -a, 0.0, a };
            w = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
            break;
        }
        case 4: {
            double a = sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0));
            double b = sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0));
            double w1 = (18.0 - sqrt(30.0)) / 36.0;
            double w2 = (18.0 + sqrt(30.0)) / 36.0;
            xi = { -b, -a, a, b };
            w = { w1, w2, w2, w1 };
            break;
        }
        case 5: {
            double a = (1.0 / 3.0) * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0));
            double b = (1.0 / 3.0) * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0));
            double w1 = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
            double w2 = 128.0 / 225.0;
            double w3 = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
            xi = { -b, -a, 0.0, a, b };
            w = { w1, w3, w2, w3, w1 };
            break;
        }
        default:
            cerr << "Nieznana liczba punktow Gaussa: " << n << "\n";
            break;
        }
    }
};


class ElemUniv {
public:
    int npc;
    GaussQuadrature quad;
    ElemUniv(int n = 2) : npc(n), quad(n) {}
};

class Node {
public:
    int id;
    double x;
    double y;
    Node(int _id = 0, double _x = 0.0, double _y = 0.0) : id(_id), x(_x), y(_y) {}
};

class Grid {
public:
    vector<Node> nodeList;
    vector<int> boundaryList;
    int nN = 0;
    int nE = 0;
    vector<struct Element> elems; 
};

class Jakobian {
public:
    double J[2][2];
    double invJ[2][2];
    double det;
    double dN_dx[4];
    double dN_dy[4];

    Jakobian() {
        det = 0.0;
        for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) { J[i][j] = invJ[i][j] = 0.0; }
        for (int k = 0; k < 4; k++) { dN_dx[k] = dN_dy[k] = 0.0; }
    }

    Jakobian(const array<double, 4>& xcoords, const array<double, 4>& ycoords,
        const array<double, 4>& dN_dE, const array<double, 4>& dN_dN)
    {

        J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
        for (int i = 0; i < 4; i++) {
            J[0][0] += dN_dE[i] * xcoords[i];
            J[0][1] += dN_dE[i] * ycoords[i];
            J[1][0] += dN_dN[i] * xcoords[i];
            J[1][1] += dN_dN[i] * ycoords[i];
        }
        det = J[0][0] * J[1][1] - J[1][0] * J[0][1];
        if (fabs(det) < 1e-14) {

        }
        invJ[0][0] = J[1][1] / det;
        invJ[0][1] = -J[0][1] / det;
        invJ[1][0] = -J[1][0] / det;
        invJ[1][1] = J[0][0] / det;

        for (int i = 0; i < 4; i++) {
            dN_dx[i] = invJ[0][0] * dN_dE[i] + invJ[0][1] * dN_dN[i];
            dN_dy[i] = invJ[1][0] * dN_dE[i] + invJ[1][1] * dN_dN[i];
        }
    }

    void dump() const {
        cout << "J = [" << J[0][0] << ", " << J[0][1] << "]\n";
        cout << "    [" << J[1][0] << ", " << J[1][1] << "]\n";
        cout << "detJ = " << det << "\n";
        cout << "invJ = [" << invJ[0][0] << ", " << invJ[0][1] << "]\n";
        cout << "       [" << invJ[1][0] << ", " << invJ[1][1] << "]\n";
        for (int i = 0; i < 4; i++) {
            cout << "dN_dx[" << i << "]=" << dN_dx[i] << ", dN_dy[" << i << "]=" << dN_dy[i] << "\n";
        }
    }
};

struct Element {
    int id = 0;
    int nodeIDs[4] = { 0,0,0,0 };
    vector<Jakobian> jacobians; 
    double Hmat[4][4];

    Element() {
        for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) Hmat[i][j] = 0.0;
    }

    void computeH(const ElemUniv& eU, const Grid& g, double conductivity,
        const vector<array<double, 4>>& naturalDerivativesE,
        const vector<array<double, 4>>& naturalDerivativesN)
    {
       
        array<double, 4> xs, ys;
        for (int i = 0; i < 4; i++) {
            int nid = nodeIDs[i] - 1; 
            xs[i] = g.nodeList[nid].x;
            ys[i] = g.nodeList[nid].y;
        }

        int np = eU.npc * eU.npc;
      
        if ((int)jacobians.size() != np) jacobians.resize(np);

        
        for (int p = 0; p < np; p++) {
          
            const array<double, 4>& dN_dE_p = naturalDerivativesE[p];
            const array<double, 4>& dN_dN_p = naturalDerivativesN[p];

            Jakobian J(xs, ys, dN_dE_p, dN_dN_p);
            jacobians[p] = J;

       
            int i = p % eU.npc;
            int j = p / eU.npc; 
            double wx = eU.quad.w[i];
            double wy = eU.quad.w[j];
            double weight = wx * wy;

            for (int a = 0; a < 4; a++) {
                for (int b = 0; b < 4; b++) {
                    Hmat[a][b] += conductivity * (J.dN_dx[a] * J.dN_dx[b] + J.dN_dy[a] * J.dN_dy[b]) * J.det * weight;
                }
            }
        }
    }

    void printH() const {
        cout << "\nMacierz H elementu " << id << ":\n";
        cout << setw(10) << " " << setw(12) << "N1" << setw(12) << "N2" << setw(12) << "N3" << setw(12) << "N4" << "\n";
        cout << string(60, '-') << "\n";
        for (int r = 0; r < 4; r++) {
            cout << "N" << (r + 1) << " ";
            for (int c = 0; c < 4; c++) {
                cout << setw(12) << fixed << setprecision(6) << Hmat[r][c];
            }
            cout << "\n";
        }
    }
};

class GlobalData {
public:
    double SimulationTime = 0.0;
    double SimulationStepTime = 0.0;
    double Conductivity = 0.0;
    double Alfa = 0.0;
    double Tot = 0.0;
    double InitialTemp = 0.0;
    double Density = 0.0;
    double SpecificHeat = 0.0;
    int nN = 0;
    int nE = 0;
    int npc = 2; 
};


void shapeFunctionsQuad(double ksi, double eta, array<double, 4>& N) {
    N[0] = 0.25 * (1.0 - ksi) * (1.0 - eta);
    N[1] = 0.25 * (1.0 + ksi) * (1.0 - eta);
    N[2] = 0.25 * (1.0 + ksi) * (1.0 + eta);
    N[3] = 0.25 * (1.0 - ksi) * (1.0 + eta);
}

void shapeDerivativesNatural(double ksi, double eta, array<double, 4>& dN_dE, array<double, 4>& dN_dN) {
    dN_dE[0] = -0.25 * (1.0 - eta);
    dN_dE[1] = 0.25 * (1.0 - eta);
    dN_dE[2] = 0.25 * (1.0 + eta);
    dN_dE[3] = -0.25 * (1.0 + eta);

    dN_dN[0] = -0.25 * (1.0 - ksi);
    dN_dN[1] = -0.25 * (1.0 + ksi);
    dN_dN[2] = 0.25 * (1.0 + ksi);
    dN_dN[3] = 0.25 * (1.0 - ksi);
}


void prepareNaturalDerivatives(const ElemUniv& eU,
    vector<array<double, 4>>& derivativesE,
    vector<array<double, 4>>& derivativesN,
    vector<array<double, 4>>& shapeFuncs)
{
    int npc = eU.npc;
    int np = npc * npc;
    derivativesE.assign(np, array<double, 4>{0, 0, 0, 0});
    derivativesN.assign(np, array<double, 4>{0, 0, 0, 0});
    shapeFuncs.assign(np, array<double, 4>{0, 0, 0, 0});

    int counter = 0;
    for (int j = 0; j < npc; ++j) {
        for (int i = 0; i < npc; ++i) {
            double ksi = eU.quad.xi[i];
            double eta = eU.quad.xi[j];

            array<double, 4> dE, dN, N;
            shapeDerivativesNatural(ksi, eta, dE, dN);
            shapeFunctionsQuad(ksi, eta, N);

            derivativesE[counter] = dE;
            derivativesN[counter] = dN;
            shapeFuncs[counter] = N;
            ++counter;
        }
    }
}



bool loadFromFile(const string& filename, GlobalData& G, Grid& grid) {
    ifstream in(filename);
    if (!in.is_open()) {
        cerr << "Nie mozna otworzyc pliku: " << filename << "\n";
        return false;
    }

    string line;
    bool inNodes = false, inElems = false, inBC = false;

    while (getline(in, line)) {
        string t = strip(line);
        if (t.empty()) continue;

       
        if (t.find("SimulationTime") != string::npos) {
            double v; if (findFirstNumberInLine(t, v)) G.SimulationTime = v; continue;
        }
        if (t.find("SimulationStepTime") != string::npos) {
            double v; if (findFirstNumberInLine(t, v)) G.SimulationStepTime = v; continue;
        }
        if (t.find("Conductivity") != string::npos) {
            double v; if (findFirstNumberInLine(t, v)) G.Conductivity = v; continue;
        }
        if (t.find("Alfa") != string::npos) {
            double v; if (findFirstNumberInLine(t, v)) G.Alfa = v; continue;
        }
        if (t.find("Tot") != string::npos) {
            double v; if (findFirstNumberInLine(t, v)) G.Tot = v; continue;
        }
        if (t.find("InitialTemp") != string::npos) {
            double v; if (findFirstNumberInLine(t, v)) G.InitialTemp = v; continue;
        }
        if (t.find("Density") != string::npos) {
            double v; if (findFirstNumberInLine(t, v)) G.Density = v; continue;
        }
        if (t.find("SpecificHeat") != string::npos) {
            double v; if (findFirstNumberInLine(t, v)) G.SpecificHeat = v; continue;
        }
        if (t.find("Nodes number") != string::npos) {
            double v; if (findFirstNumberInLine(t, v)) G.nN = static_cast<int>(v); continue;
        }
        if (t.find("Elements number") != string::npos) {
            double v; if (findFirstNumberInLine(t, v)) G.nE = static_cast<int>(v); continue;
        }

       
        if (t.find("*Node") != string::npos) { inNodes = true; inElems = false; inBC = false; continue; }
        if (t.find("*Element") != string::npos) { inElems = true; inNodes = false; inBC = false; continue; }
        if (t.find("*BC") != string::npos) { inBC = true; inNodes = false; inElems = false; continue; }

        if (inNodes) {
            auto parts = splitCSV(t);
            if (parts.size() < 3) {
                cerr << "Warning: zla linia w *Node: \"" << t << "\"\n";
                continue;
            }
            double xd, yd; int id;
            try {
                id = stoi(parts[0]);
                if (!safeParseDouble(parts[1], xd) || !safeParseDouble(parts[2], yd)) throw runtime_error("pars error");
                grid.nodeList.emplace_back(id, xd, yd);
            }
            catch (...) {
                cerr << "Warning: blad parsowania wezla: \"" << t << "\"\n";
                continue;
            }
        }
        if (inElems) {
            auto parts = splitCSV(t);
            if (parts.size() < 5) {
                cerr << "Warning: zla linia w *Element: \"" << t << "\"\n";
                continue;
            }
            try {
                Element el;
                el.id = stoi(parts[0]);
                for (int k = 0; k < 4; k++) el.nodeIDs[k] = stoi(parts[k + 1]);
                grid.elems.push_back(el);
            }
            catch (...) {
                cerr << "Warning: blad parsowania elementu: \"" << t << "\"\n";
                continue;
            }
        }
        if (inBC) {
            auto parts = splitCSV(t);
            for (auto& p : parts) {
                try {
                    int nid = stoi(p);
                    grid.boundaryList.push_back(nid);
                }
                catch (...) {
                    cerr << "Warning: blad parsowania BC: \"" << t << "\"\n";
                }
            }
        }
    } 

    grid.nN = (int)grid.nodeList.size();
    grid.nE = (int)grid.elems.size();
    if (G.nN == 0) G.nN = grid.nN;
    if (G.nE == 0) G.nE = grid.nE;

    in.close();
    return true;
}



void printGridSummary(const GlobalData& G, const Grid& grid) {
    cout << "Dane globalne:\n";
    cout << " SimulationTime: " << G.SimulationTime << "\n";
    cout << " SimulationStepTime: " << G.SimulationStepTime << "\n";
    cout << " Conductivity: " << G.Conductivity << "\n";
    cout << " Alfa: " << G.Alfa << "\n";
    cout << " Tot: " << G.Tot << "\n";
    cout << " InitialTemp: " << G.InitialTemp << "\n";
    cout << " Density: " << G.Density << "\n";
    cout << " SpecificHeat: " << G.SpecificHeat << "\n";
    cout << " Nodes (global): " << G.nN << "  Elements (global): " << G.nE << "\n\n";

    cout << "Wspolrzedne wezlow:\n";
    for (const auto& n : grid.nodeList) {
        cout << " ID: " << n.id << "  x=" << n.x << "  y=" << n.y << "\n";
    }
    cout << "\nElementy:\n";
    for (const auto& e : grid.elems) {
        cout << " Element ID: " << e.id << "  wezly: ";
        for (int k = 0; k < 4; k++) {
            cout << e.nodeIDs[k] << (k < 3 ? ", " : "");
        }
        cout << "\n";
    }
    cout << "\nWezly BC:\n";
    for (int b : grid.boundaryList) cout << b << " ";
    cout << "\n";
}



int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    GlobalData G;
    Grid grid;


    vector<string> files = { "Test1_4_4.txt", "Test2_4_4MixGrid.txt", "Test3_31_31_kwadrat.txt", "Test4_testowe.txt" };

    if (!loadFromFile(files[3], G, grid)) {
        return 1;
    }

    printGridSummary(G, grid);

 
    ElemUniv eU(G.npc);


    vector<array<double, 4>> dN_dE_all, dN_dN_all, N_all;
    prepareNaturalDerivatives(eU, dN_dE_all, dN_dN_all, N_all);

    cout << fixed << setprecision(6);
    for (auto& el : grid.elems) {

        cout << "\nElement " << el.id << "\n";
        cout << "Node coords (N1..N4):\n";
        for (int k = 0; k < 4; k++) {
            int nid = el.nodeIDs[k] - 1;
            cout << " N" << (k + 1) << " id=" << el.nodeIDs[k]
                << " (" << grid.nodeList[nid].x << ", " << grid.nodeList[nid].y << ")\n";
        }

        el.computeH(eU, grid, G.Conductivity, dN_dE_all, dN_dN_all);


        int np = eU.npc * eU.npc;
        for (int p = 0; p < np; p++) {
            double ksi = eU.quad.xi[p % eU.npc];
            double eta = eU.quad.xi[p / eU.npc];
            cout << "\npc=" << (p + 1) << " ksi=" << ksi << " eta=" << eta << "\n";
            cout << " dN_dE: ";
            for (int a = 0; a < 4; a++) cout << setw(10) << dN_dE_all[p][a];
            cout << "\n dN_dN: ";
            for (int a = 0; a < 4; a++) cout << setw(10) << dN_dN_all[p][a];
            cout << "\n";

            const Jakobian& J = el.jacobians[p];
            cout << " J:\n";
            cout << "  [" << J.J[0][0] << ", " << J.J[0][1] << "]\n";
            cout << "  [" << J.J[1][0] << ", " << J.J[1][1] << "]\n";
            cout << " detJ = " << J.det << "\n";
            cout << " J^-1:\n";
            cout << "  [" << J.invJ[0][0] << ", " << J.invJ[0][1] << "]\n";
            cout << "  [" << J.invJ[1][0] << ", " << J.invJ[1][1] << "]\n";
            cout << " dN_dx: ";
            for (int a = 0; a < 4; a++) cout << setw(10) << J.dN_dx[a];
            cout << "\n dN_dy: ";
            for (int a = 0; a < 4; a++) cout << setw(10) << J.dN_dy[a];
            cout << "\n";

          
            double manual = J.invJ[0][0] * dN_dE_all[p][0] + J.invJ[0][1] * dN_dN_all[p][0];
            cout << " manual dN_dx[0] = " << manual << "\n";
        }


        el.printH();
    }

    return 0;
}
