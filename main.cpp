#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <map>
#include <cmath>

using namespace std;


class GaussQuadrature {
public:
    int N;
    vector<double> xi;
    vector<double> w;

    GaussQuadrature(int n) {
        N = n;
        switch (n) {
        case 1: {
            xi = { 0.0 };
            w = { 2.0 };
            break;
        }
        case 2: {
            double a = sqrt(1.0 / 3.0);
            xi = { -a, a };
            w = { 1.0, 1.0 };
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
            cerr << "Nieobslugiwany schemat Gaussa: N=" << n << endl;
            break;
        }
    }
};


double f1(double x) {
    return 5 * x * x + 3 * x + 6;
}

double f2(double x, double y) {
    return 5 * x * x * y * y + 3 * x * y + 6;
}


//  1D 

double gauss1D(double (*f)(double), const GaussQuadrature& G) {
    double sum = 0.0;
    for (int i = 0; i < G.N; ++i)
        sum += G.w[i] * f(G.xi[i]);
    return sum;
}


//  2D

double gauss2D(double (*f)(double, double), const GaussQuadrature& G) {
    double sum = 0.0;
    for (int i = 0; i < G.N; ++i)
        for (int j = 0; j < G.N; ++j)
            sum += G.w[i] * G.w[j] * f(G.xi[i], G.xi[j]);
    return sum;
}


string trim(const string& str) {
    auto start = str.find_first_not_of(" \t\r\n");
    auto end = str.find_last_not_of(" \t\r\n");
    return (start == string::npos) ? "" : str.substr(start, end - start + 1);
}


vector<string> splitAndTrim(const string& s, char delim) {
    vector<string> parts;
    string token;
    stringstream ss(s);
    while (getline(ss, token, delim)) {
        string t = trim(token);
        if (!t.empty()) parts.push_back(t);
    }
    return parts;
}

class Node;
class Element;
class Grid;
class ElemUniv;
class Jakobian;


class Node {
public:
    int id;
    double x, y;
    int BC;
    Node(int id = 0, double x = 0, double y = 0, int BC = 0) : id(id), x(x), y(y), BC(BC) {}
};

class ElemUniv {
public:
    int npc;
    GaussQuadrature gauss;

    vector<vector<double>> dN_dE;
    vector<vector<double>> dN_dN;
    vector<vector<double>> N;

    ElemUniv(int npc) : npc(npc), gauss(npc) {

        int punktyCalkowania = npc * npc;

        N.resize(punktyCalkowania, vector<double>(4));
        dN_dE.resize(punktyCalkowania, vector<double>(4));
        dN_dN.resize(punktyCalkowania, vector<double>(4));

        int licznik = 0;

        for (int j = 0; j < npc; j++)
        {
            for (int i = 0; i < npc; i++)
            {


                double ksi = gauss.xi[i];
                double eta = gauss.xi[j];

                //  funkcje ksztaltu 
                N[licznik][0] = 0.25 * (1 - ksi) * (1 - eta);
                N[licznik][1] = 0.25 * (1 + ksi) * (1 - eta);
                N[licznik][2] = 0.25 * (1 + ksi) * (1 + eta);
                N[licznik][3] = 0.25 * (1 - ksi) * (1 + eta);

                //  pochodne wzglêdem E 
                dN_dE[licznik][0] = -0.25 * (1 - eta);
                dN_dE[licznik][1] = 0.25 * (1 - eta);
                dN_dE[licznik][2] = 0.25 * (1 + eta);
                dN_dE[licznik][3] = -0.25 * (1 + eta);

                //  pochodne wzglêdem N 
                dN_dN[licznik][0] = -0.25 * (1 - ksi);
                dN_dN[licznik][1] = -0.25 * (1 + ksi);
                dN_dN[licznik][2] = 0.25 * (1 + ksi);
                dN_dN[licznik][3] = 0.25 * (1 - ksi);

                licznik++;
            }

        }

    }




};

class Surface {
public:
    int npc;
    GaussQuadrature gauss;
    vector<vector<vector<double>>> N;

    Surface(int npc) : npc(npc), gauss(npc)
    {
        N.resize(4);
        for (int e = 0; e < 4; e++)
            N[e].resize(npc, vector<double>(4, 0.0));

        for (int p = 0; p < npc; p++)
        {
            double s = gauss.xi[p];


            N[0][p][0] = 0.5 * (1 - s);
            N[0][p][1] = 0.5 * (1 + s);

            N[1][p][1] = 0.5 * (1 - s);
            N[1][p][2] = 0.5 * (1 + s);

            N[2][p][2] = 0.5 * (1 - s);
            N[2][p][3] = 0.5 * (1 + s);

            N[3][p][3] = 0.5 * (1 - s);
            N[3][p][0] = 0.5 * (1 + s);
        }
    }
};



class Jakobian {
public:
    double J[2][2];
    double Jodwr[2][2];
    double detJ;

    double dN_dx[4], dN_dy[4];


    Jakobian() {
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j) {
                J[i][j] = 0.0;
                Jodwr[i][j] = 0.0;
            }
        detJ = 0.0;
        for (int i = 0; i < 4; ++i) {
            dN_dx[i] = 0.0;
            dN_dy[i] = 0.0;
        }
    }

    Jakobian(const Element& elem, const Grid& grid, const ElemUniv& eUniv, int p);



    void coutJakobian() {

        std::cout << "Macierz J:\n";
        std::cout << "[" << J[0][0] << ", " << J[0][1] << "]\n";
        std::cout << "[" << J[1][0] << ", " << J[1][1] << "]\n";
        std::cout << "Wyznacznik detJ: " << detJ << "\n";

        std::cout << "Macierz odwrotna J^-1:\n";
        std::cout << "[" << Jodwr[0][0] << ", " << Jodwr[0][1] << "]\n";
        std::cout << "[" << Jodwr[1][0] << ", " << Jodwr[1][1] << "]\n";

        std::cout << "Pochodne funkcji ksztaltu w uk³adzie globalnym:\n";
        for (int i = 0; i < 4; ++i) {
            std::cout << "dN_dx[" << i << "] = " << dN_dx[i]
                << ", dN_dy[" << i << "] = " << dN_dy[i] << "\n";
        }

        std::cout << endl;
    }


    void samJakobian() {
        std::cout << "[" << J[0][0] << ", " << J[0][1] << "]\n";
        std::cout << "[" << J[1][0] << ", " << J[1][1] << "]\n";
    }
};

struct Material {
    int id;
    std::string name;
    double conductivity;
    double density;
    double specificHeat;
    Material(int id_ = 0, const std::string& name_ = "", double k = 0.0, double rho = 0.0, double cp = 0.0)
        : id(id_), name(name_), conductivity(k), density(rho), specificHeat(cp) {
    }
};



class Element {
public:
    int id;
    int ID[4];
    vector<Jakobian> jakobiany;
    vector<vector<double>> H;
    vector<vector<double>> Hbc;
    vector<double> P;
    vector<vector<double>> C;
    int materialId = 0;

    void obliczH(const ElemUniv& eUniv, const Grid& grid, const Material& M) {
        H.assign(4, vector<double>(4, 0.0));
        C.assign(4, vector<double>(4, 0.0));

        double conductivity = M.conductivity;
        double density = M.density;
        double specificHeat = M.specificHeat;

        int npc = eUniv.npc;


        vector<vector<double>> C_consistent(4, vector<double>(4, 0.0));

        for (int p = 0; p < npc * npc; p++) {
            Jakobian& J = jakobiany[p];

            int col = p % npc;
            int row = p / npc;

            double w_x = eUniv.gauss.w[col];
            double w_y = eUniv.gauss.w[row];

            const vector<double>& fKszt = eUniv.N[p];

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    H[i][j] += conductivity * (J.dN_dx[i] * J.dN_dx[j] + J.dN_dy[i] * J.dN_dy[j]) * J.detJ * w_x * w_y;
                    C_consistent[i][j] += specificHeat * density * (fKszt[i] * fKszt[j]) * J.detJ * w_x * w_y;
                }
            }
        }

        // LUMPED MASS MATRIX 

        for (int i = 0; i < 4; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < 4; j++) {
                rowSum += C_consistent[i][j];
            }
            C[i][i] = rowSum;
        }
    }




    void obliczHbc(const Surface& surface, const Grid& grid, double alfaOut, double alfaIn);
    void obliczP(const Surface& surface, const Grid& grid, double alfaOut, double alfaIn, double TotOut, double TotIn);
};



class Grid {
public:
    int nN = 0;
    int nE = 0;
    vector<Node> nodes;
    vector<Element> elements;
    vector<int> boundaryNodes;
};


void Element::obliczHbc(const Surface& surface, const Grid& grid, double alfaOut, double alfaIn) {
    Hbc.assign(4, vector<double>(4, 0.0));
    const int edge[4][2] = { {0,1}, {1,2}, {2,3}, {3,0} };
    int npc = surface.npc;

    for (int e = 0; e < 4; e++)
    {
        int n1 = ID[edge[e][0]] - 1;
        int n2 = ID[edge[e][1]] - 1;

        if (grid.nodes[n1].BC == 0 || grid.nodes[n2].BC == 0)
            continue;

        double alfaEdge;
        int bcType = 0;

        if (grid.nodes[n1].BC == 1 && grid.nodes[n2].BC == 1) {
            alfaEdge = alfaOut;
            bcType = 1;
        }
        else if (grid.nodes[n1].BC == 2 && grid.nodes[n2].BC == 2) {
            alfaEdge = alfaIn;
            bcType = 2;
        }
        else {
            continue;
        }


        double x1 = grid.nodes[n1].x, y1 = grid.nodes[n1].y;
        double x2 = grid.nodes[n2].x, y2 = grid.nodes[n2].y;

        double dl = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
        double detJ = dl / 2.0;

        for (int p = 0; p < npc; p++)
        {
            const vector<double>& fKsz = surface.N[e][p];
            double w = surface.gauss.w[p];

            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    Hbc[i][j] += alfaEdge * fKsz[i] * fKsz[j] * detJ * w;
        }
    }
}

void Element::obliczP(const Surface& surface, const Grid& grid, double alfaOut, double alfaIn, double TotOut, double TotIn) {

    P.assign(4, 0.0);
    const int edge[4][2] = { {0,1}, {1,2}, {2,3}, {3,0} };

    int npc = surface.npc;



    for (int e = 0; e < 4; e++)
    {
        int n1 = ID[edge[e][0]] - 1;
        int n2 = ID[edge[e][1]] - 1;


        if (grid.nodes[n1].BC == 0 || grid.nodes[n2].BC == 0)
        {
            continue;
        }


        double alfaEdge;
        double totEdge;

        if (grid.nodes[n1].BC == 1 && grid.nodes[n2].BC == 1) {
            alfaEdge = alfaOut;
            totEdge = TotOut;

        }
        else if (grid.nodes[n1].BC == 2 && grid.nodes[n2].BC == 2) {
            alfaEdge = alfaIn;
            totEdge = TotIn;
        }
        else {
            continue;
        }

        double x1 = grid.nodes[n1].x, y1 = grid.nodes[n1].y;
        double x2 = grid.nodes[n2].x, y2 = grid.nodes[n2].y;

        double dl = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
        double detJ = dl / 2.0;


        for (int p = 0; p < npc; p++)
        {
            const vector<double>& fKsz = surface.N[e][p];

            double w = surface.gauss.w[p];

            for (int i = 0; i < 4; i++)

                P[i] += totEdge * alfaEdge * fKsz[i] * detJ * w;
        }



    }
}



class GlobalData {
public:
    double SimulationTime = 0;
    double SimulationStepTime = 0;
    double AlfaOut = 0;
    double AlfaIn = 0;
    double TotOut = 0;
    double TotIn = 0;
    double InitialTemp = 0;
    int nN = 0;
    int nE = 0;

    std::vector<Material> materials;
    std::map<int, size_t> materialIdToIndex;

    //ustaw
    int npc = 4;
};




Jakobian::Jakobian(const Element& elem, const Grid& grid, const ElemUniv& eUniv, int p) {


    double x[4], y[4];
    for (int i = 0; i < 4; i++) {
        int id = elem.ID[i] - 1;       // -1 bo w pliku wezly liczone od 1
        x[i] = grid.nodes[id].x;
        y[i] = grid.nodes[id].y;
    }
    J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
    for (int i = 0; i < 4; i++)
    {
        J[0][0] += eUniv.dN_dE[p][i] * x[i];                   // x/E
        J[0][1] += eUniv.dN_dE[p][i] * y[i];                   // y/E
        J[1][0] += eUniv.dN_dN[p][i] * x[i];                   // x/n
        J[1][1] += eUniv.dN_dN[p][i] * y[i];                   // y/n

    }

    detJ = (J[0][0] * J[1][1]) - (J[1][0] * J[0][1]);

    if (detJ <= 0) {
        cerr << "Zly Jacobian w elemencie " << elem.id << endl;
        exit(1);
    }

    Jodwr[0][0] = J[1][1] / detJ;
    Jodwr[0][1] = -J[0][1] / detJ;
    Jodwr[1][0] = -J[1][0] / detJ;
    Jodwr[1][1] = J[0][0] / detJ;

    //fksz globalna
    for (int i = 0; i < 4; i++) {
        dN_dx[i] = Jodwr[0][0] * eUniv.dN_dE[p][i] + Jodwr[0][1] * eUniv.dN_dN[p][i];
        dN_dy[i] = Jodwr[1][0] * eUniv.dN_dE[p][i] + Jodwr[1][1] * eUniv.dN_dN[p][i];
    }

}



void printMatrix(const vector<vector<double>>& M, const string& title) {
    cout << "\n" << title << "\n";
    for (size_t i = 0; i < M.size(); ++i) {
        for (size_t j = 0; j < M[i].size(); ++j) {
            cout << setw(12) << fixed << setprecision(6) << M[i][j];
        }
        cout << "\n";
    }
}




// pomocnik do znalezienia i sparsowania liczby po kluczu
bool parseNumberAfterKey(const string& line, double& out) {
    // znajd pierwsz¹ cyfrê, znak minus, lub kropkê
    size_t pos = line.find_first_of("0123456789-+.");
    if (pos == string::npos) return false;
    string numStr = trim(line.substr(pos));
    if (numStr.empty()) return false;
    try {
        out = stod(numStr);
        return true;
    }
    catch (...) {
        return false;
    }
}


bool isCCW(double x[4], double y[4]) {

    double area = 0.0;
    for (int i = 0; i < 4; i++) {
        int j = (i + 1) % 4;
        area += x[i] * y[j];
        area -= x[j] * y[i];
    }
    return area > 0;
}

bool loadFromFile(const string& filename, GlobalData& globalData, Grid& grid) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Nie mo¿na otworzyæ pliku: " << filename << endl;
        return false;
    }

    string line;
    bool readingNodes = false;
    bool readingElements = false;
    bool readingBC = false;
    bool readingBCOut = false;
    bool readingBCIn = false;
    bool readingMaterials = false;
    bool readingElementMat = false;

    vector<int> bcOutList;
    vector<int> bcInList;
    vector<int> bcCombinedList;

    while (getline(file, line)) {
        string trimmedLine = trim(line);
        if (trimmedLine.empty()) {
            continue;
        }

        if (trimmedLine.find("SimulationTime") != string::npos) {
            double v; if (parseNumberAfterKey(trimmedLine, v)) globalData.SimulationTime = v;
            readingNodes = readingElements = readingBC = readingBCOut = readingBCIn = readingMaterials = readingElementMat = false;
            continue;
        }
        if (trimmedLine.find("SimulationStepTime") != string::npos) {
            double v; if (parseNumberAfterKey(trimmedLine, v)) globalData.SimulationStepTime = v;
            readingNodes = readingElements = readingBC = readingBCOut = readingBCIn = readingMaterials = readingElementMat = false;
            continue;
        }

        if (trimmedLine.find("Alfa_out") != string::npos || trimmedLine.find("AlfaOut") != string::npos) {
            double v; if (parseNumberAfterKey(trimmedLine, v)) globalData.AlfaOut = v;
            readingNodes = readingElements = readingBC = readingBCOut = readingBCIn = readingMaterials = readingElementMat = false;
            continue;
        }
        if (trimmedLine.find("Alfa_in") != string::npos || trimmedLine.find("AlfaIn") != string::npos) {
            double v; if (parseNumberAfterKey(trimmedLine, v)) globalData.AlfaIn = v;
            readingNodes = readingElements = readingBC = readingBCOut = readingBCIn = readingMaterials = readingElementMat = false;
            continue;
        }

        if (trimmedLine.find("Tot_out") != string::npos || trimmedLine.find("TotOut") != string::npos) {
            double v; if (parseNumberAfterKey(trimmedLine, v)) globalData.TotOut = v;
            readingNodes = readingElements = readingBC = readingBCOut = readingBCIn = readingMaterials = readingElementMat = false;
            continue;
        }
        if (trimmedLine.find("Tot_in") != string::npos || trimmedLine.find("TotIn") != string::npos) {
            double v; if (parseNumberAfterKey(trimmedLine, v)) globalData.TotIn = v;
            readingNodes = readingElements = readingBC = readingBCOut = readingBCIn = readingMaterials = readingElementMat = false;
            continue;
        }
        if (trimmedLine.find("InitialTemp") != string::npos) {
            double v; if (parseNumberAfterKey(trimmedLine, v)) globalData.InitialTemp = v;
            readingNodes = readingElements = readingBC = readingBCOut = readingBCIn = readingMaterials = readingElementMat = false;
            continue;
        }
        if (trimmedLine.find("Nodes number") != string::npos) {
            double v; if (parseNumberAfterKey(trimmedLine, v)) globalData.nN = static_cast<int>(v);
            readingNodes = readingElements = readingBC = readingBCOut = readingBCIn = readingMaterials = readingElementMat = false;
            continue;
        }
        if (trimmedLine.find("Elements number") != string::npos) {
            double v; if (parseNumberAfterKey(trimmedLine, v)) globalData.nE = static_cast<int>(v);
            readingNodes = readingElements = readingBC = readingBCOut = readingBCIn = readingMaterials = readingElementMat = false;
            continue;
        }

        if (trimmedLine.size() >= 5 && trimmedLine.substr(0, 5) == "*Node") {
            readingNodes = true;
            readingElements = readingBC = readingBCOut = readingBCIn = readingMaterials = readingElementMat = false;
            continue;
        }
        if (trimmedLine.size() >= 8 && trimmedLine.substr(0, 8) == "*Element" && trimmedLine.find("Mat") == string::npos) {
            readingElements = true;
            readingNodes = readingBC = readingBCOut = readingBCIn = readingMaterials = readingElementMat = false;
            continue;
        }
        if (trimmedLine.find("*BC_OUT") != string::npos) {
            readingBCOut = true;
            readingNodes = readingElements = readingBC = readingBCIn = readingMaterials = readingElementMat = false;
            continue;
        }
        if (trimmedLine.find("*BC_IN") != string::npos) {
            readingBCIn = true;
            readingNodes = readingElements = readingBC = readingBCOut = readingMaterials = readingElementMat = false;
            continue;
        }
        if (trimmedLine.size() >= 3 && trimmedLine.substr(0, 3) == "*BC" && trimmedLine.find("_OUT") == string::npos && trimmedLine.find("_IN") == string::npos) {
            readingBC = true;
            readingNodes = readingElements = readingBCOut = readingBCIn = readingMaterials = readingElementMat = false;
            continue;
        }
        if (trimmedLine.find("*Materials") != string::npos) {
            readingMaterials = true;
            readingNodes = readingElements = readingBC = readingBCOut = readingBCIn = readingElementMat = false;
            continue;
        }
        if (trimmedLine.find("*ElementMat") != string::npos || trimmedLine.find("*ElementMaterial") != string::npos) {
            readingElementMat = true;
            readingMaterials = readingNodes = readingElements = readingBC = readingBCOut = readingBCIn = false;
            continue;
        }

        if (readingNodes) {
            auto parts = splitAndTrim(trimmedLine, ',');
            if (parts.size() < 3) {
                cerr << "Warning: skipping malformed node line: \"" << trimmedLine << "\"\n";
                continue;
            }
            try {
                int id = stoi(parts[0]);
                double x = stod(parts[1]);
                double y = stod(parts[2]);
                grid.nodes.emplace_back(id, x, y, 0);
            }
            catch (const exception& e) {
                cerr << "Warning: failed to parse node: \"" << trimmedLine << "\" (" << e.what() << ")\n";
            }
            continue;
        }

        if (readingElements) {
            auto parts = splitAndTrim(trimmedLine, ',');
            if (parts.size() < 5) {
                cerr << "Warning: skipping malformed element line: \"" << trimmedLine << "\"\n";
                continue;
            }
            try {
                Element e;
                e.id = stoi(parts[0]);
                for (int i = 0; i < 4; ++i) e.ID[i] = stoi(parts[i + 1]);
                e.materialId = 0;

                double x[4], y[4];
                for (int i = 0; i < 4; i++) {
                    int nid = e.ID[i] - 1;
                    x[i] = grid.nodes[nid].x;
                    y[i] = grid.nodes[nid].y;
                }

                if (!isCCW(x, y)) {
                    // swap 1 and 3 node
                    swap(e.ID[1], e.ID[3]);
                    cout << "Element " << e.id << " - zamieniono na CCW" << endl;
                }

                grid.elements.push_back(e);
            }
            catch (const exception& ex) {
                cerr << "Warning: failed to parse element: \"" << trimmedLine << "\" (" << ex.what() << ")\n";
            }
            continue;
        }

        if (readingBCOut) {
            auto parts = splitAndTrim(trimmedLine, ',');
            for (const auto& p : parts) {
                try {
                    int nodeID = stoi(p);
                    if (find(bcOutList.begin(), bcOutList.end(), nodeID) == bcOutList.end())
                        bcOutList.push_back(nodeID);

                    if (nodeID > 0 && nodeID <= static_cast<int>(grid.nodes.size())) {
                        int& bc = grid.nodes[nodeID - 1].BC;
                        if (bc != 0 && bc != 1) {
                            cerr << "ERROR: Node " << nodeID << " already has BC=" << bc
                                << " but is listed in *BC_OUT\n";
                            exit(1);
                        }
                        bc = 1; // OUT
                    }
                }
                catch (...) {
                    cerr << "Warning: bad *BC_OUT token: \"" << p << "\"\n";
                }
            }
            continue;
        }

        if (readingBCIn) {
            auto parts = splitAndTrim(trimmedLine, ',');
            for (const auto& p : parts) {
                try {
                    int nodeID = stoi(p);
                    if (find(bcInList.begin(), bcInList.end(), nodeID) == bcInList.end())
                        bcInList.push_back(nodeID);

                    if (nodeID > 0 && nodeID <= static_cast<int>(grid.nodes.size())) {
                        int& bc = grid.nodes[nodeID - 1].BC;
                        if (bc != 0 && bc != 2) {
                            cerr << "ERROR: Node " << nodeID << " already has BC=" << bc
                                << " but is listed in *BC_IN\n";
                            exit(1);
                        }
                        bc = 2; // IN
                    }
                }
                catch (...) {
                    cerr << "Warning: bad *BC_IN token: \"" << p << "\"\n";
                }
            }
            continue;
        }

        if (readingBC) {
            auto parts = splitAndTrim(trimmedLine, ',');
            for (const auto& p : parts) {
                try {
                    int nodeID = stoi(p);
                    if (find(bcCombinedList.begin(), bcCombinedList.end(), nodeID) == bcCombinedList.end())
                        bcCombinedList.push_back(nodeID);

                }
                catch (...) {
                    cerr << "Warning: bad *BC token: \"" << p << "\"\n";
                }
            }
            continue;
        }

        if (readingMaterials) {
            auto parts = splitAndTrim(trimmedLine, ',');
            if (parts.size() < 4) {
                cerr << "Warning: skipping malformed Materials line: \"" << trimmedLine << "\"\n";
                continue;
            }
            try {
                int mid = stoi(parts[0]);
                string name;
                double k = 0.0, rho = 0.0, cp = 0.0;
                if (parts.size() == 4) {
                    name = "mat" + to_string(mid);
                    k = stod(parts[1]);
                    rho = stod(parts[2]);
                    cp = stod(parts[3]);
                }
                else {
                    name = parts[1];
                    k = stod(parts[2]);
                    rho = stod(parts[3]);
                    cp = stod(parts[4]);
                }
                Material mat(mid, name, k, rho, cp);
                globalData.materialIdToIndex[mid] = globalData.materials.size();
                globalData.materials.push_back(mat);
            }
            catch (const exception& ex) {
                cerr << "Warning: failed to parse material: \"" << trimmedLine << "\" (" << ex.what() << ")\n";
            }
            continue;
        }

        if (readingElementMat) {
            auto parts = splitAndTrim(trimmedLine, ',');
            if (parts.size() < 2) {
                cerr << "Warning: skipping malformed ElementMat line: \"" << trimmedLine << "\"\n";
                continue;
            }
            try {
                int elemId = stoi(parts[0]);
                int matId = stoi(parts[1]);
                bool found = false;
                for (auto& el : grid.elements) {
                    if (el.id == elemId) {
                        el.materialId = matId;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    cerr << "Warning: Element id " << elemId << " not found for ElementMat mapping\n";
                }
            }
            catch (...) {
                cerr << "Warning: bad ElementMat tokens in line: \"" << trimmedLine << "\"\n";
            }
            continue;
        }
    }

    if (!bcCombinedList.empty()) {
        grid.boundaryNodes = bcCombinedList;
    }
    else if (!bcOutList.empty() && !bcInList.empty()) {
        grid.boundaryNodes.clear();
        size_t nOut = bcOutList.size();
        size_t nIn = bcInList.size();
        size_t nMax = max(nOut, nIn);
        for (size_t i = 0; i < nMax; ++i) {
            if (i < nOut) grid.boundaryNodes.push_back(bcOutList[i]);
            if (i < nIn)  grid.boundaryNodes.push_back(bcInList[i]);
        }
        vector<int> uniqueList;
        uniqueList.reserve(grid.boundaryNodes.size());
        for (int id : grid.boundaryNodes) {
            if (find(uniqueList.begin(), uniqueList.end(), id) == uniqueList.end())
                uniqueList.push_back(id);
        }
        grid.boundaryNodes.swap(uniqueList);
    }
    else {
        grid.boundaryNodes.clear();
        for (const auto& n : grid.nodes) {
            if (n.BC != 0) grid.boundaryNodes.push_back(n.id);
        }
    }

    grid.nN = static_cast<int>(grid.nodes.size());
    grid.nE = static_cast<int>(grid.elements.size());
    if (globalData.nN == 0) globalData.nN = grid.nN;
    if (globalData.nE == 0) globalData.nE = grid.nE;

    file.close();
    return true;
}




void printDerivativeTable(const vector<vector<double>>& derivatives, const string& label) {
    cout << "\nTabela: " << label << endl;

    // Nag³ówki kolumn
    cout << setw(6) << " "
        << setw(10) << "dN1"
        << setw(10) << "dN2"
        << setw(10) << "dN3"
        << setw(10) << "dN4" << endl;

    for (size_t i = 0; i < derivatives.size(); ++i) {
        cout << "pc" << (i + 1);
        for (size_t j = 0; j < derivatives[i].size(); ++j) {
            cout << setw(12) << fixed << setprecision(6) << derivatives[i][j];
        }
        cout << endl;
    }
}



vector<double> solveLinearSystem(vector<vector<double>> A, vector<double> b) {
    int N = A.size();
    vector<double> x(N, 0.0);




    for (int i = 0; i < N; ++i) {

        // pivot
        double maxEl = fabs(A[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < N; ++k) {
            if (fabs(A[k][i]) > maxEl) {
                maxEl = fabs(A[k][i]);
                maxRow = k;
            }
        }
        swap(A[i], A[maxRow]);
        swap(b[i], b[maxRow]);


        for (int k = i + 1; k < N; ++k) {
            double c = A[k][i] / A[i][i];
            for (int j = i; j < N; ++j)
                A[k][j] -= c * A[i][j];
            b[k] -= c * b[i];
        }
    }

    // podstawienie
    for (int i = N - 1; i >= 0; --i) {
        double sum = b[i];
        for (int j = i + 1; j < N; ++j)
            sum -= A[i][j] * x[j];
        x[i] = sum / A[i][i];
    }

    return x;
}


auto printWrapped = [&](const vector<int>& vec, const string& header, int per_line = 12) {
    std::cout << header << "\n";
    if (vec.empty()) {
        cout << " (brak)\n";
        return;
    }
    for (size_t i = 0; i < vec.size(); i += per_line) {
        size_t end = min(vec.size(), i + per_line);
        cout << " ";
        for (size_t j = i; j < end; ++j) {
            cout << vec[j];
            if (j + 1 < end) cout << ", ";
        }
        cout << "\n";
    }
    };







int main() {
    GlobalData globalData;
    Grid grid;

    vector<string> Pliki = { "Test1_4_4.txt", "Test2_4_4MixGrid.txt",
                             "Test3_31_31_kwadrat.txt", "Test4_testowe.txt", "wall.txt", "small.txt" };


    //zmien
    if (!loadFromFile(Pliki[5], globalData, grid)) return 1;

    //globalData.materials.push_back(Material(1, "tynk", 0.7, 1800.0, 840.0));
    //globalData.materialIdToIndex[1] = 0;
    //globalData.materials.push_back(Material(2, "styropian", 0.032, 250.0, 1400.0));
    //globalData.materialIdToIndex[2] = 1;
    //globalData.materials.push_back(Material(3, "cegla", 0.6, 1800.0, 840.0));
    //globalData.materialIdToIndex[3] = 2;



    cout << "TotOut=" << globalData.TotOut << " TotIn=" << globalData.TotIn << "\n";
    cout << "AlfaOut=" << globalData.AlfaOut << " AlfaIn=" << globalData.AlfaIn << "\n";
    for (int i = 0; i < globalData.materials.size(); i++)
    {
        cout << endl << globalData.materials[i].name << endl;
        cout << "conductivity " << globalData.materials[i].conductivity << endl;
        cout << "density " << globalData.materials[i].density << endl;
        cout << "specificHeat " << globalData.materials[i].specificHeat << endl << endl << endl;
    }

    for (auto& e : grid.elements) {
        if (e.materialId == 0) {
            cout << "\n\n\n\nelement nie ma materialu\n\n";
        }
    }


    vector<int> bc_out_nodes;
    vector<int> bc_in_nodes;

    for (const auto& n : grid.nodes) {
        if (n.BC == 1) bc_out_nodes.push_back(n.id);
        else if (n.BC == 2) bc_in_nodes.push_back(n.id);
    }


    cout << "Dane globalne:\n";
    cout << "SimulationTime: " << globalData.SimulationTime << "\n";
    cout << "SimulationStepTime: " << globalData.SimulationStepTime << "\n";
    cout << "Tot In: " << globalData.TotIn << "\n";
    cout << "Tot out: " << globalData.TotOut << "\n";
    cout << "InitialTemp: " << globalData.InitialTemp << "\n";
    cout << "Nodes (global): " << globalData.nN << "\tElements (global): " << globalData.nE << "\n\n";

    //  wezly i elementy
    cout << "Wspolrzedne wezlow:\n";
    for (const auto& n : grid.nodes)
        cout << "ID: " << n.id << "\tX: " << n.x << "\tY: " << n.y << "\n";

    cout << "\nElementy:\n";
    for (const auto& e : grid.elements) {
        cout << "Element ID: " << e.id << "\tWezly: ";
        for (int i = 0; i < 4; ++i) cout << e.ID[i] << (i < 3 ? ", " : "");
        cout << "\n";
    }

    printWrapped(bc_out_nodes, "\n*BC_OUT");
    printWrapped(bc_in_nodes, "\n*BC_IN");
    printWrapped(grid.boundaryNodes, "\n*BC");


    ElemUniv eUniv(globalData.npc);
    int npc = eUniv.npc;
    Surface surface(globalData.npc);





    vector<vector<double>> Hglobal(grid.nN, vector<double>(grid.nN, 0.0));
    vector<vector<double>> Cglobal(grid.nN, vector<double>(grid.nN, 0.0));

    vector<double> Pglobal(grid.nN, 0.0);
    vector<vector<double>> Hbcglobal(grid.nN, vector<double>(grid.nN, 0.0));
    vector<vector<double>> Hglobal_plus_Hbc(grid.nN, vector<double>(grid.nN, 0.0));



    for (auto& elem : grid.elements) {

        elem.jakobiany.clear();

        for (int p = 0; p < npc * npc; ++p) {
            Jakobian J(elem, grid, eUniv, p);

            if (fabs(J.detJ) < 1e-12) {
                cerr << "Warning: detJ bliski 0 w elemencie " << elem.id << " przy pc=" << p << "\n";
            }
            elem.jakobiany.push_back(J);
            // J.coutJakobian(); 
        }


        //cout << fixed << setprecision(6);
        //cout << "\n -------------------------------------------------------------------------------------------------------------------\n";
        //cout << "\n\n\nElement " << elem.id << "\n";
        //cout << "Node coords (N1..N4):\n";
        //for (int a = 0; a < 4; ++a) {
        //    int nid = elem.ID[a] - 1;
        //    cout << " N" << (a + 1) << " id=" << elem.ID[a]
        //        << " (" << grid.nodes[nid].x << ", " << grid.nodes[nid].y << ")";
        //    cout << "   BC: " << grid.nodes[nid].BC << endl;
        //}


        //for (int p = 0; p < npc * npc; ++p) {
        //    const Jakobian& J = elem.jakobiany[p];
        //    int i = p % npc;
        //    int j = p / npc;
        //    double ksi = eUniv.gauss.xi[i];
        //    double eta = eUniv.gauss.xi[j];

        //    cout << "\npc=" << (p + 1) << " ksi=" << ksi << " eta=" << eta << "\n";
        //    cout << "Jakobian: \n";
        //    std::cout << "[" << J.J[0][0] << ", " << J.J[0][1] << "]\n";
        //    std::cout << "[" << J.J[1][0] << ", " << J.J[1][1] << "]\n";
        //    cout << " detJ = " << J.detJ << "\n";
        //    cout << " dN_dx: ";
        //    for (int a = 0; a < 4; ++a) cout << setw(12) << J.dN_dx[a];
        //    cout << "\n dN_dy: ";
        //    for (int a = 0; a < 4; ++a) cout << setw(12) << J.dN_dy[a];
        //    cout << "\n";
        //}


        auto it = globalData.materialIdToIndex.find(elem.materialId);
        if (it == globalData.materialIdToIndex.end()) {
            cerr << "Brak materialu dla elementu " << elem.id << endl;
            exit(1);
        }
        const Material& mat = globalData.materials[it->second];
        elem.obliczH(eUniv, grid, mat);

        //// Wypisz macierz H elementu
        //cout << "\nMacierz H elementu " << elem.id << ":\n";
        //cout << setw(12) << " " << "N1" << setw(12) << "N2" << setw(12) << "N3" << setw(12) << "N4" << endl;
        //cout << string(55, '-') << endl;
        //for (int i = 0; i < 4; ++i) {
        //    cout << "N" << i + 1 << " ";
        //    for (int j = 0; j < 4; ++j)
        //        cout << setw(12) << elem.H[i][j];
        //    cout << endl;
        //}




        elem.obliczHbc(surface, grid, globalData.AlfaOut, globalData.AlfaIn);

        // Wypisz macierz Hbc elementu
        //cout << "\nMacierz Hbc elementu " << elem.id << ":\n";
        //cout << setw(12) << " " << "N1" << setw(12) << "N2" << setw(12) << "N3" << setw(12) << "N4" << endl;
        //cout << string(55, '-') << endl;
        //for (int i = 0; i < 4; ++i) {
        //    cout << "N" << i + 1 << " ";
        //    for (int j = 0; j < 4; ++j)
        //        cout << setw(12) << elem.Hbc[i][j];
        //    cout << endl;
        //}



        elem.obliczP(surface, grid, globalData.AlfaOut, globalData.AlfaIn, globalData.TotOut, globalData.TotIn);



        // Wypisz wektor P elementu
       // cout << "\nWektor P elementu " << elem.id << ":\n";
       // for (int i = 0; i < 4; ++i) {
           // cout << "N" << i + 1 << " : " << elem.P[i] << "\n";
       // }



        for (int a = 0; a < 4; ++a) {
            int gi = elem.ID[a] - 1;
            Pglobal[gi] += elem.P[a];

            for (int b = 0; b < 4; ++b) {
                int gj = elem.ID[b] - 1;
                Hglobal[gi][gj] += elem.H[a][b];
                Hbcglobal[gi][gj] += elem.Hbc[a][b];
                Cglobal[gi][gj] += elem.C[a][b];
                Hglobal_plus_Hbc[gi][gj] += elem.H[a][b] + elem.Hbc[a][b];
            }
        }

    }

    // cout << "\nGlobalny wektor P:\n";
    // for (int i = 0; i < grid.nN; ++i)

      //   cout << "Node " << (i + 1) << " : " << Pglobal[i] << "\n";

   //  printMatrix(Hglobal, "Globalna macierz H:");

   //  printMatrix(Cglobal, "Globalna macierz C");

   //  printMatrix(Hbcglobal, "HBCGLOBAL");

    // printMatrix(Hglobal_plus_Hbc, "Globalna macierz H + Hbc");


    vector<vector<int>> materialElements(globalData.materials.size());
    for (const auto& el : grid.elements) {
        auto it = globalData.materialIdToIndex.find(el.materialId);
        if (it == globalData.materialIdToIndex.end()) continue;
        materialElements[it->second].push_back(el.id); // save el ID
    }





    double dt = globalData.SimulationStepTime;
    double totalTime = globalData.SimulationTime;
    int steps = static_cast<int>(totalTime / dt);


    vector<double> T(grid.nN, globalData.InitialTemp);
    vector<double> Tnew(grid.nN, 0.0);

    cout << "\nTemperatura w czasie\n";

    double time = 0.0;

    for (int step = 1; step <= steps; step++) {
        time += dt;


        vector<vector<double>> H_eff = Hglobal_plus_Hbc;
        vector<double> P_eff = Pglobal;

        for (int i = 0; i < grid.nN; ++i)
            for (int j = 0; j < grid.nN; ++j)
                H_eff[i][j] += Cglobal[i][j] / dt;

        for (int i = 0; i < grid.nN; ++i) {
            double acc = 0.0;
            for (int j = 0; j < grid.nN; ++j)
                acc += (Cglobal[i][j] / dt) * T[j];
            P_eff[i] += acc;
        }


        Tnew = solveLinearSystem(H_eff, P_eff);


        T = Tnew;




        //cout << "\nCzas t = " << time << " s\n";
        //for (int i = 0; i < grid.nN; i++) {
        //    cout << "Node " << (i + 1) << " : " << T[i] << "\n";
        //}

         //50
        //if (fabs(fmod(time, step)) < 1e-9) {

        //    double Tmin = T[0];
        //    double Tmax = T[0];

        //    for (int i = 1; i < grid.nN; i++) {
        //        Tmin = min(Tmin, T[i]);
        //        Tmax = max(Tmax, T[i]);
        //    }

        //    cout << fixed << setprecision(3);
        //    cout << "\n Time[s] = " << time
        //        << " | Tmin = " << Tmin
        //        << " | Tmax = " << Tmax << endl;
        //}

        int printEvery = 1;

        if (step % printEvery == 0) {
            cout << fixed << setprecision(3);
            cout << "\n\n======================================================================================";
            cout << "\n Time[s] = " << time << "\n";
            cout << "======================================================================================\n\n";

            for (size_t m = 0; m < globalData.materials.size(); ++m) {
                const auto& mat = globalData.materials[m];
                const auto& elemIds = materialElements[m];

                cout << "----------------------------------------\n";
                cout << " Material: " << mat.name << "   (id=" << mat.id << ")\n";
                cout << "----------------------------------------\n";

                if (elemIds.empty()) {
                    cout << "(brak elementow)\n";
                    cout << "----------------------------------------\n\n";
                    continue;
                }

                double overallMin = numeric_limits<double>::infinity();
                double overallMax = -numeric_limits<double>::infinity();

                double bcMin = numeric_limits<double>::infinity();
                double bcMax = -numeric_limits<double>::infinity();
                int bcElemCount = 0;


                for (int eid : elemIds) {

                    const Element& el = grid.elements[eid - 1];


                    double Te = 0.0;
                    bool touchesBC = false;

                    for (int k = 0; k < 4; ++k) {
                        int nid = el.ID[k] - 1;
                        Te += T[nid];

                        int bc = grid.nodes[nid].BC;
                        if (mat.id == 1 && bc == 1) touchesBC = true;
                        if (mat.id == 3 && bc == 2) touchesBC = true;
                    }
                    Te /= 4.0;

                    if (isfinite(Te)) {
                        overallMin = min(overallMin, Te);
                        overallMax = max(overallMax, Te);

                        if (touchesBC) {
                            bcMin = min(bcMin, Te);
                            bcMax = max(bcMax, Te);
                            bcElemCount++;
                        }
                    }
                }

                cout << left << setw(22) << "ElemAvg Tmin:" << overallMin << "\n";
                cout << left << setw(22) << "ElemAvg Tmax:" << overallMax << "\n";

                if (mat.id == 1) {

                    if (bcElemCount > 0) {
                        // cout << left << setw(22) << "ElemAvg OUT Tmin:" << bcMin << "\n";
                       //  cout << left << setw(22) << "ElemAvg OUT Tmax:" << bcMax << "\n";
                    }
                    else {
                        cout << "Brak elementow materialu na BC_OUT\n";
                    }
                }
                else if (mat.id == 3) {

                    if (bcElemCount > 0) {
                        //  cout << left << setw(22) << "ElemAvg IN Tmin:" << bcMin << "\n";
                         // cout << left << setw(22) << "ElemAvg IN Tmax:" << bcMax << "\n";
                    }
                    else {
                        cout << "Brak elementow materialu na BC_IN\n";
                    }
                }

                cout << "----------------------------------------\n\n";
            }
        }



    }


    return 0;
}