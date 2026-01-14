#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cctype>
#include <iomanip>


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





class Element {
public:
    int id;
    int ID[4];
    vector<Jakobian> jakobiany;
    vector<vector<double>> H;
    vector<vector<double>> Hbc;
    vector<double> P;
    vector<vector<double>> C;


    void obliczH(const ElemUniv& eUniv, const Grid& grid, double conductivity, double density, double specificHeat) {


        
        H.assign(4, vector<double>(4, 0.0));
        C.assign(4, vector<double>(4, 0.0));

            int npc = eUniv.npc;


            for (int p = 0; p < npc * npc; p++)
            {
                Jakobian& J = jakobiany[p];



                int col = p % npc;   
                int row = p / npc;     
                
                double w_x = eUniv.gauss.w[col];
                double w_y = eUniv.gauss.w[row];

                const vector<double>& fKszt = eUniv.N[p];

                for (int i = 0; i < 4; i++)
                {
                    for (int j = 0; j < 4; j++)
                    {

                        H[i][j] += conductivity * (J.dN_dx[i] * J.dN_dx[j] + J.dN_dy[i] * J.dN_dy[j]) * J.detJ * w_x * w_y;
                        C[i][j] += specificHeat * density * (fKszt[i] * fKszt[j]) * J.detJ * w_x * w_y;

                    }
                }

            }

        

    }




    void obliczHbc( const Surface& surface, const Grid& grid, double conductivity);
    void obliczP(const Surface& surface, const Grid& grid, double alfa, double Tot);
};



class Grid {
public:
    int nN = 0;
    int nE = 0;
    vector<Node> nodes;
    vector<Element> elements;
    vector<int> boundaryNodes;
};


void Element::obliczHbc( const Surface& surface, const Grid& grid, double alfa) {

    Hbc.assign(4, vector<double>(4, 0.0));

    const int edge[4][2] = { {0,1}, {1,2}, {2,3}, {3,0} };

    int npc = surface.npc;

    GaussQuadrature G(npc);


    for (int e = 0; e < 4; e++)
    {
        int n1 = ID[edge[e][0]] - 1;
        int n2 = ID[edge[e][1]] - 1;


        if (grid.nodes[n1].BC == 0 || grid.nodes[n2].BC == 0)
        {
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
                    Hbc[i][j] += alfa * fKsz[i] * fKsz[j] * detJ * w;
        }


    }
}

void Element::obliczP(const Surface& surface, const Grid& grid, double alfa, double Tot) {

    P.assign(4, 0.0);
    const int edge[4][2] = { {0,1}, {1,2}, {2,3}, {3,0} };

    int npc = surface.npc;

    GaussQuadrature G(npc);


    for (int e = 0; e < 4; e++)
    {
        int n1 = ID[edge[e][0]] - 1;
        int n2 = ID[edge[e][1]] - 1;


        if (grid.nodes[n1].BC == 0 || grid.nodes[n2].BC == 0)
        {
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
                
                    P[i] += Tot * alfa * fKsz[i] * detJ * w;
        }



    }
}

class GlobalData {
public:
    double SimulationTime = 0;
    double SimulationStepTime = 0;
    double Conductivity = 0;
    double Alfa = 0;
    double Tot = 0;
    double InitialTemp = 0;
    double Density = 0;
    double SpecificHeat = 0;
    int nN = 0;
    int nE = 0;

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
    // znajdŸ pierwsz¹ cyfrê, znak minus, lub kropkê
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


    while (getline(file, line)) {
        string trimmedLine = trim(line);
        if (trimmedLine.empty()) continue;

        // GLOBAL DATA
        if (trimmedLine.find("SimulationTime") != string::npos) {
            double v; if (parseNumberAfterKey(trimmedLine, v)) globalData.SimulationTime = v;
            continue;
        }
        if (trimmedLine.find("SimulationStepTime") != string::npos) {
            double v; if (parseNumberAfterKey(trimmedLine, v)) globalData.SimulationStepTime = v;
            continue;
        }
        if (trimmedLine.find("Conductivity") != string::npos) {
            double v; if (parseNumberAfterKey(trimmedLine, v)) globalData.Conductivity = v;
            continue;
        }
        if (trimmedLine.find("Alfa") != string::npos) {
            double v; if (parseNumberAfterKey(trimmedLine, v)) globalData.Alfa = v;
            continue;
        }
        if (trimmedLine.find("Tot") != string::npos) {
            double v; if (parseNumberAfterKey(trimmedLine, v)) globalData.Tot = v;
            continue;
        }
        if (trimmedLine.find("InitialTemp") != string::npos) {
            double v; if (parseNumberAfterKey(trimmedLine, v)) globalData.InitialTemp = v;
            continue;
        }
        if (trimmedLine.find("Density") != string::npos) {
            double v; if (parseNumberAfterKey(trimmedLine, v)) globalData.Density = v;
            continue;
        }
        if (trimmedLine.find("SpecificHeat") != string::npos) {
            double v; if (parseNumberAfterKey(trimmedLine, v)) globalData.SpecificHeat = v;
            continue;
        }
        if (trimmedLine.find("Nodes number") != string::npos) {
            double v; if (parseNumberAfterKey(trimmedLine, v)) globalData.nN = static_cast<int>(v);
            continue;
        }
        if (trimmedLine.find("Elements number") != string::npos) {
            double v; if (parseNumberAfterKey(trimmedLine, v)) globalData.nE = static_cast<int>(v);
            continue;
        }

        // Sekcje
        if (trimmedLine.find("*Node") != string::npos) {
            readingNodes = true;
            readingElements = false;
            continue;
        }
        if (trimmedLine.find("*Element") != string::npos) {
            readingElements = true;
            readingNodes = false;
            continue;
        }
        if (trimmedLine.find("*BC") != string::npos) {

            readingNodes = false;
            readingElements = false;
            readingBC = true;
            continue;
        }

        // Wczytywanie wezlow (oczekujemy >=3 pól: id, x, y)
        if (readingNodes) {
            auto parts = splitAndTrim(trimmedLine, ',');
            if (parts.size() < 3) {
                // linia nie zawiera poprawnych danych - pomijamy
                cerr << "Warning: pomijam niepoprawna linie w sekcji *Node: \"" << trimmedLine << "\"\n";
                continue;
            }
            try {
                int id = stoi(parts[0]);
                double x = stod(parts[1]);
                double y = stod(parts[2]);
                grid.nodes.emplace_back(id, x, y, 0);
            }
            catch (const exception& e) {
                cerr << "Warning: blad parsowania wezla: \"" << trimmedLine << "\" (" << e.what() << ")\n";
                continue;
            }
        }

        // Wczytywanie elementów (oczekujemy >=5 pól: id, n1, n2, n3, n4)
        if (readingElements) {
            auto parts = splitAndTrim(trimmedLine, ',');
            if (parts.size() < 5) {
                cerr << "Warning: pomijam niepoprawna linie w sekcji *Element: \"" << trimmedLine << "\"\n";
                continue;
            }
            try {
                Element e;
                e.id = stoi(parts[0]);
                for (int i = 0; i < 4; ++i) e.ID[i] = stoi(parts[i + 1]);
                grid.elements.push_back(e);
            }
            catch (const exception& e) {
                cerr << "Warning: blad parsowania elementu: \"" << trimmedLine << "\" (" << e.what() << ")\n";
                continue;
            }
        }

        if (readingBC) {
            auto parts = splitAndTrim(trimmedLine, ',');
            for (const auto& p : parts) {
                try {
                    int nodeID = stoi(p);
                    grid.boundaryNodes.push_back(nodeID);

                    if (nodeID <= grid.nodes.size()) {

                        grid.nodes[nodeID - 1].BC = 1;
                    }
                }
                catch (...) {
                    cerr << "Warning: blad parsowania wezla BC: \"" << trimmedLine << "\"\n";
                }
            }
        }

    } // while getline

    // ustawienie liczników
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




int main() {
    GlobalData globalData;
    Grid grid;
    
    vector<string> Pliki = { "Test1_4_4.txt", "Test2_4_4MixGrid.txt",
                             "Test3_31_31_kwadrat.txt", "Test4_testowe.txt" };


    //zmien
    if (!loadFromFile(Pliki[0], globalData, grid)) return 1;

   
    cout << "Dane globalne:\n";
    cout << "SimulationTime: " << globalData.SimulationTime << "\n";
    cout << "SimulationStepTime: " << globalData.SimulationStepTime << "\n";
    cout << "Conductivity: " << globalData.Conductivity << "\n";
    cout << "Alfa: " << globalData.Alfa << "\n";
    cout << "Tot: " << globalData.Tot << "\n";
    cout << "InitialTemp: " << globalData.InitialTemp << "\n";
    cout << "Density: " << globalData.Density << "\n";
    cout << "SpecificHeat: " << globalData.SpecificHeat << "\n";
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

    cout << "\nWezly z warunkami brzegowymi (BC):\n";
    for (int n : grid.boundaryNodes) cout << n << " ";
    cout << "\n\n";


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

       

        elem.obliczH(eUniv, grid, globalData.Conductivity, globalData.Density, globalData.SpecificHeat);

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




        elem.obliczHbc(surface, grid, globalData.Alfa);

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

        

        elem.obliczP(surface, grid, globalData.Alfa , globalData.Tot);
        
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

        for (int i = 0; i < grid.nN; i++) {
            for (int j = 0; j < grid.nN; j++) {
                H_eff[i][j] += Cglobal[i][j] / dt;
                P_eff[i] += (Cglobal[i][j] / dt) * T[j];
            }
        }

       
        Tnew = solveLinearSystem(H_eff, P_eff);

        
        T = Tnew;

       
        //cout << "\nCzas t = " << time << " s\n";
        //for (int i = 0; i < grid.nN; i++) {
        //    cout << "Node " << (i + 1) << " : " << T[i] << "\n";
        //}

         //50
        if (fabs(fmod(time, step)) < 1e-9) {

            double Tmin = T[0];
            double Tmax = T[0];

            for (int i = 1; i < grid.nN; i++) {
                Tmin = min(Tmin, T[i]);
                Tmax = max(Tmax, T[i]);
            }

            cout << fixed << setprecision(3);
            cout << "\n Time[s] = " << time
                << " | Tmin = " << Tmin
                << " | Tmax = " << Tmax << endl;
        }



    }


    return 0;
}
