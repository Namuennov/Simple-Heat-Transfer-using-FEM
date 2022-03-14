#include <iostream>
#include <fstream>
#include <string>

int INTEGRATION_POINTS = 3;
int NO_HEIGHT = 4;
int NO_WIDTH = 4;
double HEIGHT = 0.1;
double WIDTH = 0.1;
double K_T = 25.;
double ALPHA = 300.;
double T_AMBIENT = 1200.;
double DENSITY = 7800.;
double C_P = 700.;
double T0 = 100.;
double DELTA_T = 50.;
double SIMULATION_TIME = 500.;

//struktura węzła
struct node
{
    double x, y;    //współrzędne
    bool BC;        //flaga warunku brzegowego
};

//struktura elementu
struct element
{
    int ID[4];
    double H[4][4] = {0.};
    double Hbc[4][4] = {0.};
    double P[4] = {0.};
    double C[4][4] = {0.};

    void displayH()
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++) printf("|%f\t", H[i][j]);
            printf("|\n");
        }
    }

    void displayHbc()
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++) printf("|%f\t", Hbc[i][j]);
            printf("|\n");
        }
    }

    void displayC()
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++) printf("|%f\t", C[i][j]);
            printf("|\n");
        }
    }
};

//struktura siatki
struct GRID
{
    double H, B;        //Height, Bigness
    int nH, nB;         //no. Height, no. Bigness
    int nN, nE;         //no. Nodes, no. Elements
    node* nodes;        //tablica węzłów
    element* elements;  //tablica elementów

    GRID()
    {
        H = 1.;
        B = 1.;
        nH = 2;
        nB = 2;
        nN = 4;
        nE = 1;
        nodes = new node[4];
        elements = new element[1];
    }

    GRID(double h, double b, int nh, int nb)
    {
        //przypisanie wartości pól
        H = h; B = b; nH = nh; nB = nb;
        nN = nH * nB;
        nE = (nH - 1) * (nB - 1);
        //alokacja pamięci
        nodes = new node[nN];
        elements = new element[nE];

        //zmienne pomocnicze
        int help1 = 0;
        int help2 = 0;
        //przypisanie ID elementom siatki w kierunku przeciwnym do ruchu wskazówek zegara
        for (int i = 0; i < (nB - 1); i++)
        {
            help1 = nH * i;
            for (int j = 0; j < (nH - 1); j++)
            {
                help2 = ((nH - 1) * i) + j;
                elements[help2].ID[0] = j + help1;
                elements[help2].ID[1] = j + help1 + nH;
                elements[help2].ID[2] = j + help1 + nH + 1;
                elements[help2].ID[3] = j + help1 + 1;
            }
        }

        //zmienne pomocnicze
        double dx = static_cast<double>(B / (nB - 1));
        double dy = static_cast<double>(H / (nH - 1));
        //przypisanie współrzędnych węzłom siatki
        for (int i = 0; i < (nB - 1); i++)
        {
            for (int j = 0; j < (nH - 1); j++)
            {
                help2 = (nH * i) + j;
                nodes[help2].x = dx * i;
                nodes[help2].y = dy * j;
            }
            //przypisanie współrzędnych każdemu ostatniemu węzłowi w linii w bardziej precyzyjnej manierze
            help2 = (nH * i) + nH - 1;
            nodes[help2].x = dx * i;
            nodes[help2].y = H;  //dokładne przypisanie
        }
        //przypisanie współrzędnych węzłom z ostatniej linii w bardziej precyzyjnej manierze
        for (int i = 0; i < nH; i++)
        {
            help2 = (nH * (nB - 1)) + i;
            nodes[help2].x = B;  //dokładne przypisanie
            nodes[help2].y = dy * i;
        }
        nodes[(nH * nB) - 1].y = H;  //dokładne przypisanie


        //przypisanie wartości BC (flagi warunku brzegowego)
        for (int i = 0; i < nN; i++) nodes[i].BC = false;
        for (int i = 0; i < nH; i++)
        {
            nodes[i].BC = true;
            nodes[nN - i - 1].BC = true;
        }
        for (int i = 1; i < nH; i++)
        {
            nodes[i * nH].BC = true;
            nodes[i * nH - 1].BC = true;
        }
    }

    ~GRID()
    {
        //dealokacja pamięci
        delete[] nodes;
        delete[] elements;
    }
};

//,,tablica" współczynników potrzebnych przy całkowaniu Gaussa (węzły i współczynniki)
struct tabIntegral
{
    int n;
    double* nodes;  //tablica węzłów xk
    double* coeff;  //tablica współczynników Ak

    tabIntegral(int N = 1)    //trzeba implementować dodatkową ,,sekcję" konstruktora dla każdego n z powodu korzystania ze współczynników stabelaryzowanych
    {
        n = N;

        if (n == 1)
        {
            nodes = new double[2];
            nodes[0] = -(1. / sqrt(3.));
            nodes[1] = 1. / sqrt(3.);
            coeff = new double[2];
            coeff[0] = 1.;
            coeff[1] = 1.;
        }
        if (n == 2)
        {
            nodes = new double[3];
            nodes[0] = -sqrt(0.6);
            nodes[1] = 0.;
            nodes[2] = sqrt(0.6);
            coeff = new double[3];
            coeff[0] = 5. / 9.;
            coeff[1] = 8. / 9.;
            coeff[2] = 5. / 9.;
        }
        //else printf("Not implemented yet!");
    }
};

//pomocnicza struktura punktu 2D
struct point
{
    double x;
    double y;
};

//struktura uniwersalnego czterowęzłowego dwuwymiarowego elementu
struct Element4_2D
{
    int n;                      //zmienna przydatna w destruktorze
    tabIntegral tab_help;       //pomocnicza tablica węzłów i współczynników do całkowania Gaussa
    double** pc;                //współrzędne punktów całkowania w układzie lokalnym (ksi, eta)
    double** shapeFunKsi;       //wartości pochodnych funkcji kształtu po ksi w punktach całkowania układu lokalnego
    double** shapeFunEta;       //wartości pochodnych funkcji kształtu po eta w punktach całkowania układu lokalnego

    point** pcS;                //współrzędne punktów całkowanie PO POWIERZCHNI
    double*** shapeFunWalls;    //wartości funkcji kształtu dla całkowania po powierzchni

    double*** shapeFun;         //wartości funkcji kształtu w punktach całkowania układu lokalnego

    Element4_2D(int n)
    {
        this->n = n;
        tab_help = tabIntegral(n - 1);
        pc = new double* [n * n];
        for (int i = 0; i < (n * n); i++) pc[i] = new double[2];

        //przypisanie współrzędnych punktów całkowania w odpowiedniej kolejności (>^<^>^<^>)
        for (int i = 0; i < n; i++)
        {
            if (i % 2 == 0)
            {
                for (int j = 0; j < n; j++)
                {
                    pc[(i * n) + j][0] = tab_help.nodes[j];
                    pc[(i * n) + j][1] = tab_help.nodes[i];
                }
            }
            else
            {
                for (int j = 0; j < n; j++)
                {
                    pc[(i * n) + j][0] = tab_help.nodes[n - 1 - j];
                    pc[(i * n) + j][1] = tab_help.nodes[i];
                }
            }
        }

        //dla elementu 4-węzłowego są zawsze 4 funkcje kształtu
        shapeFunKsi = new double* [4];
        for (int i = 0; i < 4; i++) shapeFunKsi[i] = new double[n * n];
        shapeFunEta = new double* [4];
        for (int i = 0; i < 4; i++) shapeFunEta[i] = new double[n * n];

        //wyliczenie wartości pochodnych funkcji kształtu we wszystkich punktach układu lokalnego
        //wartości pochodnych po zmiennej ksi
        for (int i = 0; i < (n * n); i++) shapeFunKsi[0][i] = (-0.25 * (1. - pc[i][1]));
        for (int i = 0; i < (n * n); i++) shapeFunKsi[1][i] = (0.25 * (1. - pc[i][1]));
        for (int i = 0; i < (n * n); i++) shapeFunKsi[2][i] = (0.25 * (1. + pc[i][1]));
        for (int i = 0; i < (n * n); i++) shapeFunKsi[3][i] = (-0.25 * (1. + pc[i][1]));
        //wartości pochodnych po zmiennej eta
        for (int i = 0; i < (n * n); i++) shapeFunEta[0][i] = (-0.25 * (1. - pc[i][0]));
        for (int i = 0; i < (n * n); i++) shapeFunEta[1][i] = (-0.25 * (1. + pc[i][0]));
        for (int i = 0; i < (n * n); i++) shapeFunEta[2][i] = (0.25 * (1. + pc[i][0]));
        for (int i = 0; i < (n * n); i++) shapeFunEta[3][i] = (0.25 * (1. - pc[i][0]));



        //wyznaczenie współrzędnych punktów całkowania po powierzchni
        pcS = new point* [4];
        for (int i = 0; i < 4; i++) pcS[i] = new point[n];
        //NUMERACJA ŚCIAN: dolna->prawa->górna->lewa
        //NUMERACJA PUNKTÓW NA ŚCIANIE: zgodnie z kierunkiem przechodzenia po ścianach
        for (int i = 0; i < n; i++)
        {
            pcS[0][i].x = tab_help.nodes[i];
            pcS[0][i].y = -1.;
            pcS[1][i].x = 1.;
            pcS[1][i].y = tab_help.nodes[i];
            pcS[2][i].x = tab_help.nodes[n-i-1];
            pcS[2][i].y = 1.;
            pcS[3][i].x = -1.;
            pcS[3][i].y = tab_help.nodes[n - i - 1];
        }

        shapeFunWalls = new double** [4];
        for (int i = 0; i < 4; i++)
        {
            shapeFunWalls[i] = new double* [4];
            for (int j = 0; j < 4; j++) shapeFunWalls[i][j] = new double[n];
        }
        for (int i = 0; i < 4; i++)
        {
            for (int k = 0; k < n; k++)
            {
                shapeFunWalls[i][0][k] = 0.25 * (1. - pcS[i][k].x) * (1. - pcS[i][k].y);
                shapeFunWalls[i][1][k] = 0.25 * (1. + pcS[i][k].x) * (1. - pcS[i][k].y);
                shapeFunWalls[i][2][k] = 0.25 * (1. + pcS[i][k].x) * (1. + pcS[i][k].y);
                shapeFunWalls[i][3][k] = 0.25 * (1. - pcS[i][k].x) * (1. + pcS[i][k].y);
            }
        }



        //wyznaczenie wartości funkcji kształtu w każdym z punktów całkowania
        //dla każdego z punktów całkowania zapisanie w Elemencie macierzy {N}*{N}T
        shapeFun = new double** [n*n];
        for (int i = 0; i < (n * n); i++)
        {
            shapeFun[i] = new double* [4];
            for (int j = 0; j < 4; j++) shapeFun[i][j] = new double[4];
        }
        double* shapeFunHelp = new double[4];
        for (int i = 0; i < (n * n); i++)
        {
            shapeFunHelp[0] = 0.25 * (1. - pc[i][0]) * (1. - pc[i][1]);
            shapeFunHelp[1] = 0.25 * (1. + pc[i][0]) * (1. - pc[i][1]);
            shapeFunHelp[2] = 0.25 * (1. + pc[i][0]) * (1. + pc[i][1]);
            shapeFunHelp[3] = 0.25 * (1. - pc[i][0]) * (1. + pc[i][1]);
            for (int j = 0; j < 4; j++) for (int k = 0; k < 4; k++) shapeFun[i][j][k] = shapeFunHelp[j] * shapeFunHelp[k];
        }
        delete[] shapeFunHelp;
    }

    ~Element4_2D()
    {
        //dealokacja pamięci
        for (int i = 0; i < (n * n); i++) delete[] pc[i];
        delete[] pc;
        for (int i = 0; i < 4; i++) delete[] shapeFunKsi[i];
        delete[] shapeFunKsi;
        for (int i = 0; i < 4; i++) delete[] shapeFunEta[i];
        delete[] shapeFunEta;

        for (int i = 0; i < 4; i++) delete[] pcS[i];
        delete[] pcS;
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++) delete[] shapeFunWalls[i][j];
            delete[] shapeFunWalls[i];
        }
        delete[] shapeFunWalls;

        for (int i = 0; i < (n*n); i++)
        {
            for (int j = 0; j < 4; j++) delete[] shapeFun[i][j];
            delete[] shapeFun[i];
        }
        delete[] shapeFun;
    }
};

//struktura Macierzy Jakobiego
struct Jacobian2D
{
    double jac[2][2] = { 0. };
    double reverse[2][2] = { 0. };
    double detJ;

    void display()
    {
        printf("\n| %f | %f |\t\t| %f | %f |\n| %f | %f |\t\t| %f | %f |\n", jac[0][0], jac[0][1], reverse[0][0], reverse[0][1], jac[1][0], jac[1][1], reverse[1][0], reverse[1][1]);
    }
};

//wyznaczenie Macierzy Jakobiego dla pojedynczego punktu całkowania
Jacobian2D computeJacobian2D(GRID* grid, Element4_2D* element, int no_element, int no_integrationPoint)
{
    Jacobian2D result;
    //zmienna pomocnicza
    int help;
    //obliczenie wartości elementów Macierzy Jakobiego
    for (int i = 0; i < 4; i++)
    {
        help = grid->elements[no_element].ID[i];
        result.jac[0][0] += grid->nodes[help].x * element->shapeFunKsi[i][no_integrationPoint];
        result.jac[0][1] += grid->nodes[help].y * element->shapeFunKsi[i][no_integrationPoint];
        result.jac[1][0] += grid->nodes[help].x * element->shapeFunEta[i][no_integrationPoint];
        result.jac[1][1] += grid->nodes[help].y * element->shapeFunEta[i][no_integrationPoint];
    }

    //obliczenie Jakobianu
    result.detJ = (result.jac[0][0] * result.jac[1][1]) - (result.jac[0][1] * result.jac[1][0]);
    //wyznaczenie wartości elementów macierzy odwrotnej do Macierzy Jakobiego
    result.reverse[0][0] = result.jac[1][1];
    result.reverse[0][1] = -result.jac[0][1];
    result.reverse[1][0] = -result.jac[1][0];
    result.reverse[1][1] = result.jac[0][0];
    for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) result.reverse[i][j] /= result.detJ;

    //zwrócenie wyniku
    return result;
}

//wyznaczenie macierzy H i C dla czterowęzłowego elementu
void computeH_andC_4points(GRID* grid, Element4_2D* helpElement, int no_element)
{
    //wyzerowanie liczonych zmiennych
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            grid->elements[no_element].H[i][j] = 0.;
            grid->elements[no_element].C[i][j] = 0.;
        }
    }
    //alokacja pamięci i zmienne pomocnicze
    Jacobian2D jac;
    double** dN_dx = new double*[INTEGRATION_POINTS * INTEGRATION_POINTS];
    double** dN_dy = new double*[INTEGRATION_POINTS * INTEGRATION_POINTS];
    for (int i = 0; i < (INTEGRATION_POINTS * INTEGRATION_POINTS); i++)
    {
        dN_dx[i] = new double[4];
        dN_dy[i] = new double[4];
    }
    double* detJ = new double[INTEGRATION_POINTS * INTEGRATION_POINTS];

    //obliczenie pochodnych funkcji kształtu po x, y
    for (int i = 0; i < (INTEGRATION_POINTS*INTEGRATION_POINTS); i++)
    {
        jac = computeJacobian2D(grid, helpElement, no_element, i);
        detJ[i] = jac.detJ;
        for (int j = 0; j < 4; j++)
        {
            dN_dx[i][j] = (jac.reverse[0][0] * helpElement->shapeFunKsi[j][i]) + (jac.reverse[0][1] * helpElement->shapeFunEta[j][i]);
            dN_dy[i][j] = (jac.reverse[1][0] * helpElement->shapeFunKsi[j][i]) + (jac.reverse[1][1] * helpElement->shapeFunEta[j][i]);
        }
    }

    double*** H = new double** [INTEGRATION_POINTS * INTEGRATION_POINTS];
    for (int i = 0; i < (INTEGRATION_POINTS * INTEGRATION_POINTS); i++)
    {
        H[i] = new double* [4];
        for (int j = 0; j < 4; j++) H[i][j] = new double[4];
    }

    //obliczenie wyrazów macierzy H
    //zmienne do wzięcia odpowiednich wag przy całkowaniu
    int a = 0;
    int b = 0;
    int c = 1;
    for (int i = 0; i < (INTEGRATION_POINTS * INTEGRATION_POINTS); i++)
    {
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                H[i][j][k] = ((dN_dx[i][j] * dN_dx[i][k]) + (dN_dy[i][j] * dN_dy[i][k])) * helpElement->tab_help.coeff[a] * helpElement->tab_help.coeff[b] * K_T * detJ[i];
            }
        }
        //zmiana wartości zmiennych do wybierania wag
        a+=c;
        if (a == INTEGRATION_POINTS || a==-1)
        {
            a -= c;
            b++;
            c *= -1;
        }
    }

    //zsumowanie macierzy z poszczególnych punktów całkowania i przypisanie wartości macierzy H elementu
    double help;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            help = 0.;
            for (int k = 0; k < (INTEGRATION_POINTS * INTEGRATION_POINTS); k++) help += H[k][i][j];
            grid->elements[no_element].H[i][j] = help;
        }
    }



    //obliczenie wyrazów macierzy C (bez mnożenia przez globalne współczynniki)
    //zmienne do wzięcia odpowiednich wag przy całkowaniu
    a = 0;
    b = 0;
    c = 1;
    for (int i = 0; i < (INTEGRATION_POINTS * INTEGRATION_POINTS); i++)
    {
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                grid->elements[no_element].C[j][k] += (helpElement->shapeFun[i][j][k] * helpElement->tab_help.coeff[a] * helpElement->tab_help.coeff[b] * detJ[i]);
            }
        }
        //zmiana wartości zmiennych do wybierania wag
        a += c;
        if (a == INTEGRATION_POINTS || a == -1)
        {
            a -= c;
            b++;
            c *= -1;
        }
    }
    //przemnożenie wyrazów macierzy C przez globalne współczynniki
    for (int j = 0; j < 4; j++)
        for (int k = 0; k < 4; k++)
            grid->elements[no_element].C[j][k] *= (DENSITY * C_P);

    //dealokacja pamięci
    for (int i = 0; i < (INTEGRATION_POINTS * INTEGRATION_POINTS); i++)
    {
        delete[] dN_dx[i];
        delete[] dN_dy[i];
    }
    delete[] dN_dx;
    delete[] dN_dy;
    delete[] detJ;
    for (int i = 0; i < (INTEGRATION_POINTS * INTEGRATION_POINTS); i++)
    {
        for (int j = 0; j < 4; j++) delete[] H[i][j];
        delete[] H[i];
    }
    delete[] H;
}

//wyznaczenie macierzy Hbc i wektora P dla czterowęzłowego elementu
void computeHbc_andP_4points(GRID* grid, Element4_2D* helpElement, int no_element)
{
    //wyzerowanie liczonych zmiennych
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++) grid->elements[no_element].Hbc[i][j] = 0.;
        grid->elements[no_element].P[i] = 0.;
    }
    //alokacja pamięci i zmienne pomocnicze
    int* nodes = new int[4];
    for (int i = 0; i < 4; i++) nodes[i] = grid->elements[no_element].ID[i];
    double detJ = 0.;
    double** Hbc;
    Hbc = new double* [4];
    for (int i = 0; i < 4; i++) Hbc[i] = new double[4];
    //zaalokowanie pamięci na wektor P
    double* P = new double[4];
    //pętla iterująca po ścianach
    for (int i = 0; i < 4; i++)
    {
        //sprawdzenie flag (warunku brzegowego)
        if (grid->nodes[nodes[i%4]].BC && grid->nodes[nodes[(i + 1) % 4]].BC)
        {
            //,,wyzerowanie" macierzy Hbc dla pojedynczej ściany
            for (int j = 0; j < 4; j++) for (int k = 0; k < 4; k++) Hbc[j][k] = 0.;
            //,,wyzerowanie" wektora P dla pojedynczej ściany
            for (int j = 0; j < 4; j++) P[j] = 0.;
            //obliczenie Jakobianu
            detJ = (sqrt(((grid->nodes[nodes[(i + 1) % 4]].x - grid->nodes[nodes[i % 4]].x) * (grid->nodes[nodes[(i + 1) % 4]].x - grid->nodes[nodes[i % 4]].x)) + (((grid->nodes[nodes[(i + 1) % 4]].y - grid->nodes[nodes[i % 4]].y) * (grid->nodes[nodes[(i + 1) % 4]].y - grid->nodes[nodes[i % 4]].y)))))/2.;
            //wyliczenie macierzy Hbc dla pojedynczej ściany (bez pomnożenia przez alfa, które jest później)
            for (int j = 0; j < 4; j++)
                for (int k = 0; k < 4; k++)
                {
                    for (int l = 0; l < INTEGRATION_POINTS; l++)
                        Hbc[j][k] += (helpElement->shapeFunWalls[i][j][l] * helpElement->shapeFunWalls[i][k][l] * helpElement->tab_help.coeff[l]);
                    Hbc[j][k] *= detJ;
                }
                        //Hbc[j][k] = detJ * ((helpElement->shapeFunWalls[i][j][0] * helpElement->shapeFunWalls[i][k][0]) + (helpElement->shapeFunWalls[i][j][1] * helpElement->shapeFunWalls[i][k][1]));
            //dodanie macierzy Hbc dla pojedynczej ściany do macierzy Hbc elementu (razem z przemnożeniem przez alfa)
            for (int j = 0; j < 4; j++)
                for (int k = 0; k < 4; k++) grid->elements[no_element].Hbc[j][k] += (Hbc[j][k] * ALPHA);
            //wyliczenie wektora P dla pojedynczej ściany (bez pomnożenia przez Jakobian i stałe, które są później)
            for (int j = 0; j < 4; j++)
                for (int k = 0; k < INTEGRATION_POINTS; k++)
                    P[j] += (helpElement->shapeFunWalls[i][j][k] * helpElement->tab_help.coeff[k]);
            //P[j] = helpElement->shapeFunWalls[i][j][0] + helpElement->shapeFunWalls[i][j][1];
            //dodanie wektora P dla pojedynczej ściany do wektora P elementu (razem z przemnożeniem przez Jakobian i stałe)
            for (int j = 0; j < 4; j++) grid->elements[no_element].P[j] += P[j] * detJ * ALPHA * T_AMBIENT;

        }
    }

    //dealokacja pamięci
    delete[] nodes;
    for (int i = 0; i < 4; i++) delete[] Hbc[i];
    delete[] Hbc;
    //dealokacja pamięci na wektor P
    delete[] P;
}

//wyznaczenie globalnej macierzy H
double** computeGlobalH(GRID* grid)
{
    //alokacja pamięci
    double** result = new double* [grid->nN];
    for (int i = 0; i < grid->nN; i++)
    {
        result[i] = new double[grid->nN];
        for (int j = 0; j < grid->nN; j++) result[i][j] = 0.;
    }

    //zmienna pomocnicza
    element* help;
    //obliczenie globalnej macierzy H
    //pętla iterująca po elementach
    for (int i = 0; i < grid->nE; i++)
    {
        help = &grid->elements[i];
        //dwie pętle iterujące po macierzy H elementu
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                result[help->ID[j]][help->ID[k]] += help->H[j][k];
                result[help->ID[j]][help->ID[k]] += help->Hbc[j][k];
            }
        }
    }

    //zwrócenie wyniku
    return result;
}

//wyznaczenie globalnej macierzy C
double** computeGlobalC(GRID* grid)
{
    //alokacja pamięci
    double** result = new double* [grid->nN];
    for (int i = 0; i < grid->nN; i++)
    {
        result[i] = new double[grid->nN];
        for (int j = 0; j < grid->nN; j++) result[i][j] = 0.;
    }

    //zmienna pomocnicza
    element* help;
    //obliczenie globalnej macierzy C
    //pętla iterująca po elementach
    for (int i = 0; i < grid->nE; i++)
    {
        help = &grid->elements[i];
        //dwie pętle iterujące po macierzy C elementu
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                result[help->ID[j]][help->ID[k]] += help->C[j][k];
            }
        }
    }

    //zwrócenie wyniku
    return result;
}

//wyznaczenie globalnego wektora P
double* computeGlobalP(GRID* grid)
{
    //alokacja pamięci
    double* result = new double [grid->nN];
    for (int i = 0; i < grid->nN; i++) result[i] = 0.;

    //zmienna pomocnicza
    element* help;
    //obliczenie globalnego wektora P
    //pętla iterująca po elementach
    for (int i = 0; i < grid->nE; i++)
    {
        help = &grid->elements[i];
        //pętla iterująca po wektorze P elementu
        for (int j = 0; j < 4; j++)
        {
            result[help->ID[j]] += help->P[j];
        }
    }

    //zwrócenie wyniku
    return result;
}

//=======
//rozwiązanie układu równań Metodą Eliminacji Gaussa
double* GaussElimination(double** A, double* b, int n)
{
    //walidacja danych wejściowych
    if (A == nullptr)
    {
        printf("Wrong arguments! (matrix A)\n");
        return nullptr;
    }
    if (b == nullptr)
    {
        printf("Wrong arguments! (vector b)\n");
        return nullptr;
    }
    if (n <= 0)
    {
        printf("Wrong arguments! (size lesser than 0)\n");
        return nullptr;
    }
    if (A[n - 1][n - 1] == NULL)
    {
        printf("Wrong arguments! (size of matrix A)\n");
        return nullptr;
    }
    if (b[n - 1] == NULL)
    {
        printf("Wrong arguments! (size of vector b)\n");
        return nullptr;
    }
    if (A[0][0] == 0.)
    {
        printf("Wrong arguments! (first element of matrix A equals 0)\n");
        return nullptr;
    }

    //alokacja pamięci na tablicę pomocniczą temp[][] - jako że ,,użytkownik nie widzi",
    //co się dzieje w funkcji, to można wpisać macierz A i wektor b do jednej tablicy
    double** temp = new double* [n];
    for (int i = 0; i < n; i++) temp[i] = new double[n + 1];

    //alokacja tablicy na wynik
    double* result = new double[n];

    //przepisanie do utworzonej tablicy wartości z tablic wejściowych
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) temp[i][j] = A[i][j];
    for (int i = 0; i < n; i++) temp[i][n] = b[i];

    //wyliczenie wartości znajdujących się w macierzy trójkątnej górnej
    //pierwsza pętla (ze zmienną i) odpowiada ,,etapowi", czyli eliminacji zmiennej Xi z równań
    for (int i = 0; i < (n - 1); i++)
    {
        //liczba przez którą dzielimy - wspólna dla całego etapu
        double a = temp[i][i];
        //druga pętla (ze zmienną j) odpowiada wierszowi, który jest modyfikowany
        for (int j = (n - 1); j > i; j--)
        {
            //współczynnik potrzebny do liczenia wartości elementu - wspólny dla całego wiersza
            double b = (temp[j][i] / a);
            //trzecia pętla (ze zmienną k) odpowiada elementom wiersza modyfikowanego
            for (int k = (i + 1); k < (n + 1); k++) temp[j][k] -= (b * temp[i][k]);
        }
    }

    //obliczenie wartości niewiadomych, czyli elementów tabeli result[]
    //pierwsza pętla (ze zmienną i) odpowiada liczonej niewiadomej Xi
    for (int i = (n - 1); i >= 0; i--)
    {
        //przepisanie wartości ze ,,zmodyfikowanego wektora b"
        result[i] = temp[i][n];
        //druga pętla odpowiada za odejmowanie kolejnych wartości od wyrazu wolnego
        for (int j = (n - i - 1); j > 0; j--) result[i] -= (temp[i][j + i] * result[j + i]);
        //podzielenie przez liczbę w celu uzyskania ostatecznego wyniku
        result[i] /= temp[i][i];
    }

    //dealokacja pamięci
    for (int i = 0; i < n; i++) delete[] temp[i];
    delete[] temp;
    //zwrócenie wyniku
    return result;
}
//=======

//rozwiązanie Równania Fouriera dla niestacjonarnej wymiany ciepła
double* solveFourierEquation(double** globalH, double** globalC, double* globalP, double* temperature, int size)
{
    //alokacja pamięci
    double** matrixFourier = new double* [size];
    for (int i = 0; i < size; i++) matrixFourier[i] = new double[size];
    double* vectorFourier = new double[size];

    //obliczenie macierzy do Równania Fouriera [H] + [C]/DELTA_T
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            matrixFourier[i][j] = globalH[i][j] + (globalC[i][j] / DELTA_T);
    //obliczenie wektora wyrazów wolnych do Równania Fouriera {P} + [C]/DELTA_T * {T0}
    //zmienna pomocnicza
    double help;
    //pętla po ,,długości wektora wyrazów wolnych"
    for (int i = 0; i < size; i++)
    {
        help = 0.;
        //przemnożenie macierzy 
        for (int j = 0; j < size; j++)
            help += (globalC[i][j] * temperature[j]);
        //podzielenie przez współczynnik wspólny dla całego wiersza
        help /= DELTA_T;
        //obliczenie wyrazu wektora wyrazów wolnych
        vectorFourier[i] = help + globalP[i];
    }
    //rozwiązanie układu równań Metodą Eliminacji Gaussa
    double* result = GaussElimination(matrixFourier, vectorFourier, size);

    //dealokacja pamięci
    for (int i = 0; i < size; i++) delete[] matrixFourier[i];
    delete[] matrixFourier;
    delete[] vectorFourier;
    //zwrócenie wyniku
    return result;
}

//znajdowanie minimum i maksimum w wektorze
double* findMinMaxT(double* temperature, int size)
{
    double* result = new double[2];
    //minimum
    result[0] = temperature[0];
    //maksimum
    result[1] = temperature[0];
    for (int i = 0; i < size; i++)
    {
        if (temperature[i] < result[0]) result[0] = temperature[i];
        if (temperature[i] > result[1]) result[1] = temperature[i];
    }
    return result;
}

//struktura siatki wczytywanej z pliku
struct GRID_FromFile : public GRID
{
    GRID_FromFile(std::string filename)
    {
        //inicjalizacja potrzebnych zmiennych
        std::ifstream file(filename);
        std::string line;
        int counter = 0;
        std::size_t pos;
        std::size_t pos2;
        std::string sub;
        int noNode = 0;
        int noElement = 0;
        //załadowanie danych z pliku
        while (!file.eof())
        {
            std::getline(file, line);
            switch (counter)
            {
            //załadowanie ,,parametrów globalnych"
            case 0:
                pos = line.find(" ");
                pos++;
                sub = line.substr(pos);
                SIMULATION_TIME = std::stod(sub);
                break;
            case 1:
                pos = line.find(" ");
                pos++;
                sub = line.substr(pos);
                DELTA_T = std::stod(sub);
                break;
            case 2:
                pos = line.find(" ");
                pos++;
                sub = line.substr(pos);
                K_T = std::stod(sub);
                break;
            case 3:
                pos = line.find(" ");
                pos++;
                sub = line.substr(pos);
                ALPHA = std::stod(sub);
                break;
            case 4:
                pos = line.find(" ");
                pos++;
                sub = line.substr(pos);
                T_AMBIENT = std::stod(sub);
                break;
            case 5:
                pos = line.find(" ");
                pos++;
                sub = line.substr(pos);
                T0 = std::stod(sub);
                break;
            case 6:
                pos = line.find(" ");
                pos++;
                sub = line.substr(pos);
                DENSITY = std::stod(sub);
                break;
            case 7:
                pos = line.find(" ");
                pos++;
                sub = line.substr(pos);
                C_P = std::stod(sub);
                break;
            case 8:
                pos = line.find(" ");
                pos++;
                sub = line.substr(pos);
                delete[] nodes;
                nN = std::stoi(sub);
                nodes = new node[nN];
                break;
            case 9:
                pos = line.find(" ");
                pos++;
                sub = line.substr(pos);
                delete[] elements;
                nE = std::stoi(sub);
                elements = new element[nE];
                break;
            //załadowanie współrzędnych węzłów
            case 10:
                std::getline(file, line);
                while(line.find("Element") == std::string::npos)
                {
                    pos = line.find(",");
                    sub = line.substr(0, pos);
                    noNode = std::stoi(sub) - 1;
                    pos++;
                    pos2 = line.find(",", pos);
                    sub = line.substr(pos, pos2-pos);
                    nodes[noNode].x = std::stod(sub);
                    pos2++;
                    sub = line.substr(pos2);
                    nodes[noNode].y = std::stod(sub);
                    std::getline(file, line);
                }
                break;
            //załadowanie ID elementów
            case 11:
                while (line.find("BC") == std::string::npos)
                {
                    pos = line.find(",");
                    sub = line.substr(0, pos);
                    noElement = std::stoi(sub) - 1;
                    for (int i = 0; i < 4; i++)
                    {
                        pos++;
                        pos2 = line.find(",", pos);
                        sub = line.substr(pos, pos2 - pos);
                        elements[noElement].ID[i] = std::stoi(sub) - 1;
                        pos = pos2;
                    }
                    std::getline(file, line);
                }
                break;
            //ustawienie flag warunku brzegowego
            case 12:
                for (int i = 0; i < nN; i++) nodes[i].BC = 0;
                bool check = true;
                pos = line.find(",");
                sub = line.substr(0, pos);
                nodes[std::stoi(sub) - 1].BC = 1;
                while (check)
                {
                    pos++;
                    pos2 = line.find(",", pos);
                    if (pos2 == std::string::npos) check = false;
                    sub = line.substr(pos, pos2 - pos);
                    nodes[std::stoi(sub) - 1].BC = 1;
                    pos = pos2;
                }
                break;
            }
            counter++;
        }
        file.close();
    }

    /*~GRID_FromFile()
    {
        ;
    }*/
};

int main()
{
    //==========================TEST DLA SIATKI WPROWADZANEJ W KODZIE==========================
    //Element4_2D elementTest = Element4_2D(INTEGRATION_POINTS);
    //GRID siata(HEIGHT, WIDTH, NO_HEIGHT, NO_WIDTH);
    //double** globalH;
    //double** globalC;
    //double* globalP;
    //double* temperature0 = new double[siata.nN];
    //double* temperature1 = new double[siata.nN];
    //double* minAndMax;
    //for (int i = 0; i < siata.nN; i++) temperature0[i] = T0;
    //printf("=Time 0s=\n");
    //for (int j = 0; j < siata.nN; j++) printf("%lf\n", temperature0[j]);
    //printf("\n\n\n");

    //int steps = SIMULATION_TIME / DELTA_T;
    //for (int i = 1; i <= steps; i++)
    //{
    //    for (int j = 0; j < ((NO_HEIGHT - 1) * (NO_WIDTH - 1)); j++)
    //    {
    //        computeH_andC_4points(&siata, &elementTest, j);
    //        computeHbc_andP_4points(&siata, &elementTest, j);
    //    }
    //    globalH = computeGlobalH(&siata);
    //    globalC = computeGlobalC(&siata);
    //    globalP = computeGlobalP(&siata);
    //    temperature1 = solveFourierEquation(globalH, globalC, globalP, temperature0, siata.nN);
    //    printf("=Time %fs=\n", (i*DELTA_T));
    //    for (int j = 0; j < siata.nN; j++) printf("%lf\n", temperature1[j]);
    //    minAndMax = findMinMaxT(temperature1, siata.nN);
    //    printf("T_min = %4.14lf\tT_max = %1.14lf\n", minAndMax[0], minAndMax[1]);
    //    printf("\n\n\n");
    //    delete[] temperature0;
    //    temperature0 = temperature1;
    //    delete[] minAndMax;
    //    for (int j = 0; j < siata.nN; j++)
    //    {
    //        delete[] globalH[j];
    //        delete[] globalC[j];
    //    }
    //    delete[] globalH;
    //    delete[] globalC;
    //    delete[] globalP;
    //}

    //==========================TEST DLA SIATKI WCZYTYWANEJ Z PLIKU==========================
    //Element4_2D elementTest = Element4_2D(INTEGRATION_POINTS);
    //GRID_FromFile siataPlik = GRID_FromFile("TestGrid.txt");
    //double** globalH;
    //double** globalC;
    //double* globalP;
    //double* temperature1 = new double[siataPlik.nN];
    //double* temperature0 = new double[siataPlik.nN];
    //double* minAndMax;
    //for (int i = 0; i < siataPlik.nN; i++) temperature0[i] = T0;
    //printf("=Time 0s=\n");
    //for (int j = 0; j < siataPlik.nN; j++) printf("%lf\n", temperature0[j]);
    //printf("\n\n\n");

    //int steps = SIMULATION_TIME / DELTA_T;
    //for (int i = 1; i <= steps; i++)
    //{
    //    for (int j = 0; j < siataPlik.nE; j++)
    //    {
    //        computeH_andC_4points(&siataPlik, &elementTest, j);
    //        computeHbc_andP_4points(&siataPlik, &elementTest, j);
    //    }
    //    globalH = computeGlobalH(&siataPlik);
    //    globalC = computeGlobalC(&siataPlik);
    //    globalP = computeGlobalP(&siataPlik);
    //    temperature1 = solveFourierEquation(globalH, globalC, globalP, temperature0, siataPlik.nN);
    //    printf("=Time %fs=\n", (i*DELTA_T));
    //    //for (int j = 0; j < siataPlik.nN; j++) printf("%lf\n", temperature1[j]);
    //    minAndMax = findMinMaxT(temperature1, siataPlik.nN);
    //    printf("T_min = %4.14lf\tT_max = %1.14lf\n", minAndMax[0], minAndMax[1]);
    //    printf("\n\n\n");

    //    delete[] temperature0;
    //    temperature0 = temperature1;
    //    delete[] minAndMax;
    //    for (int j = 0; j < siataPlik.nN; j++)
    //    {
    //        delete[] globalH[j];
    //        delete[] globalC[j];
    //    }
    //    delete[] globalH;
    //    delete[] globalC;
    //    delete[] globalP;
    //}
}
