#include <iostream>

#include <iostream>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <algorithm>
#include <iomanip>

using namespace std;

int Nn; //Число узлов (размерность глобальной матрицы)
int Ne; //Число элементов
double eps = 1e-15;
double delta = 1e-15;
double v = 1;
vector<double> nodes; //Координаты узлов
vector<vector<int>> elems; //Элементы разбиения(узлы и номер материала)
vector <pair<int, double>> mat; //Материал(лямбда,гамма)
vector <pair<int, double>> bc1; //Первое краевое условие (глобальный номер узла, значение ф-ции)
vector<pair<int, int>> bc2; //Второе краевое условие (глобальный номер узла, значение)
vector<int> bc3; //Третье краевое условие (глобальный номер узла)

//Матрицу храним в ленточном формате
vector<double> gl; //нижняя диагональ
vector<double> gu; //верхняя диагональ
vector<double> di; //главная диагональ
vector<double> vec; //вектор правой части
vector<double> q; //текущее приближение
vector<double> qn; //следующее приближение

double w = 1.20; //Параметр релаксации

double Func(double x) //функция правой части
{
    //return 3;
    //return x - 2;
    //return x - 1;
    //return -x;
    return x-1.0/(2*sqrt(x));
    //return x + sin(x);
}

double Lambda(int n, double q)
{
    //return 4;
    //return q * 2;
    //return q + 2;
    //return q * q;
    return sqrt(q);
    //return cos(q);
}

double Tetta(int n, double q)
{
    return q * 2;
}

double Betta(double q)
{
    return q;
}

double U_betta(double q)//функция третьего краевого условия
{
    return q;
}
double dL_dq(double q)
{
    return 1;
}

void Input()
{
    int N;
    ifstream in("node.txt");
    in >> Nn;
    q.resize(Nn);
    qn.resize(Nn);
    nodes.resize(Nn);
    di.resize(Nn);
    vec.resize(Nn);
    gl.resize(Nn);
    gu.resize(Nn);
    for (int i = 0; i < Nn; i++)
        in >> nodes[i];
    in.close();
    in.open("elem.txt");
    in >> Ne;
    elems.resize(Ne);
    for (int i = 0; i < Ne; i++)
    {
        elems[i].resize(3);
        in >> elems[i][0] >> elems[i][1] >> elems[i][2];
    }
    in.close();
    in.open("mat.txt");
    in >> N;
    mat.resize(N);
    for (int i = 0; i < N; i++)
        in >> mat[i].first >> mat[i].second;
    in.close();
    in.open("bc1.txt");
    in >> N;
    bc1.resize(N);
    for (int i = 0; i < N; i++)
        in >> bc1[i].first >> bc1[i].second;
    in.close();
    in.open("bc2.txt");
    in >> N;
    bc2.resize(N);
    for (int i = 0; i < N; i++)
        in >> bc2[i].first >> bc2[i].second;
    in.close();
    in.open("bc3.txt");
    in >> N;
    bc3.resize(N);
    for (int i = 0; i < N; i++)
    {
        in >> bc3[i];
    }
    in.close();
}

void Matrix_Build()
{
    vector<vector<double>> G; //локальная матрица жесткости
    vector<vector<double>> M; //локальная матрица массы
    vector<double> b; //локальный вектор правой части
    b.resize(2);
    G.resize(2);
    for (int i = 0; i < 2; i++)
        G[i].resize(2);
    M.resize(2);
    for (int i = 0; i < 2; i++)
        M[i].resize(2);
    double h;//Шаг
    for (int i = 0; i < Ne; i++)//Заполняем локальные матрицы и вносим в глобальную
    {
        h = nodes[elems[i][1] - 1] - nodes[elems[i][0] - 1];
        G[0][0] = G[1][1] = (Lambda(mat[elems[i][2] - 1].first, q[elems[i][0] - 1]) + Lambda(mat[elems[i][2] - 1].first, q[elems[i][1] - 1])) / (2 * h);
        G[0][1] = G[1][0] = -(Lambda(mat[elems[i][2] - 1].first, q[elems[i][0] - 1]) + Lambda(mat[elems[i][2] - 1].first, q[elems[i][1] - 1])) / (2 * h);

        M[0][0] = M[1][1] = mat[elems[i][2] - 1].second * h / 3.;
        M[0][1] = M[1][0] = mat[elems[i][2] - 1].second * h / 6.;

        b[0] = h / 6. * (2 * Func(nodes[elems[i][0] - 1]) + Func(nodes[elems[i][1] - 1]));
        b[1] = h / 6. * (Func(nodes[elems[i][0] - 1]) + 2 * Func(nodes[elems[i][1] - 1]));

        for (int j = 0; j < 2; j++)//Заполняем глобальный вектор и главную диагональ глобальной матрицы
        {
            vec[elems[i][j] - 1] += b[j];
            di[elems[i][j] - 1] += G[j][j] + M[j][j];
        }
        gu[elems[i][0] - 1] += G[0][1] + M[0][1];
        gl[elems[i][0]] += G[1][0] + M[1][0];
    }
}

void Newton_Matrix_Build()
{
    vector<vector<double>> G; //локальная матрица жесткости
    vector<vector<double>> M; //локальная матрица массы
    vector<vector<double>> A; //локальная матрица производных
    vector<double> b; //локальный вектор правой части
    b.resize(2);
    G.resize(2);
    for (int i = 0; i < 2; i++)
        G[i].resize(2);
    M.resize(2);
    for (int i = 0; i < 2; i++)
        M[i].resize(2);
    A.resize(2);
    for (int i = 0; i < 2; i++)
        A[i].resize(2);
    double h;//Шаг
    for (int i = 0; i < Ne; i++)//Заполняем локальные матрицы и вносим в глобальную
    {
        h = nodes[elems[i][1] - 1] - nodes[elems[i][0] - 1];
        G[0][0] = G[1][1] = (Lambda(mat[elems[i][2] - 1].first, q[elems[i][0] - 1]) + Lambda(mat[elems[i][2] - 1].first, q[elems[i][1] - 1])) / (2 * h);
        G[0][1] = G[1][0] = -(Lambda(mat[elems[i][2] - 1].first, q[elems[i][0] - 1]) + Lambda(mat[elems[i][2] - 1].first, q[elems[i][1] - 1])) / (2 * h);

        M[0][0] = M[1][1] = mat[elems[i][2] - 1].second * h / 3.;
        M[0][1] = M[1][0] = mat[elems[i][2] - 1].second * h / 6.;

        double dLdQ1_left, dLdQ1_right, dLdQ2_left, dLdQ2_right;
        dLdQ1_left = -dL_dq(q[elems[i][0] - 1]) / h;
        dLdQ1_right = -dL_dq(q[elems[i][1] - 1]) / h;
        dLdQ2_left = dL_dq(q[elems[i][0] - 1]) / h;
        dLdQ2_right = dL_dq(q[elems[i][1] - 1]) / h;

        double dA11dQ1, dA11dQ2, dA12dQ1, dA12dQ2, dA21dQ1, dA21dQ2, dA22dQ1, dA22dQ2;

        dA11dQ1 = dA22dQ1 = (dLdQ1_left + dLdQ1_right) / (2.0 * h);
        dA11dQ2 = dA22dQ2 = (dLdQ2_left + dLdQ2_right) / (2.0 * h);
        dA12dQ1 = dA21dQ1 = -(dLdQ1_left + dLdQ1_right) / (2.0 * h);
        dA12dQ2 = dA21dQ2 = -(dLdQ2_left + dLdQ2_right) / (2.0 * h);


        A[0][0] = dA11dQ1 * q[elems[i][0] - 1] + dA12dQ1 * q[elems[i][1] - 1];
        A[0][1] = dA11dQ2 * q[elems[i][0] - 1] + dA12dQ2 * q[elems[i][1] - 1];
        A[1][0] = dA21dQ1 * q[elems[i][0] - 1] + dA22dQ1 * q[elems[i][1] - 1];
        A[1][1] = dA21dQ2 * q[elems[i][0] - 1] + dA22dQ2 * q[elems[i][1] - 1];


        b[0] = h / 6. * (2 * Func(nodes[elems[i][0] - 1]) + Func(nodes[elems[i][1] - 1])) + v * (q[elems[i][0] - 1] * (dA11dQ1 * q[elems[i][0] - 1] + dA11dQ2 * q[elems[i][1] - 1]) + q[elems[i][1] - 1] * (dA12dQ1 * q[elems[i][0] - 1] + dA12dQ2 * q[elems[i][1] - 1]));
        b[1] = h / 6. * (Func(nodes[elems[i][0] - 1]) + 2 * Func(nodes[elems[i][1] - 1])) + v * (q[elems[i][0] - 1] * (dA21dQ1 * q[elems[i][0] - 1] + dA21dQ2 * q[elems[i][1] - 1]) + q[elems[i][1] - 1] * (dA22dQ1 * q[elems[i][0] - 1] + dA22dQ2 * q[elems[i][1] - 1]));

        for (int j = 0; j < 2; j++)//Заполняем глобальный вектор и главную диагональ глобальной матрицы
        {
            vec[elems[i][j] - 1] += b[j];
            di[elems[i][j] - 1] += G[j][j] + M[j][j] + v * A[j][j];
        }
        gu[elems[i][0] - 1] += G[0][1] + M[0][1] + v * A[0][1];
        gl[elems[i][0]] += G[1][0] + M[1][0] + v * A[1][0];
    }
}

void Edge_Сonditions() //Учет краевых условий
{
    for (int i = 0; i < bc2.size(); i++)//Второе краевое
        vec[bc2[i].first - 1] += Tetta(bc2[i].second, q[bc2[i].first - 1]);
    for (int i = 0; i < bc3.size(); i++) //Третье краевое
    {
        vec[bc3[i] - 1] += Betta(q[bc3[i] - 1]) * U_betta(q[bc3[i] - 1]);
        di[bc3[i] - 1] += Betta(q[bc3[i] - 1]);
    }
    for (int i = 0; i < bc1.size(); i++)//Первое краевое
    {
        vec[bc1[i].first - 1] = bc1[i].second;
        di[bc1[i].first - 1] = 1;
        gu[bc1[i].first - 1] = 0;
        gl[bc1[i].first - 1] = 0;
    }
}

double Norm_Vector(vector<double>& vec) //Норма вектора
{
    double sum = 0;
    for (int i = 0; i < Nn; i++)
        sum += vec[i] * vec[i];
    return sqrt(sum);
}

void Matrix_Vector_Mult(vector<double>& res)
{
    for (int i = 0; i < Nn; i++)
    {
        if (i < 0)
            res[i] += gl[i] * q[i - 1];
        res[i] += di[i] * q[i];
        if (i < Nn - 1)
            res[i] += gu[i] * q[i + 1];
    }
}

void Factorization()
{
    for (int i = 1; i < Nn; i++)
    {
        di[i] = di[i] - gl[i] * gu[i - 1];
        gu[i] = gu[i] / di[i];
    }
}

void ForwardGauss(vector<double>& res)
{
    res[0] = vec[0] / di[0];
    for (int i = 1; i < Nn; i++)
        res[i] = (vec[i] - res[i - 1] * gl[i]) / di[i];
}

void BackwardGauss(vector<double>& res, vector<double>& f)
{
    res[Nn - 1] = f[Nn - 1];
    for (int i = Nn - 2; i >= 0; i--)
        res[i] = (f[i] - res[i + 1] * gu[i]);
}

void Iteration_Method()
{
    vector<double> y;
    vector<double> Ax;
    vector<double> qm;
    int k = 0;//Число итераций
    y.resize(Nn);
    Ax.resize(Nn);
    qm.resize(Nn);
    Matrix_Build();
    Edge_Сonditions();
    do
    {
        k++;
        Factorization();
        ForwardGauss(y);
        BackwardGauss(qn, y);
        for (int i = 0; i < Nn; i++)
        {
            qm[i] = qn[i] - q[i];
            q[i] = w * qn[i] + (1 - w) * q[i];
            di[i] = 0;
            vec[i] = 0;
            gl[i] = 0;
            gu[i] = 0;
            Ax[i] = 0;
        }
        Matrix_Build();
        Edge_Сonditions();
        Matrix_Vector_Mult(Ax);
        for (int i = 0; i < Nn; i++)
            Ax[i] -= vec[i];
        cout << scientific << setprecision(15) << Norm_Vector(Ax) / Norm_Vector(vec) << " " << Norm_Vector(qm) / Norm_Vector(q)<< endl;
    } while (Norm_Vector(Ax) / Norm_Vector(vec) >= eps && Norm_Vector(qm) / Norm_Vector(q) >= delta && k < 1000);
    if (Norm_Vector(Ax) / Norm_Vector(vec) >= eps)       cout << "Выход по точности\n";
    else if (Norm_Vector(qm) / Norm_Vector(q) >= delta)  cout << "Застой\n";
    else                                                 cout << "Превышен лимит операций\n";
    cout << "k = " << k << endl;
    for (int i = 0; i < Nn; i++)
        cout << scientific << setprecision(15) << q[i] << endl;
}

void Newton_Method()
{
    vector<double> y;
    vector<double> Ax;
    vector<double> qm;
    double N = 0;
    int k = 0;//Число итераций
    y.resize(Nn);
    Ax.resize(Nn);
    qm.resize(Nn);
    Newton_Matrix_Build();
    Edge_Сonditions();
    do
    {
        k++;
        Factorization();
        ForwardGauss(y);
        BackwardGauss(qn, y);
        for (int i = 0; i < Nn; i++)
        {
            qm[i] = qn[i] - q[i];
            q[i] = w * qn[i] + (1 - w) * q[i];
            di[i] = 0;
            vec[i] = 0;
            gl[i] = 0;
            gu[i] = 0;
            Ax[i] = 0;
        }
        Newton_Matrix_Build();
        Edge_Сonditions();
        Matrix_Vector_Mult(Ax);
        for (int i = 0; i < Nn; i++)
            Ax[i] -= vec[i];
        if (k % 5 == 1)
            cout << scientific << setprecision(15) << Norm_Vector(Ax) / Norm_Vector(vec) << endl;
        if (abs(N - Norm_Vector(Ax) / Norm_Vector(vec)) < 1e-15)
        {
            v = v / 2;
            for (int i = 0; i < Nn; i++)
            {
                di[i] = 0;
                vec[i] = 0;
                gl[i] = 0;
                gu[i] = 0;
                Ax[i] = 0;
            }
            Newton_Matrix_Build();
            Edge_Сonditions();
            Matrix_Vector_Mult(Ax);
            for (int i = 0; i < Nn; i++)
                Ax[i] -= vec[i];
        }
        N = Norm_Vector(Ax) / Norm_Vector(vec);
    } while (Norm_Vector(Ax) / Norm_Vector(vec) >= eps && Norm_Vector(qm) / Norm_Vector(q) >= delta && k < 1000);
    cout << "k = " << k << endl;
    for (int i = 0; i < Nn; i++)
        cout << scientific << setprecision(15) << q[i] << endl;
}

int main()
{
    setlocale(LC_ALL, "");
    Input();
    for (int i = 0; i < Nn; i++)
        q[i] = 1;
    cout << w << " Newton " << endl;
    Newton_Method();
    //cout << w << " Iteration_Method \n" << endl;
    //Iteration_Method();
}
