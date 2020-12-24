#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <initializer_list>

const int QUANTITY_Y = 6;
const int QUANTITY_X = 3;

double* GaussMethod(double** Matrix);
void AproksimFunction();
class Matrix
{
	double Y[QUANTITY_Y] = { 0.1f, 0.4f, 0.9f, 1.6f, 2.5f, 3.6f };
	double CalcElement(unsigned int n, unsigned int m);
	double CalcB(unsigned short n, double Y[QUANTITY_Y]);
public:
	double** mat;
	Matrix(unsigned short height = 3, unsigned short width = 3)
	{
		mat = new double* [height];
		for (unsigned short step(0); step < width; step++)
		{
			mat[step] = new double[width];
		}
		FillMatrix();
	}

	~Matrix()
	{
		delete[] mat;
	}
	void FillMatrix();
	void Show();
};

class BazF
{
	double f1(double);
	double f2(double);
	double f3(double);
public:
	double** basic_f;
	BazF(std::initializer_list<double> x_arr)
	{
		basic_f = new double* [QUANTITY_Y];

		for (unsigned short step(0); step < QUANTITY_Y; step++)
		{
			basic_f[step] = new double[QUANTITY_X];
		}
		for (unsigned short step_y(0); step_y < QUANTITY_Y; step_y++)
		{
			for (unsigned short step_x(0); step_x < QUANTITY_X; step_x++)
			{
				if (step_x == 0)
				{
					basic_f[step_y][step_x] = f1(*(x_arr.begin() + step_y));
				}
				else if(step_x == 1)
				{
					basic_f[step_y][step_x] = f2(*(x_arr.begin() + step_y));
				}
				else
				{
					basic_f[step_y][step_x] = f3(*(x_arr.begin() + step_y));
				}
			}
		}
	}
	~BazF()
	{
		delete[] basic_f;
	}
};

int main()
{
	setlocale(LC_ALL, "russian");
	AproksimFunction();
	system("pause");
}

double BazF::f1(double x)
{
	return 1;
}

double BazF::f2(double x)
{
	return x;
}

double BazF::f3(double x)
{
	return pow(M_E, -x);
}

double Matrix::CalcElement(unsigned int n, unsigned int m)
{
	BazF basic_f({1.1f,2.1f,3.1f,4.1f,5.1f,6.1f });
	double result = 0;
	for (unsigned short i(0); i < QUANTITY_Y; i++)
	{
		result += (basic_f.basic_f[i][n-1] * basic_f.basic_f[i][m - 1]);
	}
	return result;
}

double Matrix::CalcB(unsigned short n, double Y[QUANTITY_Y])
{
	BazF basic_f({ 1.1f,2.1f,3.1f,4.1f,5.1f,6.1f });
	double result = 0;
	for (unsigned short i(0); i < QUANTITY_Y; i++)
	{
		result += (Y[i] * basic_f.basic_f[i][n - 1]);
	}
	return result;
}

void Matrix::FillMatrix()
{
	for (unsigned short i(0); i < 3; i++)
	{
		for (unsigned short j(0); j < 3; j++)
		{
			mat[i][j] = CalcElement(i+1, j+1);
		}
		mat[i][3] = CalcB(i + 1, Y);
	}
}

void Matrix::Show()
{
	for (unsigned short i(0); i < 3; i++)
	{
		for (unsigned short j(0); j < 4; j++)
		{
			std:: cout << mat[i][j] << " ";
		}
		std::cout << "\n";
	}
}
double* GaussMethod(double** Matrix)
{
	double* result = new double[3];
	for (unsigned short j(1); j < 3; j++)
	{
		double subdiv = Matrix[j][0] / Matrix[0][0];
		for (unsigned short i(0); i < 4; i++)
		{
			Matrix[j][i] = Matrix[j][i] - (Matrix[0][i] * subdiv);
		}
	}
	for (unsigned short j(2); j < 3; j++)
	{
		double subdiv = Matrix[j][1] / Matrix[1][1];
		for (unsigned short i(1); i < 4; i++)
		{
			Matrix[j][i] = Matrix[j][i] - (Matrix[1][i] * subdiv);
		}
	}

	result[2] = Matrix[2][3] / Matrix[2][2];
	result[1] = (Matrix[1][3] - Matrix[1][2] * result[2]) / Matrix[1][1];
	result[0] = (Matrix[0][3] - Matrix[0][1] * result[1] - Matrix[0][2] * result[2]) / Matrix[0][0];
	return result;
}
void AproksimFunction()
{
	double Y[QUANTITY_Y] = { 0.1f, 0.4f, 0.9f, 1.6f, 2.5f, 3.6f };
	Matrix mat_object;
	BazF basic_f({ 1.1f,2.1f,3.1f,4.1f,5.1f,6.1f });
	//Вывод изначальной матрицы
	std::cout << "Изначальная матрица\n";
	for (unsigned short i(0); i < 3; i++)
	{
		for (unsigned short j(0); j < 4; j++)
		{
			std::cout << mat_object.mat[i][j] << " ";
		}
		std::cout << std::endl;
	}
	double* answers = GaussMethod(mat_object.mat);
	double f[QUANTITY_Y] = {};
	double sigma[QUANTITY_Y] = {};
	double J = 0;
	for (unsigned short i(0); i < QUANTITY_Y; i++)
	{
		f[i] = *answers * basic_f.basic_f[i][0] + *(answers + 1) * basic_f.basic_f[i][1] + *(answers + 2) * basic_f.basic_f[i][2];
		sigma[i] = abs(f[i] - Y[i]);
		J += pow(sigma[i], 2);
	}

	//Отсюда идет вывод данных
	std::cout << "\nМатрица после обработки\n";
	for (unsigned short i(0); i < 3; i++)
	{
		for (unsigned short j(0); j < 4; j++)
		{
			std::cout << mat_object.mat[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "\nX,Y и Z в результате выполнения метода Гаусса\n";
	for (unsigned short i(0); i < 3; i++)
	{
		std::cout << *(answers + i) << " ";
	}
	std::cout << "\n";
	std::cout << "\nЗначения апрокс. функций\n";
	for (unsigned short i(0); i < QUANTITY_Y; i++)
	{
		std::cout << f[i] << " ";
	}
	std::cout << "\n";
	std::cout << "\nЗначения сигмы\n";
	for (unsigned short i(0); i < QUANTITY_Y; i++)
	{
		std::cout << sigma[i] << " ";
	}
	std::cout << "\n";
	std::cout << "\nJ = " << J << "\n";
	delete[] answers;
}