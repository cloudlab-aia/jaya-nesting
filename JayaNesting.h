/*----------------------------------------------------------------------------*/
/*  FICHERO:       simutornoCPU.h									          */
/*  AUTOR:         Antonio Jimeno											  */
/*													                          */
/*  RESUMEN												                      */
/*  ~~~~~~~												                      */
/* Fichero de definiciones y estructuras                                      */
/*    						                                                  */
/*----------------------------------------------------------------------------*/

#ifndef _JAYA_H_
#define _JAYA_H_

/*============================================================================ */
/* Constantes											                       */
/*============================================================================ */
#define ERRORSIM 1
#define OKSIM    0
#define MAXDOUBLE 1e40

#define METHOD 1
#define METHODNAME "NESTING" 
//#define POPULATION 8

#include <random>

int POPULATION = 8;
int Numiter = 1;
int Iterbat[1] {8000};
int Numrun = 1;
int Runbat[6]  {8, 16, 32, 64, 128, 256};


/*============================================================================ */
/* Variables Globales										                   */
/*============================================================================ */
int Runs;
int Iterations;
int Population;
int AdaptPop;
double BestSol;
double InitialSol;
double **Solutions;
double *ObjetivoCPU;
double *ObjetivoGPU;
int Evaluations;
int Hits=0;
int RunsDone = 0;
int TotalIterations = 0;
bool AdaptativePopulation = false;
double PatienceImprovementRateLimit = 0.01;
double cpu_start_time, cpu_end_time;
// Propias del nesting
int NumPoly;
void*** Polygons = NULL;
void*** CPPolygons = NULL;
DL_Dxf* Dxf = NULL;
Area* AArea = NULL;
dxfFilter* Filter = NULL;

double MyObjective(double* vars);

/*============================================================================ */
/* CONSTANTES NESTING: Dimensiones del Tablero				                   */
/*============================================================================ */
#define MINVARVALUEX -100.0//-120.0//-250.0
#define MAXVARVALUEX 100.0//120.0//250.0
#define MINVARVALUEY -100.0//-130.0//-250.0
#define MAXVARVALUEY 100.0//130.0//250.0
#define PI_2 6.283185307179586476925286766559

/*============================================================================ */
/* Funciones Estadísticas		    						                   */
/*============================================================================ */
void CalculateStatistics(double* value, double& mean, double& minval, double& maxval, double& stddev)

{

	double sum = 0.0;
	double temp = 0.0;
	minval = MAXDOUBLE;
	maxval = -MAXDOUBLE;
	for (int i = 0; i < Runs; i++)
	{
		sum += value[i];
		if (value[i]>maxval) maxval = value[i];
		if (value[i]<minval) minval = value[i];
	}
	mean = sum / (double)Runs;
	for (int i = 0; i < Runs; i++)	temp += (value[i] - mean) * (value[i] - mean);
	temp /= Runs;
	stddev = (temp>1e-6) ? sqrt(temp) : 0.0;
}



/*============================================================================ */
/* Funciones Aleatorias										                   */
/*============================================================================ */
std::random_device rd;
std::mt19937 mt(rd());
std::uniform_real_distribution<double> dist(0, std::nextafter(1, DBL_MAX));
inline double var_randX()
{
	//return MINVARVALUEX + (MAXVARVALUEX - MINVARVALUEX) * (double)rand() / ((double)RAND_MAX);
	return dist(mt);
}

inline double var_randY()
{
	//return MINVARVALUEY + (MAXVARVALUEY - MINVARVALUEY) * (double)rand() / ((double)RAND_MAX);
	return dist(mt);
}

inline double var_randA()
{
	//return PI_2 * (double)rand() / ((double)RAND_MAX);
	return PI_2 * dist(mt);
}

inline double coef_rand()
{
	//return (double)rand() / (double)RAND_MAX;
	return dist(mt);
}

/*============================================================================ */
/* Funciones de tratamiento de memoria							 */
/*============================================================================ */
void DeleteSolutions(int runs)
{
	if (Solutions != NULL)
	{
		for (int i = 0; i < runs; i++)
		if (Solutions[i] != NULL) free(Solutions[i]);
		free(Solutions);
		Solutions = NULL;
	}
	return;
}

/*
* Crea la matriz de soluciones
*/
int CreateSolutions(int runs, int vars)
{
	Solutions = (double**)malloc(runs*sizeof(void*));
	if (Solutions == NULL) return ERRORSIM;
	for (int j = 0; j < runs; j++)
	{
		Solutions[j] = (double*)malloc(vars*(int)sizeof(double));
		if (Solutions[j] == NULL)
		{
			DeleteSolutions(runs);
			return ERRORSIM;
		}
	}
	return OKSIM;
}

void DeletePopulation(double** x)
{
	if (x != NULL)
	{
		for (int i = 0; i < Population; i++)
		if (x[i] != NULL) free(x[i]);
		free(x);
		x = NULL;
	}
	return;
}


/*
 Crea la poblacion inicial	con valores aleatorios.
 La poblacion es una matriz de doubles de tamaño Population donde cada fila es un individuo.
 Cada individuo tiene NumPoly*3 + 1 elementos. Es decir, cada individuo tiene NumPoly poligonos
 y cada poligono tiene 3 variables (x, y, a) y un elemento mas para la evaluacion.
 Devuelve la poblacion.
*/
double** CreatePopulation(int& imin, int& imax, double* sol = NULL)
{
	// Se crea la poblacion
	double **x = NULL;
	x = (double**)malloc(Population*sizeof(void*));
	if (x == NULL) return NULL;
	for (int j = 0; j < Population; j++)
	{
		x[j] = (double*)malloc((NumPoly*3 + 1)*(int)sizeof(double)); // se crea hueco para variable y evaluacion
		if (x[j] == NULL)
		{
			DeletePopulation(x);
			return NULL;
		}
	}
	// Se evalúa
	imin = imax = 0;
	BestSol = MAXDOUBLE;
	double minVal = MAXDOUBLE;
	double maxVal = -MAXDOUBLE;
	for (int i = 0; i < Population; i++)
	{
		// Se crea cada individuo
		for (int j = 0; j < NumPoly; j ++)
		{
			// x[i][j] y x[i][j + 1] hay que ajustarlo en el rango del área
			x[i][j * 3] = (i != 0) ? var_randX() : ((sol != NULL) ? sol[j * 3] : 0.0);
			x[i][j*3 + 1] = (i != 0) ? var_randY() : ((sol != NULL) ? sol[j * 3 + 1] : 0.0);
			x[i][j*3 + 2] = (i != 0) ? var_randA() : ((sol != NULL) ? sol[j * 3 + 2] : 0.0);
		}
		// Se evalua cada individuo
		double eval = MyObjective(x[i]);

		if (eval>maxVal) { maxVal = eval, imax = i; }
		if (eval<minVal) { minVal = eval, imin = i; }
		BestSol = minVal;
		x[i][NumPoly*3] = eval;
	}
	return x;
}




#endif // _JAYA_H_