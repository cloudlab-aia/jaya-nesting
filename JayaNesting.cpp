/*----------------------------------------------------------------------------*/
/*  FICHERO:       simutorno.cu									          */
/*  AUTOR:         Antonio Jimeno											  */
/*													                          */
/*  RESUMEN												                      */
/*  ~~~~~~~												                      */
/* Ejercicio grupal para simulación del movimiento de una herramienta         */
/* tipo torno utilizando GPUs                                                 */
/*----------------------------------------------------------------------------*/
//#define PRINTDEBUG

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <Windows.h>
#include <iostream>
#include <fstream>
#include <iomanip>


// includes, project

#include "Colisionador.h"
#include "JayaNesting.h"



typedef LARGE_INTEGER timeStamp;

double getTime();


void drawFile2(dxfFilter* filter, void*** polygons, const int& P, DL_Dxf* dxf, Point *cvPoints, const int &numCVPoints, const double& convexHullArea, const std::string& nameFile, Area* area) {
	DL_Codes::version exportVersion = DL_Codes::AC1009;
	DL_WriterA* dw = dxf->out(nameFile.c_str(), exportVersion);

	if (polygons == NULL)
		return;

	if (dw == NULL) {
		std::cerr << "NO se pudo abrir el archivo" << nameFile << std::endl;
		exit(1);
	}

	// Header
	dxf->writeHeader(*dw);

	// Aqui podemos poner Variables (si las necesitamos)

	// Tipo de unidades para dibujar por bloque (usamos mm)
	// Group - Variable name
	dw->dxfString(9, "$INSUNITS");
	// Group code - Value
	dw->dxfInt(70, 4);

	// Dimension line color (usamos ByLayer)
	dw->dxfString(9, "$DIMCLRD");
	dw->dxfInt(70, 256);

	// Dimension extension line color (usamos ByLayer)
	dw->dxfString(9, "$DIMCLRE");
	dw->dxfInt(70, 256);

	// Sets units for all dimension types except Angular (usamos decimal)
	dw->dxfString(9, "$DIMLUNIT");
	dw->dxfInt(70, 2);

	// Dimension line lineweight (usamos ByLayer) NOTA: en el pdf ByLayer = -2, en la libreria pone -1 creo
	dw->dxfString(9, "$DIMLWD");
	dw->dxfInt(70, -1);

	// XY drawing limits lower-left corner (esquina inferior en 0,0)
	dw->dxfString(9, "$LIMMIN");
	dw->dxfReal(10, 0.0);
	dw->dxfReal(20, 0.0);

	// Default polyline width (0?)
	dw->dxfString(9, "$PLINEWID");
	dw->dxfInt(40, 0);

	dw->sectionEnd();

	// Tables
	dw->sectionTables();

	// ViewPorts
	dxf->writeVPort(*dw);

	// Linetypes (almacenamos todos los linetypes admitidos por dxflib)
	dw->tableLineTypes(25);

	dxf->writeLineType(*dw, DL_LineTypeData("BYBLOCK", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("BYLAYER", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("CONTINUOUS", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("ACAD_ISO02W100", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("ACAD_ISO03W100", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("ACAD_ISO04W100", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("ACAD_ISO05W100", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("BORDER", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("BORDER2", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("BORDERX2", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("CENTER", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("CENTER2", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("CENTERX2", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("DASHDOT", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("DASHDOT2", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("DASHDOTX2", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("DASHED", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("DASHED2", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("DASHEDX2", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("DIVIDE", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("DIVIDE2", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("DIVIDEX2", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("DOT", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("DOT2", 0));
	dxf->writeLineType(*dw, DL_LineTypeData("DOTX2", 0));

	dw->tableEnd();

	dw->tableLayers(P);

	// Dibujamos 2 capas
	// Capa 0 (predeterminada)
	dxf->writeLayer(
		*dw,
		DL_LayerData("0", 0), // En la pagina pone el flag = 0 (segiundo parametro), en el pdf pone 2
		DL_Attributes(
			std::string(""), // LO dejamos vacio
			DL_Codes::white, // Color defecto
			100,             // Ancho defecto
			"CONTINUOUS"     // Tipo de linea por defecto
		)
	);

	// Capa Shapes (para las piezas)
	dxf->writeLayer(
		*dw,
		DL_LayerData("Polygons", 0),
		DL_Attributes(
			std::string(""),
			DL_Codes::white,
			100,
			"CONTINUOUS"
		)
	);

	// Capa ConvexHull (contorno del area)
	dxf->writeLayer(
		*dw,
		DL_LayerData("ConvexHull", 0),
		DL_Attributes(
			std::string(""),
			DL_Codes::yellow,
			100,
			"CONTINUOUS"
		)
	);

	// Capa Area
	dxf->writeLayer(
		*dw,
		DL_LayerData("Area", 0),
		DL_Attributes(
			std::string(""),
			DL_Codes::red,
			100,
			"CONTINUOUS"
		)
	);

	dw->tableEnd();

	// Other tables (son necesarias)
	// Tabla de simbolos de estilos
	dxf->writeStyle(*dw);
	// Tabla de símbolos View
	dxf->writeView(*dw);
	// Tabla de símbolos UCS
	dxf->writeUcs(*dw);

	dw->tableAppid(1);
	dw->tableAppidEntry(0x12);
	dw->dxfString(2, "ACAD");
	dw->dxfInt(70, 0);

	dw->tableEnd();

	// Estilos de dimension (definen el aspecto de las dimensiones)
	// NOTA: Esta sección es necesaria en VER_R13. Tener en cuenta que este método actualmente solo escribe
	// una sección DIMSTYLE falsa para hacer que el archivo sea legible por Autocad
	dxf->writeDimStyle(
		*dw,
		1, // arrowSize
		1, // extensionLineExtension
		1, // extensionLineOffset
		1, // dimensionGap
		1  // dimensionTextSize
	);

	// Tabla de siimbolos BLOCK
	dxf->writeBlockRecord(*dw);

	// Bloque minimo requerido (1)
	dxf->writeBlockRecord(*dw, "myblock1");
	//dxf.writeBlockRecord(*dw, "myblock2");

	dw->tableEnd();

	dw->sectionEnd();

	// Blocks (define las entidades de cada bloque)
	dw->sectionBlocks();

	dxf->writeBlock(
		*dw,
		DL_BlockData("*Model_Space", 0, 0.0, 0.0, 0.0));
	dxf->writeEndBlock(*dw, "*Model_Space");

	dxf->writeBlock(*dw,
		DL_BlockData("*Paper_Space", 0, 0.0, 0.0, 0.0));
	dxf->writeEndBlock(*dw, "*Paper_Space");

	dxf->writeBlock(*dw,
		DL_BlockData("*Paper_Space0", 0, 0.0, 0.0, 0.0));
	dxf->writeEndBlock(*dw, "*Paper_Space0");

	dxf->writeBlock(*dw,
		DL_BlockData("myblock1", 0, 0.0, 0.0, 0.0));

	dxf->writeEndBlock(*dw, "myblock1");

	dw->sectionEnd();

	// Entities Section (define las entidades del dibujo)
	dw->sectionEntities();

	// Dibujamos las piezas en la capa Polygons
	// Num pol + centroide + lowerPoint + BB + Num convex poly * (Num pts + puntos + axes + proyecciones)
	for (size_t i = 0; i < P; ++i) {
		// Dibujamos aristas
		for (int j = 0; j < *(int*)(polygons[i][4]); ++j) {
			Point* points = (Point*)(polygons[i][5]);
			dxf->writeLine(
				*dw,
				DL_LineData(
					points[j].x,
					points[j].y,
					0.0,
					points[(j + 1) % *(int*)(polygons[i][4])].x,
					points[(j + 1) % *(int*)(polygons[i][4])].y,
					0.0),
				DL_Attributes(
					"Polygons",
					256,
					-1,
					"BYLAYER"
				)
			);

			double* BB = (double*)(polygons[i][3]);
#ifdef PRINTDEBUG
			std::cout << "??" << std::endl;
			std::cout << BB[0] << std::endl;
#endif
			dxf->writeLine(
				*dw,
				DL_LineData(
					BB[2],
					BB[3],
					0.0,
					BB[0],
					BB[3],
					0.0
				),
				DL_Attributes(
					"BB",
					155,
					-1,
					"BYLAYER"
				)
			);

			dxf->writeLine(
				*dw,
				DL_LineData(
					BB[0],
					BB[3],
					0.0,
					BB[0],
					BB[1],
					0.0
				),
				DL_Attributes(
					"BB",
					155,
					-1,
					"BYLAYER"
				)
			);

			dxf->writeLine(
				*dw,
				DL_LineData(
					BB[0],
					BB[1],
					0.0,
					BB[2],
					BB[1],
					0.0
				),
				DL_Attributes(
					"BB",
					155,
					-1,
					"BYLAYER"
				)
			);

			dxf->writeLine(
				*dw,
				DL_LineData(
					BB[2],
					BB[1],
					0.0,
					BB[2],
					BB[3],
					0.0
				),
				DL_Attributes(
					"BB",
					155,
					-1,
					"BYLAYER"
				)
			);
		}

		// BB
		//for (int j = 0; j < 4; j++) {
			
		//}
	}

	if (area) {
		// Dibujamos el Area
		for (size_t i = 0; i < area->getNumPoints(); ++i) {
			dxf->writeLine(
				*dw,
				DL_LineData(
					area->getPoint(i).x,
					area->getPoint(i).y,
					0.0,
					area->getPoint((i + 1) % area->getNumPoints()).x,
					area->getPoint((i + 1) % area->getNumPoints()).y,
					0.0),
				DL_Attributes(
					"Area",
					256,
					-1,
					"BYLAYER"
				)
			);
		}
	}

	if (convexHullArea >= 0) {
		// Area del convex Hull
		Point lowerPoint(-50.0, 20.0);

		dxf->writeText(
			*dw,
			DL_TextData(
				lowerPoint.x + 30.0,
				lowerPoint.y - 20.0,
				0.0,
				lowerPoint.x + 30.0,
				lowerPoint.y - 20.0,
				0,
				8.0,
				1.0,
				0,
				0,
				0,
				std::to_string(convexHullArea) + " u",
				"standard",
				0
			),
			DL_Attributes(
				"ConvexHull",
				256,
				-1,
				"BYLAYER"
			)
		);

		// Dibujamos el contorno del Convex Hull
		// Dibujamos aristas
		for (int i = 0; i < numCVPoints; ++i) {
			dxf->writeLine(
				*dw,
				DL_LineData(
					cvPoints[i].x,
					cvPoints[i].y,
					0.0,
					cvPoints[(i + 1) % numCVPoints].x,
					cvPoints[(i + 1) % numCVPoints].y,
					0.0),
				DL_Attributes(
					"ConvexHull",
					256,
					-1,
					"BYLAYER"
				)
			);
		}
	}

	dw->sectionEnd();

	// Objects
	dxf->writeObjects(*dw);
	dxf->writeObjectsEnd(*dw);

	dw->dxfEOF();
	dw->close();

	delete dw;
}

Projection project3(const Point& axis, Point* points, const int& numPoints) {
	double p = points[0].x * axis.x + points[0].y * axis.y;
	double min = p;
	double max = min;

	//std::cout<<"Inicial: "<<min<<" "<<max<<std::endl;

	for (/*const auto &c : components*/int i = 0; i < numPoints - 1; ++i) {
		//p = c->getVertex1().getX() * axis.getX() + c->getVertex1().getY() * axis.getY();
		p = points[i].x * axis.x + points[i].y * axis.y;

		if (p < min) {
			min = p;
		}
		else if (p > max) {
			max = p;
		}

		//p = c->getVertex2().getX() * axis.getX() + c->getVertex2().getY() * axis.getY();
		p = points[i + 1].x * axis.x + points[i + 1].y * axis.y;

		if (p < min) {
			min = p;
		}
		else if (p > max) {
			max = p;
		}
	}

	//std::cout<<"AXIS: "<<axis.getX()<<" - "<<axis.getY()<<std::endl;

	return Projection(min, max);
}

// Num pol + centroide + lowerPoint + BB + Num convex poly * (Num pts + puntos + axes + proyecciones)
void move3(void** polygon, const double& x, const double& y/*const Point& p*/) {
	Point p1, p2;
	double vX, vY;
	Point p(x, y);
	// Número polígonos 
	int numPols = *(int*)(polygon[0]);
	// Centroide
	Point* cent = (Point*)(polygon[1]);
	cent->x += p.x;
	cent->y += p.y;
	// Lower Point
	Point* lw = (Point*)(polygon[2]);
	lw->x += p.x;
	lw->y += p.y;
	// Bounding Box
	double* BB = (double*)(polygon[3]);
	BB[0] += x;
	BB[1] += y;
	BB[2] += x;
	BB[3] += y;
	// Resto
	int index = 4;
	int numPts = 0;
	for (int i = 0; i < numPols; ++i) {
		// Puntos
		numPts = *(int*)(polygon[index++]);
		Point* points = (Point*)(polygon[index++]);
		for (int j = 0; j < numPts; ++j) {
			points[j].x += p.x;
			points[j].y += p.y;
		}
		// Axes
		Point* axes = (Point*)(polygon[index++]);
		// Projecciones
		Projection* proj = (Projection*)(polygon[index++]);
		Projection aux;
		for (int j = 0; j < numPts; ++j) {
			// Recalculamos los ejes de proyección y proyecciones
			p1.x = points[j].x;
			p1.y = points[j].y;
			p2.x = points[(j + 1) % numPts].x;
			p2.y = points[(j + 1) % numPts].y;
			vX = p2.x - p1.x;
			vY = p2.y - p1.y;
			axes[j] = Point(-vY, vX);

			aux = project3(axes[j], points, numPts);
			proj[j].min = aux.min;
			proj[j].max = aux.max;
		}
	}
}

void updateBB(double *BB, double x, double y) {
	if (x > BB[0]) BB[0] = x;
	if (y > BB[1]) BB[1] = y;
	if (x < BB[2]) BB[2] = x;
	if (y < BB[3]) BB[3] = y;
}

void rotate3(void** polygon, const double& angle) {
	Point p, newP, p1, p2;
	double vX, vY;
	// Num polígonos
	int numPols = *(int*)(polygon[0]);
	// Centroide
	Point* cent = (Point*)(polygon[1]);
	// Lower Point
	Point* lw = (Point*)(polygon[2]);
	// Bounding Box
	double* BB = (double*)(polygon[3]);
	int index = 4;
	int numPuntos = 0;
	//printf("------ >> ANGLEEE: %f", angle);
	for (int i = 0; i < numPols; ++i) {
		// Num puntos
		numPuntos = *(int*)(polygon[index++]);
		// Puntos
		Point* points = (Point*)(polygon[index++]);
		// Axes
		Point* axes = (Point*)(polygon[index++]);
		// Proyecciones
		Projection* proj = (Projection*)(polygon[index++]);
		for (int j = 0; j < numPuntos; ++j) {
			p.x = points[j].x - cent->x;
			p.y = points[j].y - cent->y;

			newP.x = p.x * cos(angle) - p.y * sin(angle);
			newP.y = p.y * cos(angle) + p.x * sin(angle);

			points[j].x = newP.x + cent->x;
			points[j].y = newP.y + cent->y;

			if (i == 0) {
				if (j == 0) {
					lw->x = points[j].x;
					lw->y = points[j].y;
				}
				else {
/*#ifdef PRINTDEBUG
					printf("------------------- \n");
					printf("%f %f - %f %f\n", points[j].x, points[j].y, lw->x, lw->y);
#endif*/
					if (points[j].y < lw->y || (points[j].y == lw->y && points[j].x < lw->x)) {
						lw->x = points[j].x;
						lw->y = points[j].y;
					}
				}
			}
		}

		for (int j = 0; j < numPuntos; j++) {
			// Recalculamos los ejes de proyección
			p1.x = points[j].x;
			p1.y = points[j].y;
			p2.x = points[(j + 1) % numPuntos].x;
			p2.y = points[(j + 1) % numPuntos].y;
			vX = p2.x - p1.x;
			vY = p2.y - p1.y;
			axes[j] = Point(-vY, vX);

			// Recalculamos projecciones
			proj[j] = project3(axes[j], points, numPuntos);
		}
	}

	// Actualizamos Bounding Box
	BB[0] = -MAXDOUBLE;
	BB[1] = -MAXDOUBLE;
	BB[2] = MAXDOUBLE;
	BB[3] = MAXDOUBLE;

	for (int i = 0; i < *(int*)(polygon[4]); i++) {
		updateBB(BB, ((Point*)polygon[5])[i].x, ((Point*)polygon[5])[i].y);
	}
}


double convexHull3(void*** polygons, const int& P, dxfFilter *filter, Point *&points, int &numPoints) {
	// Num pol + centroide + lowerPoint + BB + Num convex poly * (Num pts + puntos + axes + proyecciones)
	/*int numPoints = 0;
	for (int i = 0; i < P; ++i) {
		numPoints += *(int*)(polygons[i][3]);
	}*/

	if (numPoints < 3)
		return -1.0;

	//Point* points = new Point[numPoints];
	int lw = 0;
	int index = 0;
	// Num pol + centroide + lowerPoint + Num convex poly * (Num pts + puntos + axes + proyecciones)
	if (P == 0)
		return 0;

	// Cogemos el punto mas bajo y a la izquierda del conjunto de putnos que forman todos los polígonos
	//Point *point = calculateLowerPoint2(polygons, P);
	Point* point = (Point*)(polygons[0][2]);
	for (int i = 0; i < P; ++i) {
		//std::cout << point->x << ", " << point->y << " --- " << ((Point*)(polygons[i][2]))->x << ", " << ((Point*)(polygons[i][2]))->y << std::endl;
		if (((Point*)(polygons[i][2]))->y < point->y || (((Point*)(polygons[i][2]))->y == point->y && ((Point*)(polygons[i][2]))->x < point->x))
			point = (Point*)(polygons[i][2]);
	}

	// std::cout << "Punto mas bajo: " << point->x << ", " << point->y << std::endl;

	// Metemos todos los puntos en el array points, y nos quedamos con el índice (lw) que le corresponde al punto más bajo y a la 
	// izquierda encontrado anteriormente (point)
	for (size_t p = 0; p < P; ++p) {
		for (int i = 0; i < *(int*)(polygons[p][4]); ++i) {
			if (((Point*)(polygons[p][5]))[i].x == point->x && ((Point*)(polygons[p][5]))[i].y == point->y)
				lw = index;

			points[index++] = ((Point*)(polygons[p][5]))[i];
#ifdef PRINTDEBUG
			printf(" %f, %f\n", ((Point*)(polygons[p][5]))[i].x, ((Point*)(polygons[p][5]))[i].y);
#endif
			filter->setPointConvexHullPolygon(((Point*)(polygons[p][5]))[i]);
		}
	}
#ifdef PRINTDEBUG
	std::cout << "Puntos de donde sacamos CONVEXHULL:\n";
	for (int i = 0; i < numPoints; i++) {
		std::cout << i<<": " << points[i].x << ", " << points[i].y << "\n";
	}
	std::cout << std::endl;
#endif
	Point* cv = new Point[numPoints];
	index = 0;
	cv[index++] = points[lw];
	int nextPoint;
	int actualPoint = lw;

	do {
		nextPoint = (actualPoint + 1) % numPoints;
#ifdef PRINTDEBUG		
		std::cout << actualPoint << " -- " << nextPoint << std::endl;
#endif

		for (int i = 0; i < numPoints; ++i) {
			//Point *a = &points[actualPoint];
			//Point *b = &points[nextPoint];
			//Point *c = &points[i];
			double z = (points[i].y - points[actualPoint].y) * (points[nextPoint].x - points[i].x) - (points[i].x - points[actualPoint].x) * (points[nextPoint].y - points[i].y);

			if (z > 0.0)
				nextPoint = i;
		}
		//printf("---> %f, %f\n", points[nextPoint].x, points[nextPoint].y);
#ifdef PRINTDEBUG
		std::cout << "---> " << index << " -- " << actualPoint <<" ---" << lw << " != " << nextPoint << std::endl;
#endif
		actualPoint = nextPoint;
#ifdef PRINTDEBUG
		if (index > numPoints) {
			std::cout << "MAL" << std::endl;
		}
#endif
		if (actualPoint != lw && index<numPoints)
			cv[index++] = points[nextPoint];

	} while ((actualPoint != lw) && (index<numPoints));
#ifdef PRINTDEBUG
	std::cout << index << std::endl;
	for (int i = 0; i < index; i++)
		std::cout << cv[i].x << " * " << cv[i].y << std::endl;

	std::cout << index - 1 << " " << numPoints << std::endl;

	printf("------\n");
	for (int i = 0; i < index; ++i) {
		printf("Point: %f, %f\n", cv[i].x, cv[i].y);
	}
#endif
	// Calculamos area
	//int j = numPoints - 1;
	int j = index - 1;
	double area = .0;
	for (int i = 0; i < index; ++i) {
		area += (cv[j].x + cv[i].x) * (cv[j].y - cv[i].y);
		j = i;
	}

	// numPoints = número de puntos del convex hull final
	numPoints = index;
	// points = array puntos del convex hull final
	delete[] points;
	points = NULL;
	points = cv;

	return std::abs(area / 2);
}

void printPolygons(void*** polygons, const int *P) {
	printf("Número polígonos:%d \n", *P);
	// Num pol + centroide + lowerPoint + BB + Num convex poly * (Num pts + puntos + axes + proyecciones)
	for (int j = 0; j < *P; ++j) {
		int* nP = (int*)(polygons[j][0]);
		printf("CENTROIDE: %f %f\n", ((Point*)(polygons[j][1]))->x, ((Point*)(polygons[j][1]))->y);
		printf("LOWER POINT: %f %f\n", ((Point*)(polygons[j][2]))->x, ((Point*)(polygons[j][2]))->y);
		printf("Bounding Box:\n");
		for (int i = 0; i < 4; i++) {
			printf("%f\n", ((double*)(polygons[j][3]))[i]);
		}
		for (int i = 0; i < *nP; ++i) {
			printf("POLIGONO %d\n", j);
			int* numP = (int*)(polygons[j][4 * i + 4]);

			printf("Num puntos: %d\n", *numP);
			Point* pts = (Point*)(polygons[j][4 * i + 5]);

			for (int k = 0; k < *numP; ++k) {
				printf("Punto %d: %f, %f\n", k, pts[k].x, pts[k].y);
			}
			//int *numAxes = (int*)(buffer[0][j]);
			//j++;
			Point* axes = (Point*)(polygons[j][4 * i + 6]);

			for (int k = 0; k < *numP; ++k) {
				printf("Axes %d: %f, %f\n", k, axes[k].x, axes[k].y);
			}
			Projection* proj = (Projection*)(polygons[j][4 * i + 7]);

			for (int k = 0; k < *numP; ++k) {
				printf("Projections %d: %f, %f\n", k, proj[k].min, proj[k].max);
			}
		}

		printf("\n");
	}
}

const std::string nameOutPutCPUFile = "salidaCPU.dxf";

// Todo lo externo a vars, mejor moverlo fuera a mem global
double MyObjective(double* vars)
{
	double f = -1.0; Evaluations++;
	// Evaluación del nesting para las variables involucradas
	int numPoints = 0;
	Point *cvPoints = NULL;
	int collision;

	if (Polygons == NULL) return DBL_MAX;

	// Hacemos copia de los polígono originales, para trabajr sobre esta copia
	Filter->transferData(CPPolygons);

	int nnn = Filter->getNumPolygons();
	int* nP = NULL;
	nP = &nnn;
#ifdef PRINTDEBUG
	std::cout << "POLÍGONOS LEÍDOS\n\n";
	printPolygons(CPPolygons, nP);
	std::cout << "Movimientos:\n ";
	std::cout << NumPoly << std::endl;
#endif
	for (int i = 0; i < NumPoly; i++) {
#ifdef PRINTDEBUG
		std::cout << vars[i*3] << " " << vars[i*3 + 1] << " " << vars[i*3 + 2] << std::endl;
#endif
		move3(CPPolygons[i], vars[i * 3], vars[i * 3 + 1]);
		rotate3(CPPolygons[i], vars[i * 3 + 2]);
	}
	//move3(CPPolygons[0], 62.6453, 37.4859);
	//rotate3(CPPolygons[0], 0.142281);
#ifdef PRINTDEBUG
	std::cout << "YA MOVIDOS\n";
	printPolygons(CPPolygons, nP);
	drawFile2(Filter, CPPolygons, Filter->getNumPolygons(), Dxf, cvPoints, numPoints, f, nameOutPutCPUFile, AArea);
#endif

	// Comprobamos si hay colisión entre alguno de los polígonos
	// Calculamos convex hull
	for (int i = 0; i < NumPoly; ++i) {
		numPoints += *(int*)(CPPolygons[i][4]);
	}
	//drawFile2(Filter, CPPolygons, Filter->getNumPolygons(), Dxf, cvPoints, numPoints, f, nameOutPutCPUFile, AArea);
	cvPoints = new Point[numPoints];
	f = convexHull3(CPPolygons, NumPoly, Filter, cvPoints, numPoints);
	collision = Filter->collisionDetection2(CPPolygons, NumPoly);
#ifdef PRINTDEBUG
	if (collision == 0 && vars[0] != 0.0 && vars[1] != 0.0) {
		std::cout << "PARA "<< collision << std::endl;
		//while(true);
	}
#endif
	if (collision) f = BestSol + (fabs(f - BestSol) + 0.01*BestSol)*(double)collision;

	//if (AArea != NULL && !AArea->isInsideArea(cvPoints, numPoints))
	//	f *= 2.0; //f = DBL_MAX;

	// Fichero de salida con ´la posición final de las piezas
#ifdef PRINTDEBUG
	if (Dxf && !collision)
		drawFile2(Filter, CPPolygons, Filter->getNumPolygons(), Dxf, cvPoints, numPoints, f, nameOutPutCPUFile, AArea);
#endif
	// Liberamos memoria de CPPolygons
	Filter->freeData(CPPolygons);

	if (cvPoints) delete[] cvPoints;

	return f;
}

double MyObjective2(double* vars, const int n)
{
	double f = -1.0; Evaluations++;
	// Evaluación del nesting para las variables involucradas
	int numPoints = 0;
	Point *cvPoints = NULL;
	int collision;

	if (Polygons == NULL) return DBL_MAX;

	// Hacemos copia de los polígono originales, para trabajr sobre esta copia
	Filter->transferData(CPPolygons);

	int nnn = Filter->getNumPolygons();
	int* nP = NULL;
	nP = &nnn;
	//for (int i = 0; i < n; i++) {
	move3(CPPolygons[n-1], vars[0], vars[1]);
	rotate3(CPPolygons[n-1], vars[2]);
	//}
	for (int i = 0; i < n; ++i) {
		numPoints += *(int*)(CPPolygons[i][4]);
	}
	//drawFile2(Filter, CPPolygons, Filter->getNumPolygons(), Dxf, cvPoints, numPoints, f, nameOutPutCPUFile, AArea);
	cvPoints = new Point[numPoints];
	f = convexHull3(CPPolygons, n, Filter, cvPoints, numPoints);
	collision = Filter->collisionDetection2(CPPolygons, n);
	if (collision) f = BestSol + (fabs(f - BestSol) + 0.01 * BestSol) * (double)collision;

	// Liberamos memoria de CPPolygons
	Filter->freeData(CPPolygons);

	if (cvPoints) delete[] cvPoints;

	return f;
}


/*----------------------------------------------------------------------------*/
/*  FUNCION A PARALELIZAR  (versión secuencial-CPU)
/* %% Jaya algorithm */
/*----------------------------------------------------------------------------*/

/*
Actualiza la población de individuos x en función de los individuos
de la población actual y de los mejores y peores individuos de la población
actual.
Guarda la evaluacion del mejor individuo en la variable global BestSol y
su índice en la variable imin.
Guarda el índice del peor individuo en la variable imax.
*/
void UpdatePopulation(double** x, int& imin, int& imax, int& newbest, const int iter)
{
	double* xnew = new double[NumPoly * 3];
	double minVal = x[imin][NumPoly * 3];
	double maxVal = x[imax][NumPoly * 3];
	bool updateMax = false;
	for (int i = 0; i < AdaptPop; i++)
	{
		// Se crea un nuevo individuo en funcion del algoritmo Jaya
		// xnew(i, j) = x(i, j) + r(1)*(Best(j) - abs(x(i, j))) - r(2)*(worst(j) - abs(x(i, j)));
		// Se evalua y actualiza si el individuo es mejor\



		
		
		/*for (int j = 0; j < NumPoly * 3; j += 3)
		{
			xnew[j] = x[i][j] + coef_rand()*(x[imin][j] - fabs(x[i][j])) - coef_rand()*(x[imax][j] - fabs(x[i][j]));
			xnew[j + 1] = x[i][j + 1] + coef_rand()*(x[imin][j + 1] - fabs(x[i][j + 1])) - coef_rand()*(x[imax][j + 1] - fabs(x[i][j + 1]));
			xnew[j + 2] = x[i][j + 2] + coef_rand()*(x[imin][j + 2] - fabs(x[i][j + 2])) - coef_rand()*(x[imax][j + 2] - fabs(x[i][j + 2]));
		}*/
		

		
		for (int j = 0; j < NumPoly * 3; j += 3)
		{
			xnew[j] = x[i][j] + coef_rand() * (x[imin][j] - x[i][j]) - coef_rand() * (x[imax][j] - x[i][j]);
			xnew[j + 1] = x[i][j + 1] + coef_rand() * (x[imin][j + 1] - x[i][j + 1]) - coef_rand() * (x[imax][j + 1] - x[i][j + 1]);
			xnew[j + 2] = x[i][j + 2] + coef_rand() * (x[imin][j + 2] - x[i][j + 2]) - coef_rand() * (x[imax][j + 2] - x[i][j + 2]);
		}
		

		/*
		for (int j = 0; j < NumPoly * 3; j += 3)
		{
			xnew[j] = x[i][j] + coef_rand() * (x[imin][j] - x[i][j]) - (coef_rand() * (x[imax][j] - x[i][j]) / fabs(x[imax][j] - x[i][j] / (x[imin][j] - x[i][j])));
			xnew[j + 1] = x[i][j + 1] + coef_rand() * (x[imin][j + 1] - x[i][j + 1]) - (coef_rand() * (x[imax][j + 1] - x[i][j + 1]) / fabs(x[imax][j + 1] - x[i][j + 1] / (x[imin][j + 1] - x[i][j + 1])));
			xnew[j + 2] = x[i][j + 2] + coef_rand() * (x[imin][j + 2] - x[i][j + 2]) - (coef_rand() * (x[imax][j + 2] - x[i][j + 2]) / fabs(x[imax][j + 2] - x[i][j + 2] / (x[imin][j + 2] - x[i][j + 2])));
		}
		*/

		double eval = MyObjective(xnew);
		// Print eval and wait for key
		//std::cout << "Eval: " << eval << std::endl;
		//std::cin.get()
		if (eval<x[i][NumPoly * 3])
		{
			for (int j = 0; j < NumPoly * 3; j++) 
				x[i][j] = xnew[j];
			x[i][NumPoly * 3] = eval;
			if (eval < minVal) { 
				// Print previous worst and best when new best is found
				//std::cout << "=====================  HIT  ====================" << std::endl;
				//std::cout << "Iteration: " << iter << std::endl;
				//std::cout << "Prev Worst: " << maxVal << std::endl;
				//std::cout << "Prev Best: " << minVal << std::endl;
				//std::cout << "================================================" << std::endl;
				minVal = eval;
				imin = i; 
				BestSol = minVal; 
				newbest++;
			}
			if (imax == i) 
				updateMax = true;
			/*if (i == 0)
				std::cout << "Updated 0" << std::endl;
			else
				std::cout << "Hit" << std::endl;*/
		}
	}
	if (updateMax) {
		maxVal = x[0][NumPoly * 3];
		imax = 0;
		for (int i = 1; i < AdaptPop; i++)
		{
			if (x[i][NumPoly * 3] > maxVal) 
			{ 
				maxVal = x[i][NumPoly * 3]; 
				imax = i; 
			}
		}
		/*if (imax == 0) {
			std::cout << "0 is worst" << std::endl;
			std::cout << "newbest" << newbest << std::endl;
		}*/
	}

	// Decrementamos el número de individuos cada 500 iteraciones hasta llegar a 10. Si el ultimo individuo es el mejor o el peor, se copia al anterior (siempre que el anterior no sea el mejor o el peor)
	// Si el último individuo es el mejor o el peor y el anterior es el mejor o el peor, se copia al anterior del anterior
	
	if (AdaptativePopulation && iter % 50 == 0 && AdaptPop > 16) {
		if (imin == AdaptPop - 1 && imin != 0) {
			if (imax == AdaptPop - 1) {
				if (imin != AdaptPop - 2) {
					for (int j = 0; j < NumPoly * 3; j++)
						x[AdaptPop - 2][j] = x[AdaptPop - 1][j];
					x[AdaptPop - 2][NumPoly * 3] = x[AdaptPop - 1][NumPoly * 3];
					imin = AdaptPop - 2;
					AdaptPop--;
				}
			}
			else {
				for (int j = 0; j < NumPoly * 3; j++)
					x[AdaptPop - 1][j] = x[AdaptPop - 2][j];
				x[AdaptPop - 1][NumPoly * 3] = x[AdaptPop - 2][NumPoly * 3];
				imin = AdaptPop - 1;
				AdaptPop--;
			}
		}
		else if (imax == AdaptPop - 1 && imax != 0) {
			if (imin == AdaptPop - 1) {
				if (imax != AdaptPop - 2) {
					for (int j = 0; j < NumPoly * 3; j++)
						x[AdaptPop - 2][j] = x[AdaptPop - 1][j];
					x[AdaptPop - 2][NumPoly * 3] = x[AdaptPop - 1][NumPoly * 3];
					imax = AdaptPop - 2;
					AdaptPop--;
				}
			}
			else {
				for (int j = 0; j < NumPoly * 3; j++)
					x[AdaptPop - 1][j] = x[AdaptPop - 2][j];
				x[AdaptPop - 1][NumPoly * 3] = x[AdaptPop - 2][NumPoly * 3];
				imax = AdaptPop - 1;
				AdaptPop--;
			}
		}
		else {
			for (int j = 0; j < NumPoly * 3; j++)
				x[AdaptPop - 1][j] = x[AdaptPop - 2][j];
			x[AdaptPop - 1][NumPoly * 3] = x[AdaptPop - 2][NumPoly * 3];
			AdaptPop--;
		}
		std::cout << "AdaptPop: " << AdaptPop << std::endl;
	}
	
	// Se actualizan los nuevos mínimos y máximos
	//minVal = x[imin][NumPoly * 3];
	//maxVal = x[imax][NumPoly * 3];
	/*
	for (int i = 0; i < Population; i++)
	{
		if (x[i][NumPoly * 3] > maxVal) { maxVal = x[i][NumPoly * 3]; imax = i; }
		if (x[i][NumPoly * 3] < minVal) { 
			minVal = x[i][NumPoly * 3]; imin = i; 
			BestSol = minVal;
		}
	}
	*/
}

/*
 Ejecución de la simulación en CPU. Se ejecutan Runs simulaciones independientes
 de forma secuencial. Cada simulación se ejecuta durante Iterations iteraciones
 y se almacena la mejor solución obtenida en el array Solutions.
 La solución se almacena en un array de NumPoly*3+1 elementos, donde los tres
 primeros elementos de cada polígono son las coordenadas x,y de la posición
 y el tercer elemento es el ángulo de rotación. El último elemento es el valor
 de evaluacion (convex hull) de la solución.
 Devuelve ERRORSIM si hay algún error en la simulación o OKSIM si la simulación
 se ha ejecutado correctamente.
*/
int JayaCPU()
{
	double** x = NULL;
	int irun = 0;
	double lastRunBest = 0.0;
	while (irun < Runs)
	{
		std::cout << "RUN: " << irun << std::endl;
		int iter = 0;
		int imin, imax;
		int patience = 0;
		if (irun == 0) {
			x = CreatePopulation(imin, imax);
			InitialSol = BestSol;
		}
		else
			x = CreatePopulation(imin, imax, Solutions[irun-1]);
		
		int newbestsol = 0;
		fprintf(stdout, "Worst:%f\n", x[imax][NumPoly * 3]);
		fprintf(stdout, "Best:%f\n", x[imin][NumPoly * 3]);
		double lastWorst = x[imax][NumPoly * 3];
		AdaptPop = Population;
		while (iter < Iterations)
		{
			UpdatePopulation(x, imin, imax, newbestsol, iter);
			if ((lastWorst - x[imax][NumPoly * 3]) > (lastWorst * PatienceImprovementRateLimit)) {
				lastWorst = x[imax][NumPoly * 3];
				patience = 0;
			}
			else patience++;
			if (patience > 5000) {
				fprintf(stdout, "Early stopping\n");
				break;
			}
			// Imprimir el peor cada 100 iteraciones
			/*if (iter % 100 == 0) {
				fprintf(stdout, "Worst:%f\n", x[imax][NumPoly * 3]);
			}*/

			iter++;
		}
		TotalIterations += iter;
		fprintf(stdout, "Times newbest:%d\n", newbestsol);
		Hits += newbestsol;
		fprintf(stdout, "Worst:%f\n", x[imax][NumPoly * 3]);
		fprintf(stdout, "Worst index:%i\n", imax);
		fprintf(stdout, "Best:%f\n", x[imin][NumPoly * 3]);
		fprintf(stdout, "Best index:%i\n", imin);
		memcpy(Solutions[irun], x[imin], (NumPoly * 3 + 1)*sizeof(double)); // se copia la solucion
		irun++;
		if (irun == 1 || x[imin][NumPoly * 3] < lastRunBest)
			lastRunBest = x[imin][NumPoly * 3];
		else break;
		
	}
	RunsDone = irun;
	if (x != NULL) DeletePopulation(x);

	return OKSIM;	// Simulación CORRECTA
}


// ---------------------------------------------------------------
// ---------------------------------------------------------------
// ---------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////
//PROGRAMA PRINCIPAL
////////////////////////////////////////////////////////////////////////////////
int
runTest(int runs, int iterations, char* fDxf, char* fArea)
{


	Dxf = new DL_Dxf();
	AArea = new Area(fArea, Dxf);
	Filter = new dxfFilter(fDxf);
	// Leemos polígonos del DXF y se almacenan en las variables globales
	// IMPORTANTE: La estructura que debe tener el DXF es la siguiente:
	// 1. Cada polígono debe estar en una capa distinta
	// 2. Cada polígono/capa tiene una serie de entidades tipo LINE (DL_LineData) que forman el polígono
	// Cualquier otra entidad que no sea LINE no se lee
	Filter->readPolygons(fDxf, Dxf);
	Filter->transferData(Polygons);
	NumPoly = Filter->getNumPolygons();

	// Evaluaciones de funcion
	Evaluations = 0;
	// Parametros Jaya
	Runs = runs;
	Iterations = iterations;
	Population = POPULATION; // Antes era definible desde fuera ahora no

	// Crea el array de soluciones para cada ejecución y añade una variable más donde se almacena el valor objetivo
	if (CreateSolutions(Runs, NumPoly * 3 + 1) == ERRORSIM)
	{
		fprintf(stderr, "No se puede ubicar la memoria de soluciones\n");
		return ERRORSIM;
	}
	// Creación buffer resultados e inicializaciones para versiones CPU y GPU
	ObjetivoCPU = (double*)malloc(Runs*sizeof(double));
	srand(time(NULL));
	cpu_start_time = getTime();
	if (JayaCPU() == ERRORSIM)
	{
		fprintf(stderr, "Simulacion CPU incorrecta\n");
		DeleteSolutions(Runs);
		if (ObjetivoCPU != NULL) free(ObjetivoCPU);
		if (ObjetivoGPU != NULL) free(ObjetivoGPU);
		return ERRORSIM;
	}
	cpu_end_time = getTime();
	// Generamos DXF del mejor CPU
	int bestindex = RunsDone-1;
	double bestsol = MAXDOUBLE;
	/*for (int i = 0; i < Runs; i++)
	{
		ObjetivoCPU[i] = Solutions[i][NumPoly * 3];
		if (ObjetivoCPU[i] < bestsol) {
			bestsol = ObjetivoCPU[i]; bestindex = i;
			std::cout << "Bestsol: " << bestsol << std::endl;
		}
	}*/
	// Hacemos copia de los polígono originales, para trabajr sobre esta copia
	Filter->transferData(CPPolygons);  int numPoints = 0;

	/*int auxNumPoints = numPoints;
	for (int i = 0; i < NumPoly; i++) {
		auxNumPoints += *(int*)(CPPolygons[i][4]);
	}
	Point* auxcvPoints = new Point[auxNumPoints];
	double initialHull = convexHull3(CPPolygons, NumPoly, Filter, auxcvPoints, auxNumPoints);*/

	for (int i = 0; i < NumPoly; i++) {
		move3(CPPolygons[i], Solutions[bestindex][i * 3], Solutions[bestindex][i * 3 + 1]);
		rotate3(CPPolygons[i], Solutions[bestindex][i * 3 + 2]);
		numPoints += *(int*)(CPPolygons[i][4]);
	}
	Point* cvPoints = new Point[numPoints];
	bestsol = convexHull3(CPPolygons, NumPoly, Filter, cvPoints, numPoints);
	//printf("bestsol: %f\n", bestsol);
	drawFile2(Filter, CPPolygons, Filter->getNumPolygons(), Dxf, cvPoints, numPoints, bestsol, nameOutPutCPUFile, AArea);
	// Liberamos memoria de CPPolygons
	Filter->freeData(CPPolygons);
	if (cvPoints) delete[] cvPoints;
	
	// Computar estadísticos
	//double cpuMean, cpuBest, cpuWorst, cpuStdev;
	//CalculateStatistics(ObjetivoCPU, cpuMean, cpuBest, cpuWorst, cpuStdev);
	// Impresión de resultados
	//fprintf(fout, "%d\t%d\t%d\t%.3lf\t%f\t%f\t%1.4le\t%1.4le\t%1.4le\t%1.4le\t%1.4le\t%1.4le\t%1.4le\t%1.4le\t%d\n", Runs, Iterations, Population, cpu_end_time - cpu_start_time, cpuMean, cpuStdev, cpuBest, cpuWorst, Evaluations);
	//fout << std::setprecision(3) << std::fixed;
	//fout << std::scientific;
	//fout << std::setw(20) << Runs << std::setw(20) << Iterations << std::setw(20) << Population << std::setw(20) << cpu_end_time - cpu_start_time << std::setw(20) << InitialSol << std::setw(20) << bestsol << std::setw(20) << Evaluations << std::endl;

	// Guardar InitialSol y bestsol en fichero aparte llamado resultados.txt en una linea cada uno
	// Sobreescribimos el fichero
	std::ofstream fout2;
	fout2.open("resultados.txt", std::ofstream::out | std::ofstream::trunc);
	fout2 << cpu_end_time - cpu_start_time << std::endl << InitialSol << std::endl << bestsol << std::endl << Hits << std::endl << TotalIterations << std::endl << RunsDone << std::endl << Filter->getNumPolygons() << std::endl;
	fout2.close();

	// Limpieza de buffers
	DeleteSolutions(Runs);
	if (ObjetivoCPU != NULL) free(ObjetivoCPU);
	Filter->freeData(Polygons);
	return OKSIM;
}



int main(int argc, char** argv)
{

	// Parametros Jaya
	if (argc != 5) { fprintf(stderr, "Uso: %s figuras.dxf area.dxf population iterations\n", argv[0]); exit(1); }
	// Modify population and iterations
	POPULATION = atoi(argv[3]);
	Iterations = atoi(argv[4]);
	//char filename[255];
	//sprintf(filename, "%s_Pop%d.txt", METHODNAME, POPULATION);
	//std::ofstream fout;
	//fout.open(filename);
	//fprintf(fout, "Runs\tIter\tPopul\tSpeed\tTimeG\tTimeC\tMeanG\tMeanC\tDevG\tDevC\tBestG\tBestC\tWorstG\tWorstC\tNumeval\n");
	if (true) {
		//fout << std::setprecision(3) << std::fixed;
		//fout << std::setw(20) << "Runs" << std::setw(20) << "Iter" << std::setw(20) << "Popul" << std::setw(20) << "TimeC" << std::setw(20) << "Initial" << std::setw(20) << "BestC" << std::setw(20) << "Numeval" << std::endl;
		for (int runs = 0; runs < Numrun; runs++)
			for (int iterations = 0; iterations < Numiter; iterations++)
			{
				fprintf(stderr, "Computando run:%d iter:%d popul:%d...", Runbat[runs], Iterations, POPULATION);
				if (runTest(Runbat[runs], Iterations, argv[1], argv[2]) != OKSIM) fprintf(stderr, "Error!\n"); else fprintf(stderr, "Ok!\n");
				//fout.flush();
			}
		//fout.close();
	}
	fprintf(stdout, "\n------------- Fin de la computacion-----------------\n");
	//getchar();
	return 0;
}
	
////////// Funciones auxiliares
double getTime()
{
	timeStamp start;
	timeStamp dwFreq;
	QueryPerformanceFrequency(&dwFreq);
	QueryPerformanceCounter(&start);
	return double(start.QuadPart) / double(dwFreq.QuadPart);
}