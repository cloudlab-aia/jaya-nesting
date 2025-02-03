#pragma once
#ifndef COLISIONADOR_H
#define COLISIONADOR_H

#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <stack>
#include "dl_dxf.h"
#include "dl_creationadapter.h"



// Array dinamico
template <class T>
class DynamicArray {
    friend std::ostream& operator<<(std::ostream& os, const DynamicArray& d) {
        for (auto e = d.arr; e != d.arr + d.capacity; ++e) {
            os << *e << " ";
        }

        return os;
    }

private:
    T* arr;
    // Capacidad actual
    size_t capacity;
    // Num elementos actuales en el array
    size_t numElements;

public:
    DynamicArray(size_t);
    DynamicArray();
    ~DynamicArray();

    T* getArr();
    T* get(int) const;
    void set(int, T);
    void pushBack(T);
    bool remove(int);
    void clear();
    T* back() const;
    T* front() const;
    size_t capty() const;
    size_t size() const;

    DynamicArray<T>& operator=(const DynamicArray<T>&);
};

// Clase DynamicArray
template <class T> 
DynamicArray<T>::DynamicArray(size_t N) {
    capacity = N;
    arr = new T[capacity];
    numElements = 0;
}

template <class T>
DynamicArray<T>::DynamicArray() {
    capacity = 5;
    arr = new T[capacity];
    numElements = 0;
}

template <class T>
DynamicArray<T>::~DynamicArray() {
    delete[] arr;
}

template <class T> 
T* DynamicArray<T>::get(int i) const {
    if (i < 0 || i >(int)capacity) {
        return NULL;
    }

    return &arr[i];
}

template <class T>
T* DynamicArray<T>::getArr() {
    return arr;
}

template <class T>
void DynamicArray<T>::set(int i, T value) {
    arr[i] = value;
}

template <class T>
void DynamicArray<T>::pushBack(T value) {
    if (numElements == capacity) {
        // No hay mas espacio
        // Doblamos el tamaño del array anterior
        T* new_array = new T[2 * capacity];

        for (int i = 0; i < (int)capacity; ++i) {
            new_array[i] = arr[i];
        }

        // Desaignamos la memoria de arr
        delete[] arr;

        arr = new_array;
        capacity = 2 * capacity;
    }

    set(numElements, value);
    ++numElements;
}

template <class T>
bool DynamicArray<T>::remove(int i) {
    if (i < 0 || i >(int)capacity) {
        return false;
    }

    int j = i;
    while (j < (int)numElements - 1) {
        arr[j] = arr[j + 1];

        ++j;
    }

    // ??
    arr[j].~Polygon();
    --numElements;

    return true;
}

template <class T>
void DynamicArray<T>::clear() {
    delete[] arr;
    arr = NULL;

    capacity = 5;
    numElements = 0;
    arr = new T[capacity];
}

template <class T>
T* DynamicArray<T>::back() const {
    return &arr[numElements - 1];
}

template <class T>
T* DynamicArray<T>::front() const {
    return &arr[0];
}

template <class T>
DynamicArray<T>& DynamicArray<T>::operator=(const DynamicArray<T>& d) {
    if (this == &d)
        return *this;

    if (numElements != d.numElements) {
        delete[] arr;
        arr = NULL;
        capacity = d.capacity;
        arr = new T[capacity];
        numElements = d.numElements;
    }

    for (int i = 0; i < (int)d.numElements; ++i)
        arr[i] = d.arr[i];

    return *this;
}

template <class T>
size_t DynamicArray<T>::capty() const {
    return capacity;
}

template <class T>
size_t DynamicArray<T>::size() const {
    return numElements;
}

// Puntos
struct Point {
    double x;
    double y;

    Point() {
        x = 0;
        y = 0;
    }

    Point(const double& a, const double& b) {
        x = a;
        y = b;
    }

    Point(const Point& p) {
        x = p.x;
        y = p.y;
    }

    Point(const Point &a, const Point &b) {
        x = (a.x + b.x) * 0.5;
        y = (a.y + b.y) * 0.5;
    }

    void operator=(const Point& p) {
        x = p.x;
        y = p.y;
    }

    Point operator+(const Point& p) const {
        Point res;
        res.x = x + p.x;
        res.y = y + p.y;
        return res;
    }

    /*Point operator-(const Point& p) const {
        return Point(x - p.x, y - p.y);
    }*/

    Point minus(const Point& p) const {
        return Point(x - p.x, y - p.y);
    }

    double operator*(const Point& p) const {
        return (x * p.x + y * p.y);
    }

    bool operator==(const Point& p) const {
        return ((x == p.x) && (y == p.y));
    }

    bool operator!=(const Point& p) const {
        return !((x == p.x) && (y == p.y));
    }

    Point operator/(const double& a) const {
        return Point(x / a, y / a);
    }

    bool operator<(const Point& p) const {
        return ((y < p.y) || (y == p.y && x < p.x));
    }

    int crossSign(const Point& b) const {
        return x * b.y - y * b.x;
    }
};

// Vector
class Vector {
private:
    //double x, y;

public:
    double x, y;

    Vector() { x = 0; y = 0; }
    Vector(const Vector& v) { this->x = v.x; this->y = v.y; }
    Vector(const double& x, const double& y) { this->x = x; this->y = y; }

    double getX() const { return x; }
    double getY() const { return y; }
    void setX(const double& x) { this->x = x; }
    void setY(const double& y) { this->y = y; }
};

// p - p = v
inline Vector operator-(const Point& a, const Point& b) {
    return Vector(a.x - b.x, a.y - b.y);
}

// v - v = v
inline Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a.getX() - b.getX(), a.getY() - b.getY());
}

// p + v = p
inline Point operator+(const Point* a, const Vector b) {
    return Point(a->x + b.getX(), a->y + b.getY());
}

// p + v = p
inline Point operator+(const Point a, const Vector b) {
    return Point(a.x + b.getX(), a.y + b.getY());
}

// v + v = v
inline Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a.getX() + b.getX(), a.getY() + b.getY());
}

// s * v = v
inline Vector operator*(const double& s, const Vector& v) {
    return Vector(s * v.getX(), s * v.getY());
}

// v * v = s (dot product)
inline double operator*(const Vector& a, const Vector& b) {
    return (a.getX() * b.getX() + a.getY() * b.getY());
}

inline Point mid(const Point &a, const Point &b) {
    return Point((a.x + b.x) * 0.5, (a.y + b.y) * 0.5);
}

inline Vector vabs(const Vector& a) {
    return Vector(fabs(a.getX()), fabs(a.getY()));
}

// Bezier
class Bezier {
    friend std::ostream& operator<<(std::ostream&, const Bezier&);

private:
    //double alignCurveX(const double&) const;
    //double alignCurveY(const double&) const;

public:
    //Point* p0, * p1, * p2, * p3;
    Point p0, p1, p2, p3;

    //Bezier(Point*, Point*, Point*, Point*);
    Bezier(Point, Point, Point, Point);
    Bezier();
    ~Bezier();

    void splitCurve(Bezier &a, Bezier &b) const;
    //bool intersectBB(const Bezier& a, const Bezier& b) const;

    bool intersect(Bezier b);
    //double** findIntersections(Bezier a, Bezier b);
    void parameterSplitLeft(double, Bezier&);
    //void recursivelyIntersect(Bezier, double, double, int, Bezier, double, double, int, double**, int&);

    bool findIntersection(Point, Point);

    Bezier& operator=(const Bezier&);
};

// Proyeccion
struct Projection {
    double min;
    double max;

    Projection() {
        min = 0;
        max = 0;
    }

    Projection(const double& a, const double& b) {
        min = a;
        max = b;
    }

    // Comprobamos si esta proyeccion solapa con la proyeccion pasada como parametro
    bool overlap(const Projection& p) const {
        if (p.min > max || p.max < min) {
            return false;
        }

        return true;
    }
};

// Area
class Area : public DL_CreationAdapter {
    friend std::ostream& operator<<(std::ostream&, const Area&);

private:
    std::ifstream file;
    std::vector<Point> points;

public:
    Area(const std::string&, DL_Dxf*);
    ~Area();
    virtual void addLine(const DL_LineData&);
    bool isInsideArea(const DynamicArray<Point>&) const;
    bool isInsideArea(const Point*, const int&) const;
    bool isConvex(const Point&, const Point&, const Point&) const;

    size_t getNumPoints() const { return points.size(); }
    const Point& getPoint(const int& i) const { return points[i]; }
};

// Poligono
class Polygon {
    friend std::ostream& operator<<(std::ostream&, const Polygon&);

private:
    std::string layer;
    // Aristas
    std::vector<Point> points;
    // Splines
    std::vector<Point> pointsBC;
    // Las subConvexPolygons tienen el centroide = (0, 0)
    Point centroid;
    Point lowerPoint;
    //int numPoints;
    int numPolygons;
    // x_min y_min x_max y_max
    std::vector<double> BB;
    
    std::vector<Polygon> subConvexPolygons;

    bool findIntersection(const int&, const Point&, const Point&) const;

public:
    std::vector<Point> axes;
    std::vector<Projection> projections;
    Polygon();
    Polygon(const int&, const std::string&);
    Polygon(const std::string&);
    ~Polygon();
    void clear();
    void initBB();
    Polygon(const Polygon&);

    Polygon& operator=(const Polygon&);

    Point& operator[](int i) {
        return points[i];
    }

    const Point& operator[](int i) const {
        return points[i];
    }

    const Point& getPoint(int i) const {
        return points[i];
    }

    const Point & getBCPoint(int i) const {
        return pointsBC[i];
    }

    void setPoint(Point p) {
        points.push_back(p);
    }

    int getNumPointsBC() const {
        return (int)pointsBC.size();
    }

    const std::vector<Point>& getPointsBC() const {
        return pointsBC;
    }

    int getNumPoints() const {
        return (int)points.size();
    }

    Point& getPoint(int i) {
        return points[i];
    }

    const Point& getCentroid() const {
        return centroid;
    }

    Point* getCentroid2() {
        return &centroid;
    }

    std::string getLayer() const {
        return layer;
    }

    const Point& getLowerPoint() const {
        return lowerPoint;
    }

    Point* getLowerPoint2() {
        return &lowerPoint;
    }

    void insertVertice(const Point& p) {
        points.push_back(p);
    }

    void insertVerticeBC(const Point& p) {
        pointsBC.push_back(p);
    }

    Projection getProjection(const int& i) const {
        return projections[i];
    }

    void setProjection(const Projection& p) {
        projections.push_back(p);
    }

    void setProjection(const Projection& p, const int& i) {
        projections[i] = p;
    }

    void setAxis(const Point& a) {
        axes.push_back(a);
    }

    int getNumAxis() const {
        return (int)axes.size();
    }

    Point getAxis(const int& i) const {
        return axes[i];
    }

    Projection project(const Point&) const;

    Point getAxisComponent(const int& i) const {
        return axes[i];
    }

    void setNumPolygons(int n) {
        numPolygons = n;
    }

    int* getNumPolygons() {
        return &numPolygons;
    }

    int getNumsubConvexPolygons() const {
        return (int)subConvexPolygons.size();
    }

    Polygon* getConvexShape2(const int& i) {
        return &(subConvexPolygons[i]);
    }

    const Polygon& getConvexShape(const int& i) const {
        return subConvexPolygons[i];
    }

    bool isConcavePolygon() const;
    void addConvexPolygon(const Polygon*);
    bool checkCollision(const Polygon*) const;
    bool checkCollisionBC(const Polygon*) const;

    void move(const Point&);
    void rotate(const double&);
    void calculateCentroid();
    void updateLowerVertice(const Point&);

    void updateBB(const Point&);
    double* getBB() {
        return &(BB[0]);
    }
    //int pointSideCurve(const Point&, const Point&, const Point&, const Point&, const Point&);
};

// Particion
class Partition {
private:
    struct VertexPartition {
        bool active;
        bool convex;
        bool isEar;

        Point p;
        double angle;
        VertexPartition* prev;
        VertexPartition* next;

        VertexPartition();
    };

    bool isReflex(Point&, Point&, Point&);
    bool isInsideTriangle(Point&, Point&, Point&, Point&);

    // ??
    Point normalize(const Point&);
    void updateVertex(VertexPartition*, VertexPartition*, const int&);

    bool triangulation(class::Polygon*, std::vector<class::Polygon>*, const std::string&);
    bool convexPartition(class::Polygon*, std::vector<class::Polygon>*, const std::string&);

public:
    bool isConvex(const Point&, const Point&, const Point&);
    void concaveToConvex(class::Polygon*);
};

// Leer fichero dxf
class dxfFilter : public DL_CreationAdapter {
private:
    std::ifstream file;
    // Poligonos originales
    //DynamicArray<class::Polygon> polys;
    // Copia del vector de poligonos, a los que se le aplica mover y rotar piezas
    DynamicArray<class::Polygon> polysCP;
    //std::vector<Polygon*> polys;
    class::Polygon convexHullPolygon;
    //DynamicArray<Point> verticesConvexHull;
    //std::vector<Point> verticesConvexHull;

    std::string readLayer();
    double convexHullArea(std::stack<Point>);
    Point calculateLowerPoint() const;

public:
    DynamicArray<class::Polygon> polys;

    // ----

    dxfFilter(const std::string&);
    ~dxfFilter();

    int getNumPolygons() const {
        return (int)polys.size();
    }

    class::Polygon* getPolygon(const int& i) const {
        // return polys[i];
        return polys.get(i);
    }

    class::Polygon* getPolygons() {
        return polys.getArr();
    }

    void setPointConvexHullPolygon(Point p) {
        convexHullPolygon.setPoint(p);
    }

    Point getPointConvexHullPolygon(const int &i) const {
        return convexHullPolygon[i];
    }

    int getNumPointsConvexHullPolygon() const {
        return convexHullPolygon.getNumPoints();
    }
 
    virtual void addLine(const DL_LineData&);
    virtual void addControlPoint(const DL_ControlPointData&);
    bool collisionDetection() const;
    int collisionDetection2(void***, const int&);
    bool collisionDetectionBB(void**, void**);
    bool collisionDetection1Polygon(void***, const int&, const int&);
    double convexHull();

    void move2(void**, const double &x, const double &y/*const Point&*/);
    void rotate2(void**, const double&);

    /*DynamicArray<Point> getVerticesConvexHull() const {
        return verticesConvexHull;
    }

    int getNumVerticesCH() const {
        return (int)verticesConvexHull.size();
    }

    void addVerticeConvexHull(const Point& p) {
        //verticesConvexHull.push_back(p);
        verticesConvexHull.pushBack(p);
    }

    Point* getVerticeConvexHull(const int& i) const {
        //return verticesConvexHull[i];
        return verticesConvexHull.get(i);
    }*/

    void readPolygons(const std::string&, DL_Dxf*);

    void drawFile(DL_Dxf*, const double&, const std::string&, Area*) const;

    void printPolygons() const;

    // Funcion de evaluacion
    double evFunction(DL_Dxf*, const std::string&, Area*, const Point&, const double&);
    double evFunction2(void***, const int&, const double&, const double&, const double&);
    void transferData(void***&) const;
    void freeData(void***&) const;
};

#endif