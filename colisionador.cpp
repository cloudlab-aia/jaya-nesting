#include <stdio.h>
#include <algorithm>    
#include <complex>
#define _USE_MATH_DEFINES
#define MAXDOUBLE 1e40
#include "Colisionador.h"
/* The value of 1.0 / (1L<<23) is float machine epsilon. */
#ifdef FLOAT_ACCURACY
#define INV_EPS (1L<<23)
#else
/* The value of 1.0 / (1L<<14) is enough for most applications */
#define INV_EPS (1L<<14)
#endif

// Class Area
Area::Area(const std::string& nameFile, DL_Dxf* dxf) {
    std::string line;

    file.open(nameFile, std::ifstream::binary);

    while (getline(file, line)) {
        if (!line.empty())
            line.pop_back();

        if (line == "ENTITIES")
            break;
        //std::cout<<line<<std::endl;
    }

    if (!dxf->in(nameFile, this)) {
        std::cerr << "NO se ha podido leer el archivo" << nameFile << std::endl;
        exit(1);
    }
}

bool Area::isConvex(const Point& a, const Point& b, const Point& c) const {
    Point ab = b.minus(a);
    Point ac = c.minus(a);
    //Point ab = b - a;
    //Point ac = c - a;

    if (ab.crossSign(ac) > 0) {
        return true;
    }

    return false;
}

bool Area::isInsideArea(const Point* cvPoints, const int& numPoints) const {
    for (int i = 0; i < numPoints; i++) {
        for (int j = 0; j < (int)points.size(); j++) {
            if (isConvex(points[j], points[(j - 1 < 0) ? points.size() - 1 : j - 1], cvPoints[i]))
                return false;
        }
    }

    return true;
}

bool Area::isInsideArea(const DynamicArray<Point>& verticesConvexHull) const {
    Point p;
    //Partition pp;

    for (size_t i = 0; i < verticesConvexHull.size(); ++i) {
        p = *verticesConvexHull.get(i);
        //std::cout<<"Punto: "<<p.x<<", "<<p.y<<std::endl;

        for (int j = 0; j < (int)points.size(); ++j) {
            //std::cout<<"p1: "<<points[j].x<<", "<<points[j].y<<" "<<points[(j - 1 < 0) ? points.size() - 1 : j - 1].x<<", "<<points[(j - 1 < 0) ? points.size() - 1 : j - 1].y<<std::endl;
            if (isConvex(points[j], points[(j - 1 < 0) ? points.size() - 1 : j - 1], p))
                return false;
        }
    }

    return true;
}

std::ostream& operator<<(std::ostream& os, const Area& a) {
    os << "Area:\n";

    for (const auto& p : a.points) {
        os << "(" << p.x << ", " << p.y << ") ";
    }

    return os;
}

Area::~Area() {
    file.close();
}

void Area::addLine(const DL_LineData& data) {
    points.push_back(Point(data.x1, data.y1));
}

// Class Polygon

void Polygon::initBB() {
    BB.resize(4);
    BB[0] = MAXDOUBLE;
    BB[1] = MAXDOUBLE;
    BB[2] = -MAXDOUBLE;
    BB[3] = -MAXDOUBLE;
}

Polygon::Polygon() {
    layer = "";
    initBB();
}

Polygon::Polygon(const std::string& layer) {
    this->layer = layer;
    initBB();
}

Polygon::Polygon(const int& numPoints, const std::string& layer) {
    this->layer = layer;

    points.resize(numPoints);
    initBB();
}

void Polygon::clear() {
    points.clear();
    layer = "";
    BB.clear();
}

Polygon::~Polygon() {
    clear();
}

Polygon::Polygon(const Polygon& p) {
    if (this != &p) {
        this->layer = p.layer;

        points.clear();
        this->points = p.points;
        this->projections = p.projections;
        this->axes = p.axes;
        this->lowerPoint = p.lowerPoint;
        this->subConvexPolygons = p.subConvexPolygons;
        this->BB = p.BB;
    }
}

Projection Polygon::project(const Point& axis) const {
    auto p = points[0].x * axis.x + points[0].y * axis.y;
    double min = p;
    double max = min;

    //std::cout<<"Inicial: "<<min<<" "<<max<<std::endl;
    for (/*const auto &c : components*/size_t i = 0; i < points.size() - 1; ++i) {
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
    //std::cout<<"PROJECTION: "<<min<<" - "<<max<<std::endl;

    return Projection(min, max);
}

bool Polygon::isConcavePolygon() const {
    // Si la direccion del vector resultado del producto vectorial de cada par de aristas que
    // conforman la figura es la misma, es convexo
    Point a, b, c;

    /*for(auto v : points)
        std::cout<<v.x<<", "<<v.y<<" ";
    std::cout<<std::endl;*/

    for (int i = 0; i < (int)points.size(); ++i) {
        b = points[i];
        //std::cout<<i - 1<<std::endl;
        //std::cout<<points.back().x<<", "<<points.back().y<<" ?? "<<points[i - 1].x<<", "<<points[i - 1].y<<std::endl;
        (i - 1 < 0) ? a = points.back() : a = points[i - 1];
        c = points[(i + 1) % points.size()];

        Point ab = b.minus(a);
        Point ac = c.minus(a);
        //Point ab = b - a;
        //Point ac = c - a;

        if (ab.crossSign(ac) < 0) {
            // Es concavo
            return true;
        }
    }

    return false;
}

void Polygon::addConvexPolygon(const Polygon* p) {
    subConvexPolygons.push_back(*p);

    for (int i = 0; i < p->getNumPoints(); ++i) {
        // Calculamos para cada airsta, su eje de proyeccion (igual que al final de addLine())
        double vX = p->getPoint((i + 1) % p->getNumPoints()).x - p->getPoint(i).x;
        double vY = p->getPoint((i + 1) % p->getNumPoints()).y - p->getPoint(i).y;
        // Obtenemos el vector perpendicular a este (= eje de proyeccion)
        // NO estan normalizados
        subConvexPolygons.back().setAxis(Point(-vY, vX));
        // Obtenemos la proyeccion
        subConvexPolygons.back().setProjection(
            subConvexPolygons.back().project(
                subConvexPolygons.back().getAxisComponent(i)
            )
        );
    }
}

void Polygon::updateBB(const Point &p) {
    if (BB[0] == MAXDOUBLE || p.x > BB[0]) {
        BB[0] = p.x;
    }

    if (BB[1] == MAXDOUBLE || p.y > BB[1]) {
        BB[1] = p.y;
    }

    if (BB[2] == -MAXDOUBLE || p.x < BB[2]) {
        BB[2] = p.x;
    }

    if (BB[3] == -MAXDOUBLE || p.y < BB[3]) {
        BB[3] = p.y;
    }
}

void sortSpecial(double* a) {
    bool flip;
    double temp;

    do {
        flip = false;
        for (int i = 0; i < 3 - 1; ++i) {
            if ((a[i + 1] >= 0 && a[i] > a[i + 1]) || (a[i] < 0 && a[i + 1] >= 0)) {
                flip = true;
                temp = a[i];
                a[i] = a[i + 1];
                a[i + 1] = temp;

            }
        }
    } while (flip);

    //return a;
}

int sgn(double x) {
    if (x < 0.0)
        return -1;

    return 1;
}

double* cubicRoots(double* P) {
    double a, b, c, d, A, B, C;
    double Q, R, D, S, T, Im;
    double* t;
    double DQ;
    double th;

    a = P[0];
    b = P[1];
    c = P[2];
    d = P[3];

    A = b / a;
    B = c / a;
    C = d / a;

    Q = (3 * B - std::pow(A, 2.0)) / 9;
    R = (9 * A * B - 27 * C - 2 * std::pow(A, 3.0)) / 54;
    D = std::pow(Q, 3.0) + std::pow(R, 2.0);        // polynomial discriminant
#ifdef PRINTDEBUG
    std::cout<<A <<" " << B << " " << C << " " << Q << " " << R << " " << D<<std::endl;
#endif
    t = new double[3];

    if (P[0] == 0) {
        if (P[1] == 0) {
#ifdef PRINTDEBUG
            std::cout << "Linear formula detected" << std::endl;
#endif
            //var t = Array();
            t[0] = -1 * (P[3] / P[2]);
            t[1] = -1;
            t[2] = -1;

            /*discard out of spec roots*/
            for (int i = 0; i < 1; ++i)
                if (t[i] < 0 || t[i] > 1.0)
                    t[i] = -1;

            /*sort but place -1 at the end*/
            sortSpecial(t);

            //console.log(t[0]);
            return t;
        }

        //console.log("Quadratic formula detected");
#ifdef PRINTDEBUG
        std::cout << "Quadratic formula detected" << std::endl;
#endif

        DQ = std::pow(P[2], 2.0) - 4 * P[1] * P[3]; // quadratic discriminant
        if (DQ >= 0) {
            DQ = std::sqrt(DQ);
            //var t = Array();
            t[0] = -1 * ((DQ + P[2]) / (2 * P[1]));
            t[1] = ((DQ - P[2]) / (2 * P[1]));
            t[2] = -1;

            /*discard out of spec roots*/
            for (int i = 0; i < 2; ++i)
                if (t[i] < 0 || t[i] > 1.0)
                    t[i] = -1;

            /*sort but place -1 at the end*/
            sortSpecial(t);

            // console.log(t[0]+" "+t[1]);
#ifdef PRINTDEBUG
            std::cout << t[0] << " " << t[1] << std::endl;
#endif
            return t;
        }
    }

    if (D >= 0) {                                 // complex or duplicate roots
        S = sgn(R + std::sqrt(D)) * std::pow(std::abs(R + std::sqrt(D)), (1.0 / 3.0));
        T = sgn(R - std::sqrt(D)) * std::pow(std::abs(R - std::sqrt(D)), (1.0 / 3.0));

        t[0] = -A / 3 + (S + T);                    // real root
        t[1] = -A / 3 - (S + T) / 2;                  // real part of complex root
        t[2] = -A / 3 - (S + T) / 2;                  // real part of complex root
        Im = std::abs(std::sqrt(3) * (S - T) / 2);    // complex part of root pair   

        /*discard complex roots*/
        if (Im != 0) {
            t[1] = -1;
            t[2] = -1;
        }

    }
    else {                                          // distinct real roots
        th = std::acos(R / std::sqrt(-std::pow(Q, 3.0)));

        t[0] = 2 * std::sqrt(-Q) * std::cos(th / 3) - A / 3;
        t[1] = 2 * std::sqrt(-Q) * std::cos((th + 2 * M_PI) / 3) - A / 3;
        t[2] = 2 * std::sqrt(-Q) * std::cos((th + 4 * M_PI) / 3) - A / 3;
        Im = 0.0;
    }

    /*discard out of spec roots*/
    for (int i = 0; i < 3; ++i)
        if (t[i] < 0 || t[i] > 1.0)
            t[i] = -1;

    /*sort but place -1 at the end*/
    sortSpecial(t);

    // console.log(t[0]+" "+t[1]+" "+t[2]);
#ifdef PRINTDEBUG
    std::cout << t[0] << " " << t[1] << " " << t[2] << std::endl;
#endif
    return t;
}

double* bezierCoeffs(double P0, double P1, double P2, double P3) {
    double* Z = new double[4];

    Z[0] = -P0 + 3 * P1 + -3 * P2 + P3;
    Z[1] = 3 * P0 - 6 * P1 + 3 * P2;
    Z[2] = -3 * P0 + 3 * P1;
    Z[3] = P0;

    return Z;
}

bool Polygon::findIntersection(const int &indexBC, const Point &pp1, const Point &pp2) const {
    double A, B, C;
    double P[4], X[2];
    double* bx, * by;
    double* r;
    double s, t;
    double* I = new double[3];
    bool intersection = false;

    A = pp2.y - pp1.y;
    B = pp1.x - pp2.x;
    C = pp1.x * (pp1.y - pp2.y) + pp1.y * (pp2.x - pp1.x);

    bx = bezierCoeffs(pointsBC[indexBC].x, pointsBC[indexBC + 1].x, pointsBC[indexBC + 2].x, pointsBC[indexBC + 3].x);
    by = bezierCoeffs(pointsBC[indexBC].y, pointsBC[indexBC + 1].y, pointsBC[indexBC + 2].y, pointsBC[indexBC + 3].y);
#ifdef PRINTDEBUG
    std::cout << bx[0] << ", " << bx[1] << ", " << bx[2] << ", " << bx[3] << " - " << by[0] << ", " << by[1] << ", " << by[2] << ", " <<
        by[3] << ", " << std::endl;
#endif
    P[0] = A * bx[0] + B * by[0];		/*t^3*/
    P[1] = A * bx[1] + B * by[1];		/*t^2*/
    P[2] = A * bx[2] + B * by[2];		/*t*/
    P[3] = A * bx[3] + B * by[3] + C;   /*1*/
#ifdef PRINTDEBUG
    std::cout << P[0] << " " << P[1] << " " << P[2] << " " << P[3] << " " << std::endl;
#endif
    r = cubicRoots(P);

    // Pueden existir como máximo 3 intersecciones entre una línea y una curba cúbica de bezier
    for (int i = 0; !intersection && i < 3; ++i) {
        t = r[i];

        X[0] = bx[0] * t * t * t + bx[1] * t * t + bx[2] * t + bx[3];
        X[1] = by[0] * t * t * t + by[1] * t * t + by[2] * t + by[3];

        if ((pp2.x - pp1.x) != 0)           /*if not vertical line*/
            s = (X[0] - pp1.x) / (pp2.x - pp1.x);
        else
            s = (X[1] - pp1.y) / (pp2.y - pp1.y);

        /*in bounds?*/
        if (t < 0 || t > 1.0 || s < 0 || s > 1.0 || isnan(X[0]) || isnan(X[1])) {
            //X[0] = -100;  /*move off screen*/
            //X[1] = -100;
            X[0] = DBL_MAX;
            X[1] = DBL_MAX;
        }
        else {
            intersection = true;
        }

        /*move intersection point*/
        // I[i].setAttributeNS(NULL, "cx", X[0]);
        // I[i].setAttributeNS(NULL, "cy", X[1]);
#ifdef PRINTDEBUG
        std::cout << "Puntos de intersección:\n";
        std::cout << "X: " << X[0] << " Y: " << X[1] << std::endl;
#endif
    }

    return intersection;
}

// Clase Bezier
Bezier::Bezier(Point p0, Point p1, Point p2, Point p3) {
    this->p0 = p0;
    //this->p0->refcount += 8;
    this->p1 = p1;
    //this->p1->refcount += 8;
    this->p2 = p2;
    //this->p2->refcount += 8;
    this->p3 = p3;
    //this->p3->refcount += 8;
}

Bezier::Bezier() {

}
 
Bezier::~Bezier() {

}

Bezier& Bezier::operator=(const Bezier& b) {
    if (this != &b) {
        p0 = b.p0;
        p1 = b.p1;
        p2 = b.p2;
        p3 = b.p3;
    }

    return *this;
}

double alignCurveX(const double& t, const Bezier& curve) {
    return (
        curve.p0.x * (1 - t) * (1 - t) * (1 - t) +
        3 * curve.p1.x * t * (1 - t) * (1 - t) +
        3 * curve.p2.x * t * t * (1 - t) +
        curve.p3.x * t * t * t
        );
}

double alignCurveY(const double& t, const Bezier& curve) {
    return (
        curve.p0.y * (1 - t) * (1 - t) * (1 - t) +
        3 * curve.p1.y * t * (1 - t) * (1 - t) +
        3 * curve.p2.y * t * t * (1 - t) +
        curve.p3.y * t * t * t
        );
}

void computeBB(const Bezier& curve, double& xl, double& xh, double& yl, double& yh) {
    // P'(t) = (-3*p0 + 9p1 - 9p2 + 3p3)*t^2 + (6*p0 - 12p1 + 6p2)*t + (-3p0 + 3p1)
    // t = (-b +- raiz(b^2 - 4ac))/2a -> t0, t1, t2, t3 -> [0, 1]
    //double xl, xh, yl, yh;
    double t0 = -1, t1 = -1, t2 = -1, t3 = -1;
    // X
    double a = -3 * curve.p0.x + 9 * curve.p1.x - 9 * curve.p2.x + 3 * curve.p3.x;
    double b = 6 * curve.p0.x - 12 * curve.p1.x + 6 * curve.p2.x;
    double c = -3 * curve.p0.x + 3 * curve.p1.x;
    double r = b * b - 4 * a * c;

    if (curve.p0.x < curve.p3.x) {
        xl = curve.p0.x;
        xh = curve.p3.x;
    }
    else {
        xl = curve.p3.x;
        xh = curve.p0.x;
    }

    if (r >= 0) {
        t0 = (-b + std::sqrt(r)) / (2 * a);

        if (t0 > 0 && t0 < 1) {
            // Alineamos al eje la curva (apartado Aligning curves)
            double x = alignCurveX(t0, curve);

            if (x < xl)
                xl = x;

            if (x > xh)
                xh = x;
        }

        t1 = (-b - std::sqrt(r)) / (2 * a);

        if (t1 > 0 && t1 < 1) {
            double x = alignCurveX(t1, curve);

            if (x < xl)
                xl = x;

            if (x > xh)
                xh = x;
        }
    }

    // Y
    a = -3 * curve.p0.y + 9 * curve.p1.y - 9 * curve.p2.y + 3 * curve.p3.y;
    b = 6 * curve.p0.y - 12 * curve.p1.y + 6 * curve.p2.y;
    c = -3 * curve.p0.y + 3 * curve.p1.y;
    r = b * b - 4 * a * c;

    if (curve.p0.y < curve.p3.y) {
        yl = curve.p0.y;
        yh = curve.p3.y;
    }
    else {
        yl = curve.p3.y;
        yh = curve.p0.y;
    }

    if (r >= 0) {
        t2 = (-b + std::sqrt(r)) / (2 * a);

        if (t2 > 0 && t2 < 1) {
            double y = alignCurveY(t2, curve);

            if (y < yl)
                yl = y;

            if (y > yh)
                yh = y;
        }

        t3 = (-b - std::sqrt(r)) / (2 * a);

        if (t3 > 0 && t3 < 1) {
            double y = alignCurveY(t3, curve);

            if (y < yl)
                yl = y;

            if (y > yh)
                yh = y;
        }
    }
}

bool intersectBB(const Bezier& a, const Bezier& b) {
    double axl, axh, ayl, ayh, bxl, bxh, byl, byh;

    computeBB(a, axl, axh, ayl, ayh);
    computeBB(b, bxl, bxh, byl, byh);

    if (axl > bxh || ayl > byh || bxl > axh || byl > ayh) {
        return false;
    }

    return true;
}

double log4(double x) {
    return 0.5 * log2(x);
}

void Bezier::splitCurve(Bezier &a, Bezier &b) const {
    // a -> parte izquierda
    // b -> parte derecha
    //std::cout<<"es aqui??"<<std::endl;
    a.p0 = p0;
    b.p3 = p3;
    a.p1 = Point(p0, p1);
    //std::cout<<"sii "<<p2->x<<" "<<p2->y<<" "<<p3->x<<" "<<p3->y<<std::endl;
    b.p1 = Point(p1, p2);
    //std::cout<<"nose"<<std::endl;
    b.p2 = Point(p2, p3);
    //std::cout<<"????"<<std::endl;
    a.p2 = Point(a.p1, b.p1);
    b.p1 = mid(b.p1, b.p2);
    a.p3 = Point(a.p2, b.p1);
    b.p0 = Point(a.p2, b.p1);
}

bool ok = false;

bool recursivelyIntersect(Bezier a, double t0, double t1, int deptha, Bezier b, double u0, double u1,
    int depthb/*, double** parameters, int& index*/) {
    Bezier a1, a2, b1, b2;
#ifdef PRINTDEBUG
    std::cout << "empiezo" << std::endl;
    std::cout << deptha << " " << depthb << std::endl;
#endif
    if (deptha > 0) {
        a.splitCurve(a1, a2);

        double tmid = (t0 + t1) * 0.5;

        --deptha;

        if (depthb > 0) {
            b.splitCurve(b1, b2);
            double umid = (u0 + u1) * 0.5;

            --depthb;

            if (!ok && intersectBB(a1, b1)) {
                recursivelyIntersect(a1, t0, tmid, deptha, b1, u0, umid, depthb/*, parameters, index*/);
            }
#ifdef PRINTDEBUG
            std::cout << "????--?" << std::endl;
#endif
            if (!ok && intersectBB(a2, b1))
                recursivelyIntersect(a2, tmid, t1, deptha, b1, u0, umid, depthb/*, parameters, index*/);
#ifdef PRINTDEBUG
            std::cout << ":...." << std::endl;
#endif
            if (!ok && intersectBB(a1, b2))
                recursivelyIntersect(a1, t0, tmid, deptha, b2, umid, u1, depthb/*, parameters, index*/);

            if (!ok && intersectBB(a2, b2))
                recursivelyIntersect(a2, tmid, t1, deptha, b2, umid, u1, depthb/*, parameters, index*/);
#ifdef PRINTDEBUG
            std::cout << "....??" << std::endl;
#endif
        }
        else {
#ifdef PRINTDEBUG
            std::cout << "00" << std::endl;
#endif
            if (!ok && intersectBB(a1, b))
                recursivelyIntersect(a1, t0, tmid, deptha, b, u0, u1, depthb/*, parameters, index*/);

            if (!ok && intersectBB(a2, b))
                recursivelyIntersect(a2, tmid, t1, deptha, b, u0, u1, depthb/*, parameters, index*/);
        }
    }
    else {
#ifdef PRINTDEBUG
        std::cout << "fin?????" << std::endl;
#endif
        if (depthb > 0) {
            b.splitCurve(b1, b2);
            double umid = (u0 + u1) * 0.5;

            --depthb;

            if (!ok && intersectBB(a, b1))
                recursivelyIntersect(a, t0, t1, deptha, b1, u0, umid, depthb/*, parameters, index*/);

            if (!ok && intersectBB(a, b2))
                recursivelyIntersect(a, t0, t1, deptha, b2, umid, u1, depthb/*, parameters, index*/);
        }
        else {
            // Ya no podemos subdividir mas las curvas
            double xlk = a.p3.x - a.p0.x;
            double ylk = a.p3.y - a.p0.y;
            double xnm = b.p3.x - b.p0.x;
            double ynm = b.p3.y - b.p0.y;
            double xmk = b.p0.x - a.p0.x;
            double ymk = b.p0.y - a.p0.y;
            double det = xnm * ylk - ynm * xlk;

            if (1.0 + det == 1.0) {
                return false;
            }
            else {
                double detinv = 1.0 / det;
                double s = (xnm * ymk - ynm * xmk) * detinv;
                double t = (xlk * ymk - ylk * xmk) * detinv;

                // Esto no se si puede ocurrir, pero serian valores NO validos
                if ((s < 0.0) || (s > 1.0) || (t < 0.0) || (t > 1.0))
                    return false;

                // Hay interseccion
                //parameters[0][index] = t0 + s * (t1 - t0);
                //parameters[1][index] = u0 + t * (u1 - u0);
                //++index;
                ok = true;
                return ok;
            }
        }
    }
#ifdef PRINTDEBUG
    std::cout << "......... " << deptha << " " << depthb << std::endl;
#endif
    return ok;
}

bool findIntersections(const Bezier a, const Bezier b) {
    //double** parameters = new double* [2];
    //int index = 0;

    // Maximo nueve intersecciones por curva
    //parameters[0] = new double[9];
    //parameters[1] = new double[9];

    if (intersectBB(a, b)) {
        // No se solapan los BB
        // Arreglar warnning
        Vector la1 = vabs((a.p2 - a.p1) - (a.p1 - a.p0));
        Vector la2 = vabs((a.p3 - a.p2) - (a.p2 - a.p1));
        Vector la;

        if (la1.getX() > la2.getX())
            la.setX(la1.getX());
        else
            la.setX(la2.getX());

        if (la1.getY() > la2.getY())
            la.setY(la1.getY());
        else
            la.setY(la2.getY());

        Vector lb1 = vabs((b.p2 - b.p1) - (b.p1 - b.p0));
        Vector lb2 = vabs((b.p3 - b.p2) - (b.p2 - b.p1));
        Vector lb;

        if (lb1.getX() > lb2.getX())
            lb.setX(lb1.getX());
        else
            lb.setX(lb2.getX());

        if (lb1.getY() > lb2.getY())
            lb.setY(lb1.getY());
        else
            lb.setY(lb2.getY());

        double l0;

        if (la.getX() > la.getY())
            l0 = la.getX();
        else
            l0 = la.getY();

        int ra;

        if (l0 * 0.75 * M_SQRT2 + 1.0 == 1.0)
            ra = 0;
        else
            ra = (int)ceil(log4(M_SQRT2 * 6.0 / 8.0 * INV_EPS * l0));

        if (lb.getX() > lb.getY())
            l0 = lb.getX();
        else
            l0 = lb.getY();

        int rb;

        if (l0 * 0.75 * M_SQRT2 + 1.0 == 1.0)
            rb = 0;
        else
            rb = (int)ceil(log4(M_SQRT2 * 6.0 / 8.0 * INV_EPS * l0));


        // (cutva a, inicio curva a, final curva a, profundidad de subdivision restante de a, ..., puntos de 
        // interseccion, indice para parameters)
        return recursivelyIntersect(a, 0.0, 1.0, ra, b, 0.0, 1.0, rb/*, parameters, index*/);
    }

    //if (index > 9)
    //    parameters[0][index] = parameters[1][index] = -1.0;

    // Ordenamos las intersecciones siguiendo la curva desde t = 0, hasta t = 1
    //sort(parameters[0], index);
    //sort(parameters[1], index);

    //return parameters;
    return false;
}

bool Bezier::intersect(Bezier B) {
    /*B.p0->refcount++;
    B.p1->refcount++;
    B.p2->refcount++;
    B.p3->refcount++;*/

    double xl, xh, yl, yh;
    computeBB(*this, xl, xh, yl, yh);
#ifdef PRINTDEBUG
    std::cout << xl << " " << xh << " " << yl << " " << yh << std::endl;
#endif
    //Bezier** rvalue = new Bezier * [2];

    // Como maximo pueden existir 9 intersecciones de una curva con otra (sin contar las intersecciones con ella 
    // misma), por tanto, se pueden generar como maximo 10 segmentos de una curva
    // Esta curva
    //rvalue[0] = new Bezier[10];
    // Curva B
    //rvalue[1] = new Bezier[10];

    // Devuelve las intersecciones de en cada una de las dos curvas (por separado) en el rango [0, 1] 
    // t[2][9] (2 curvas / 9 intersecciones en cada curva como maximo)
    return findIntersections(*this, B);
}

int pointSideCurve(const Point &p0, const Point &p1, const Point &p2, const Point &p3, const Point &p) {
    double a3 = p3.x - p0.x - 3 * p2.x + 3 * p1.x;
    double a2 = -6 * p1.x + 3 * p0.x + 3 * p2.x;
    double a1 = -3 * p0.x + 3 * p1.x;
    double a0 = p0.x;
    double b3 = p3.y - p0.y - 3 * p2.y + 3 * p1.y;
    double b2 = -6 * p1.y + 3 * p0.y + 3 * p2.y;
    double b1 = -3 * p0.y + 3 * p1.y;
    double b0 = p0.y;

    double vxxx = b3 * b3 * b3;
    double vxxy = -3 * a3 * b3 * b3;
    double vxyy = 3 * b3 * a3 * a3;
    double vyyy = -a3 * a3 * a3;
    double vxx = -3 * a3 * b1 * b2 * b3 + a1 * b2 * b3 * b3 - a2 * b3 * b2 * b2 + 2 * a2 * b1 * b3 * b3 + 3 * a3 * b0 * b3 * b3 +
        a3 * b2 * b2 * b2 - 3 * a0 * b3 * b3 * b3;
    double vxy = a1 * a3 * b2 * b3 - a2 * a3 * b1 * b3 - 6 * b0 * b3 * a3 * a3 - 3 * a1 * a2 * b3 * b3 - 2 * a2 * a3 * b2 * b2 +
        2 * b2 * b3 * a2 * a2 + 3 * b1 * b2 * a3 * a3 + 6 * a0 * a3 * b3 * b3;
    double vyy = 3 * a1 * a2 * a3 * b3 + a3 * b2 * a2 * a2 - a2 * b1 * a3 * a3 - 3 * a0 * b3 * a3 * a3 - 2 * a1 * b2 * a3 * a3 -
        b3 * a2 * a2 * a2 + 3 * b0 * a3 * a3 * a3;
    double vx = a2 * a3 * b0 * b1 * b3 - a1 * a2 * b1 * b2 * b3 - a1 * a3 * b0 * b2 * b3 + 6 * a0 * a3 * b1 * b2 * b3 +
        b1 * a1 * a1 * b3 * b3 + b3 * a2 * a2 * b1 * b1 + 3 * b3 * a3 * a3 * b0 * b0 + a1 * a3 * b1 * b2 * b2 -
        a2 * a3 * b2 * b1 * b1 - 6 * a0 * a3 * b0 * b3 * b3 - 4 * a0 * a2 * b1 * b3 * b3 - 3 * b0 * b1 * b2 * a3 * a3 -
        2 * a0 * a1 * b2 * b3 * b3 - 2 * a1 * a3 * b3 * b1 * b1 - 2 * b0 * b2 * b3 * a2 * a2 + 2 * a0 * a2 * b3 * b2 * b2 +
        2 * a2 * a3 * b0 * b2 * b2 + 3 * a1 * a2 * b0 * b3 * b3 + a3 * a3 * b1 * b1 * b1 + 3 * a0 * a0 * b3 * b3 * b3 -
        2 * a0 * a3 * b2 * b2 * b2;
    double vy = a0 * a2 * a3 * b1 * b3 + a1 * a2 * a3 * b1 * b2 - a0 * a1 * a3 * b2 * b3 - 6 * a1 * a2 * a3 * b0 * b3 -
        a1 * a1 * a1 * b3 * b3 - 3 * a3 * a3 * a3 * b0 * b0 - a1 * a3 * a3 * b1 * b1 - a3 * a1 * a1 * b2 * b2 -
        3 * a3 * a0 * a0 * b3 * b3 + a2 * b2 * b3 * a1 * a1 - a1 * b1 * b3 * a2 * a2 - 3 * a0 * b1 * b2 * a3 * a3 -
        2 * a0 * b2 * b3 * a2 * a2 - 2 * a3 * b0 * b2 * a2 * a2 + 2 * a0 * a2 * a3 * b2 * b2 + 2 * a2 * b0 * b1 * a3 * a3 +
        2 * a3 * b1 * b3 * a1 * a1 + 3 * a0 * a1 * a2 * b3 * b3 + 4 * a1 * b0 * b2 * a3 * a3 + 6 * a0 * b0 * b3 * a3 * a3 +
        2 * b0 * b3 * a2 * a2 * a2;
    double v0 = a0 * a1 * a2 * b1 * b2 * b3 + a0 * a1 * a3 * b0 * b2 * b3 - a0 * a2 * a3 * b0 * b1 * b3 - a1 * a2 * a3 * b0 * b1 * b2 +
        b0 * a1 * a1 * a1 * b3 * b3 - b3 * a2 * a2 * a2 * b0 * b0 + a1 * b0 * a3 * a3 * b1 * b1 + a1 * b2 * a0 * a0 * b3 * b3 +
        a3 * b0 * a1 * a1 * b2 * b2 + a3 * b2 * a2 * a2 * b0 * b0 - a0 * b1 * a1 * a1 * b3 * b3 - a0 * b3 * a2 * a2 * b1 * b1 -
        a2 * b1 * a3 * a3 * b0 * b0 - a2 * b3 * a0 * a0 * b2 * b2 - 3 * a0 * b3 * a3 * a3 * b0 * b0 - 2 * a1 * b2 * a3 * a3 * b0 * b0 +
        2 * a2 * b1 * a0 * a0 * b3 * b3 + 3 * a3 * b0 * a0 * a0 * b3 * b3 + a0 * a2 * a3 * b2 * b1 * b1 + a1 * b0 * b1 * b3 * a2 * a2 -
        a0 * a1 * a3 * b1 * b2 * b2 - a2 * b0 * b2 * b3 * a1 * a1 - 3 * a0 * a1 * a2 * b0 * b3 * b3 - 3 * a3 * b1 * b2 * b3 * a0 * a0 -
        2 * a0 * a2 * a3 * b0 * b2 * b2 - 2 * a3 * b0 * b1 * b3 * a1 * a1 + 2 * a0 * a1 * a3 * b3 * b1 * b1 +
        2 * a0 * b0 * b2 * b3 * a2 * a2 + 3 * a0 * b0 * b1 * b2 * a3 * a3 + 3 * a1 * a2 * a3 * b3 * b0 * b0 + a3 * a3 * a3 * b0 * b0 * b0 -
        a0 * a0 * a0 * b3 * b3 * b3 + a3 * a0 * a0 * b2 * b2 * b2 - a0 * a3 * a3 * b1 * b1 * b1;

    double implicitCurveEq = vxxx * p.x * p.x * p.x + vxxy * p.x * p.x * p.y + vxyy * p.x * p.y * p.y + vyyy * p.y * p.y * p.y +
        vxx * p.x * p.x + vxy * p.x * p.y + vyy * p.y * p.y + vx * p.x + vy * p.y + v0;

#ifdef PRINTDEBUG
	std::cout << p0.x << ", " << p0.y << " - " << p1.x << ", " << p1.y << " - " << p2.x << ", " << p2.y << " - " << p3.x << ", " << 
        p3.y << std::endl;
    std::cout << p.x<<","<<p.y << " - " << implicitCurveEq << std::endl;
#endif

    return (implicitCurveEq >= 0) ? -1 : 0;
}

bool Polygon::checkCollisionBC(const Polygon* p) const {
    // Colisión entre splines y aristas u otros splines
    for (size_t i = 0; i < pointsBC.size(); i += 4) {
        // Spline con aristas
        for (size_t j = 0; j < p->getNumPoints(); ++j) {
            Point pp1 = p->getPoint(j);
            Point pp2 = p->getPoint((j + 1) % p->getNumPoints());
            if ( findIntersection(i, p->getPoint(j), p->getPoint((j + 1) % p->getNumPoints())) || 
                !findIntersection(i, p->getPoint(j), centroid) ) {
                return true;
            }
        }
#ifdef PRINTDEBUG
        std::cout << *this << std::endl;
#endif
        // Colisión entre splines
        Bezier c1(getBCPoint(i), getBCPoint(i + 1), getBCPoint(i + 2), getBCPoint(i + 3));

        for (size_t j = 0; j < p->getNumPointsBC(); j += 4) {
            Bezier c2(p->getBCPoint(j), p->getBCPoint(j + 1), p->getBCPoint(j + 2), p->getBCPoint(j + 3));

            if (c1.intersect(c2)) {
                return true;
            }
        }
    }

    return false;
}

bool isSpline(const Polygon &p, const Point& p1) {
    for (int i = 0; i < p.getNumPointsBC(); i += 4) {
        if (p1 == p.getBCPoint(i))
            return true;
    }

    return false;
}

Projection project2(const Point& axis, Point* points, const int& numPoints) {
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

bool checkCollision2(void** pol1, const int &index1, void** pol2, const int& index2) {
    // Num pol + centroide + lowerPoint + BB + Num convex poly * (Num pts + puntos + axes + proyecciones)
    int numPointsPol1 = *(int*)(pol1[index1 + 4]);
    Projection* p1 = (Projection*)(pol1[index1 + 7]);
    Point* axes = (Point*)(pol1[index1 + 6]);
    int numPointsPol2 = *(int*)(pol2[index2 + 4]);
    Point* points = (Point*)(pol2[index2 + 5]);
    for (int i = 0; i < numPointsPol1; ++i) {
        Projection p2 = project2(axes[i], points, numPointsPol2);

        if (p2.min > p1[i].max || p2.max < p1[i].min)
            return false;
    }

    return true;
}

bool Polygon::checkCollision(const Polygon* p) const {
    bool collision = true;
    //std::cout<<"???"<<std::endl;
    // NOTA: para que existe colision, debe existir solapamiento sobre todos los ejes de proyeccion (de ambas piezas)
    
    // Colisión entre aristas
    if (points.size() > 2) {
        for (size_t i = 0; i < points.size(); ++i) {
            if (!isSpline(*this, this->getPoint(i))) {
                // Obtenemos la proyeccion de la pieza actual sobre el eje de proyeccion de c
                //std::cout<<projections.size()<<std::endl;
                Projection p1 = projections[i];

                // Proyectamos s sobre el eje de proyeccion de c
                //std::cout<<"Axis: "<<vertices[i].getX()<<" "<<vertices[i].getY()<<std::endl;
                Projection p2 = p->project(axes[i]/*components[c]->getAxis()*/);

                // Comprobamos si ambas proyecciones se solapan
                if (!p1.overlap(p2)) {
                    return false;
                }
            }
        }
    }
    else {
        return false;
    }

    return collision;
}

void Polygon::updateLowerVertice(const Point& p) {
    if (points.size() == 0) {
        lowerPoint = p;

        return;
    }

    if (p < lowerPoint /*(p.y < lowerPoint.y) || (p.y == lowerPoint.y && p.x < lowerPoint.y)*/) {
        lowerPoint = p;
    }
}

// Movemos un solo polígono
// numPolygons = Polígono + Nº sub polígonos (si los tiene)
void dxfFilter::move2(void** polygon, const double &x, const double &y/*const Point& p*/) {
    Point p(x, y);
    // Número polígonos 
    int numPols = *(int*)(polygon[0]);
    // Centroide
    Point* cent = (Point*)(polygon[1]);
    cent->x += p.x;
    cent->y += p.y;
    // Resto
    int index = 2;
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
            aux = project2(axes[j], points, numPts);
            proj[j].min = aux.min;
            proj[j].max = aux.max;
        }
    }
}

void Polygon::move(const Point& p) {
    // Actualizamos el array de puntos
    for (size_t i = 0; i < points.size(); ++i) {
        //std::cout << points[i].x << " " << points[i].y << " + " << p.x << " " << p.y << std::endl;
        points[i] = points[i] + p;
    }

    // Actualizamos el array de puntos de splines
    for (size_t i = 0; i < pointsBC.size(); ++i) {
        pointsBC[i] = pointsBC[i] + p;
    }

    // Actualizamos el array de proyecciones del poligono principal
    for (size_t i = 0; i < points.size(); ++i) {
        setProjection(project(axes[i]), i);
    }

    // Actualizamos el centroide
    centroid = centroid + p;

    // Actualizamos el punto mas bajo
    //lowerPoint = lowerPoint + p;

    // Si es concavo, actualizamos el array de puntos y proyecciones tambien de sus subpiezas convexas
    if (subConvexPolygons.size() > 0) {
        for (size_t i = 0; i < subConvexPolygons.size(); ++i) {
            for (size_t j = 0; j < subConvexPolygons[i].points.size(); ++j) {
                subConvexPolygons[i].points[j] = subConvexPolygons[i].points[j] + p;

                subConvexPolygons[i].setProjection(subConvexPolygons[i].project(subConvexPolygons[i].axes[j]), j);
            }
        }
    }
}

void Polygon::calculateCentroid() {
    Point ct;

    for (size_t i = 0; i < points.size(); ++i) {
        ct.x += points[i].x;
        ct.y += points[i].y;
    }

    ct.x /= points.size();
    ct.y /= points.size();

    centroid = ct;
}

void dxfFilter::rotate2(void** polygon, const double& angle) {
    Point p, newP, p1, p2;
    double vX, vY;
    // Num polígonos
    int numPols = *(int*)(polygon[0]);
    // Centroide
    Point* cent = (Point*)(polygon[1]);
    int index = 2;
    int numPuntos = 0;
    for (int i = 0; i < numPols; ++i) {
        // Num puntos
        numPuntos = *(int*)(polygon[index++]);
        // Puntos
        Point* points = (Point*)(polygon[index++]);
        // Axes
        Point* axes = (Point*)(polygon[index++]);
        // Proyecciones
        Projection* proj = (Projection*)(polygon[index++]);
        for(int j = 0; j < numPuntos; ++j){
            p.x = points[j].x - cent->x;
            p.y = points[j].y - cent->y;

            newP.x = p.x * cos(angle) - p.y * sin(angle);
            newP.y = p.y * cos(angle) + p.x * sin(angle);

            points[j].x = newP.x + cent->x;
            points[j].y = newP.y + cent->y;

            // Recalculamos los ejes de proyeccion (axes)
            p1.x = points[j].x;
            p1.y = points[j].y;
            p2.x = points[(j + 1) % numPuntos].x;
            p2.y = points[(j + 1) % numPuntos].y;
            vX = p2.x - p1.x;
            vY = p2.y - p1.y;
            axes[j] = Point(-vY, vX);

            // Recalculamos projecciones
            proj[j] = project2(axes[j], points, numPuntos);
        }
    }
}

void Polygon::rotate(const double& angle) {
    // x1 = x0*cos(a) - y0*sen(a)
    // y1 = y0*cos(a) + x0*sen(a)
    Point p, newP, p1, p2;
    double vX, vY;
    // Recalculamos lowerPoint
    lowerPoint = centroid;

    // Recalculamos puntos splines (si hay splines en el polígono)
    for (size_t i = 0; i < pointsBC.size(); ++i) {
        p = pointsBC[i].minus(centroid);
        //p = pointsBC[i] - centroid;
        newP = Point(p.x * cos(angle) - p.y * sin(angle), p.y * cos(angle) + p.x * sin(angle));
        pointsBC[i] = newP + centroid;
    }

    for (size_t i = 0; i < points.size(); ++i) {
        // Recalculamos los puntos
        p = points[i].minus(centroid);
        //p = points[i] - centroid;
        newP = Point(p.x * cos(angle) - p.y * sin(angle), p.y * cos(angle) + p.x * sin(angle));
        points[i] = newP + centroid;

        if (points[i] < lowerPoint)
            lowerPoint = points[i];

        // Recalculamos los ejes de proyeccion (axes)
        p1 = points[i];
        p2 = points[(i + 1) % (int)points.size()];

        vX = p2.x - p1.x;
        vY = p2.y - p1.y;

        axes[i] = Point(-vY, vX);

        // Recalculamos los ejes de proyeccion
        setProjection(project(axes[i]), i);
    }

    if (subConvexPolygons.size() > 0) {
        for (size_t i = 0; i < subConvexPolygons.size(); ++i) {
            for (int j = 0; j < (int)subConvexPolygons[i].points.size(); ++j) {
                // Recalculamos los puntos
                p = subConvexPolygons[i].points[j].minus(centroid);
                // p = subConvexPolygons[i].points[j] - centroid;
                newP = Point(p.x * cos(angle) - p.y * sin(angle), p.y * cos(angle) + p.x * sin(angle));
                subConvexPolygons[i].points[j] = newP + centroid;

                // Recalculamos los ejes de proyeccion (axes)
                p1 = subConvexPolygons[i].points[j];
                p2 = subConvexPolygons[i].points[(j + 1) % (int)subConvexPolygons[i].points.size()];

                vX = p2.x - p1.x;
                vY = p2.y - p1.y;

                subConvexPolygons[i].axes[j] = Point(-vY, vX);

                // Recalculamos los ejes de proyeccion
                subConvexPolygons[i].setProjection(subConvexPolygons[i].project(subConvexPolygons[i].axes[j]), j);
            }
        }
    }
}

/*int Polygon::pointSideCurve(const Point& p) const {

}*/

Polygon& Polygon::operator=(const Polygon& p) {
    clear();

    this->layer = p.layer;
    this->points = p.points;
    this->pointsBC = p.pointsBC;
    this->projections = p.projections;
    this->axes = p.axes;
    this->centroid = p.centroid;
    this->lowerPoint = p.lowerPoint;
    this->subConvexPolygons = p.subConvexPolygons;
    this->numPolygons = p.numPolygons;
    this->BB = p.BB;

    return *this;
}

std::ostream& operator<<(std::ostream& os, const Polygon& p) {
    os << "Polygon (Layer: " << p.layer << "):\n\tVertices:\n";

    for (size_t i = 0; i < p.points.size(); ++i) {
        os << "\t\t(" << p[i].x << ", " << p[i].y << "), (" << p[(i + 1) % (int)p.points.size()].x <<
            ", " << p[(i + 1) % (int)p.points.size()].y << ")\n";
    }

    os << "\tSplines:\n";

    for (size_t i = 0; i < p.pointsBC.size(); ++i) {
        if (i > 0 && i % 4 == 0)
            os << "\n";
        os << "\t\t(" << p.getBCPoint(i).x << ", " << p.getBCPoint(i).y << ")\n";
    }

    os << "\n\tProjections:\n";

    for (const auto& pr : p.projections) {
        os << "\t\t\t(" << pr.min << ", " << pr.max << ")\n";
    }

    os << "\n\tAxes:\n";
    for (const auto& ax : p.axes) {
        os << "\t\t\t(" << ax.x << ", " << ax.y << ")\n";
    }

    os << "\tCentroid: (" << p.centroid.x << ", " << p.centroid.y << ")\n";

    if (p.getNumsubConvexPolygons() > 0) {
        os << "\tSubConvexxPolygons:\n";

        for (const auto& a : p.subConvexPolygons) {
            os << "\t\tConvexShape:\n";

            for (int i = 0; i < a.getNumPoints(); ++i) {
                os << "\t\t\t(" << a.points[i].x << ", " << a.points[i].y << "), (" << a.points[(i + 1) % a.getNumPoints()].x <<
                    ", " << a.points[(i + 1) % a.getNumPoints()].y << ")\n";
            }

            os << "\t\tProjections:\n";
            for (const auto& p : a.projections) {
                os << "\t\t\t(" << p.min << ", " << p.max << ")\n";
            }
        }
    }

    return os;
}

// Class Partition

Partition::VertexPartition::VertexPartition() {
    prev = NULL;
    next = NULL;
}

bool Partition::isConvex(const Point& a, const Point& b, const Point& c) {
    Point ab = b.minus(a);
    Point ac = c.minus(a);
    //Point ab = b - a;
    //Point ac = c - a;

    if (ab.crossSign(ac) >= 0) {
        return true;
    }

    return false;
}

bool Partition::isReflex(Point& a, Point& b, Point& c) {
    return !isConvex(a, b, c);
}

bool Partition::isInsideTriangle(Point& p1, Point& p2, Point& p3, Point& p) {
    if (isConvex(p1, p, p2)) {
        return false;
    }

    if (isConvex(p2, p, p3)) {
        return false;
    }

    if (isConvex(p3, p, p1)) {
        return false;
    }

    return true;
}

Point Partition::normalize(const Point& p) {
    Point res;
    double mod = sqrt(p.x * p.x + p.y * p.y);

    if (mod != 0) {
        return p / mod;
    }

    return Point(0, 0);
}

void Partition::updateVertex(VertexPartition* v, VertexPartition* vertices, const int& numVertices) {
    // Comprobamos si es concavo o convexo v
    v->convex = isConvex(v->prev->p, v->p, v->next->p);

    // Angulo entre los 2 vectores
    Point ab = v->prev->p.minus(v->p);
    Point ac = v->next->p.minus(v->p);
    //Point ab = v->prev->p - v->p;
    //Point ac = v->next->p - v->p;

    Point abn = normalize(ab);
    Point acn = normalize(ac);

    v->angle = abn.x * acn.x + abn.y * acn.y;

    //std::cout<<"PT: "<<v->p.x<<", "<<v->p.y<<" "<<v->angle<<std::endl;

    // Si es convexo es Ear, por tanto, comprobamos si en el triangulo que forman v, v->prev y v->next,
    // no existe ningun vertice restante de la figura
    if (v->convex) {
        v->isEar = true;

        for (int i = 0; i < numVertices; ++i) {
            // Si es uno de los 3 vertices, lo descartamos ya
            if ((v->p.x == vertices[i].p.x && v->p.y == vertices[i].p.y) ||
                (v->prev->p.x == vertices[i].p.x && v->prev->p.y == vertices[i].p.y) ||
                (v->next->p.x == vertices[i].p.x && v->next->p.y == vertices[i].p.y)) {
                continue;
            }

            // Comprobamos que no este dentro
            if (isInsideTriangle(v->prev->p, v->p, v->next->p, vertices[i].p)) {
                v->isEar = false;
                break;
            }
        }
    }
    else {
        v->isEar = false;
    }
}

bool Partition::triangulation(Polygon* pl, std::vector<Polygon>* result, const std::string& layer) {
    if (pl->getNumPoints() < 3) {
        return false;
    }

    if (pl->getNumPoints() == 3) {
        result->push_back(*pl);

        return true;
    }

    VertexPartition* vertices = new VertexPartition[pl->getNumPoints()];
    // Inicializamos los vertices del poligono: activo, coordenadas, 
    // vertice prev y next
    for (int i = 0; i < pl->getNumPoints(); ++i) {
        // Pasamos el vertice
        vertices[i].p = pl->getPoint(i);
        // Activamos vertice
        vertices[i].active = true;
        // Vertice siguiente al i
        vertices[i].next = &(vertices[(i + 1) % pl->getNumPoints()]);
        // Vertice a nterior al i
        vertices[i].prev = (i - 1 < 0) ? &(vertices[pl->getNumPoints() - 1]) : &(vertices[i - 1]);
    }

    // Inicializamos atributos es convexo, angulo e isEar
    for (int i = 0; i < pl->getNumPoints(); ++i) {
        updateVertex(&(vertices[i]), vertices, pl->getNumPoints());
    }

    // Comienza algoritmo Ear Clipping
    Polygon triangle;
    VertexPartition* vertexEar;
    bool ear;
    vertexEar = NULL;

    for (int i = 0; i < pl->getNumPoints() - 3; ++i) {
        ear = false;

        // Buscamos el vertice con un angulo mayor y que sea ear
        for (int j = 0; j < pl->getNumPoints(); ++j) {
            if (vertices[j].active && vertices[j].isEar) {
                if (!ear) {
                    ear = true;
                    vertexEar = &(vertices[j]);
                }
                else {
                    // Si ya teniamos un ear, comprobamos si el vertice actual tiene un mayor angulo
                    if (vertices[j].angle > vertexEar->angle) {
                        vertexEar = &(vertices[j]);
                    }
                }
            }
        }

        //std::cout<<vertexEar->p.x<<", "<<vertexEar->p.y<<" "<<vertexEar->angle<<std::endl;

        // Si no hemos encontrado vertices ear, ya no podemos realizar una triangulacion a la pieza
        if (!ear) {
            delete[] vertices;

            return false;
        }

        // Anaydimos el triangulo al resultado
        result->push_back(Polygon(3, layer));
        result->back()[0] = vertexEar->prev->p;
        result->back()[1] = vertexEar->p;
        result->back()[2] = vertexEar->next->p;

        // Actualizamos los tres vertices
        vertexEar->active = false;
        vertexEar->prev->next = vertexEar->next;
        vertexEar->next->prev = vertexEar->prev;

        if (i == pl->getNumPoints() - 4) {
            break;
        }

        updateVertex(vertexEar->prev, vertices, pl->getNumPoints());
        updateVertex(vertexEar->next, vertices, pl->getNumPoints());
    }

    // Quedaran 3 vertices (1 triangulo) con solo uno de estos vertices activo
    for (int i = 0; i < pl->getNumPoints(); ++i) {
        if (vertices[i].active) {
            result->push_back(Polygon(3, layer));
            result->back()[0] = vertices[i].prev->p;
            result->back()[1] = vertices[i].p;
            result->back()[2] = vertices[i].next->p;

            break;
        }
    }

    delete[] vertices;


    return true;
}

bool Partition::convexPartition(Polygon* pl, std::vector<Polygon>* result, const std::string& layer) {
    // Poner funcion isConcavePolygon para comprobar si el poligono es concavo
    if (pl->getNumPoints() < 3) {
        return false;
    }

    // Triangulamos el poligono
    std::vector<Polygon> triangles;

    if (!triangulation(pl, &triangles, layer)) {
        return false;
    }

    // Algoritmo de Hertel-Mehlhorn
    std::vector<Polygon>::iterator iter1, iter2;
    Polygon* pl1, * pl2;
    Point d1, d2;
    int v11 = -1, v12 = -1, v21 = -1, v22 = -1;
    bool diagonal = false;

    pl1 = NULL;
    pl2 = NULL;

    for (iter1 = triangles.begin(); iter1 != triangles.end(); ++iter1) {
        pl1 = &(*iter1);

        for (v11 = 0; v11 < pl1->getNumPoints(); ++v11) {
            d1 = pl1->getPoint(v11);
            v12 = (v11 + 1) % pl1->getNumPoints();
            d2 = pl1->getPoint(v12);
            //std::cout<<"=== "<<d1.x<<", "<<d1.y<<" "<<d2.x<<", "<<d2.y<<std::endl;

            // Comprobamos si la arista (d1, d2) de la pieza p1, hace de diagonal con otra pieza p2
            diagonal = false;

            for (iter2 = iter1; iter2 != triangles.end(); ++iter2) {
                if (iter1 == iter2)
                    continue;

                pl2 = &(*iter2);

                // Buscamos en la pieza pl2 si tiene la arista (d1, d2)
                for (v21 = 0; v21 < pl2->getNumPoints(); ++v21) {
                    if ((pl2->getPoint(v21).x != d2.x) || (pl2->getPoint(v21).y != d2.y))
                        continue;

                    v22 = (v21 + 1) % pl2->getNumPoints();

                    if ((pl2->getPoint(v22).x != d1.x) || (pl2->getPoint(v22).y != d1.y))
                        continue;

                    // La arista (v12, v22) == (d1, d2)
                    diagonal = true;

                    break;
                }

                if (diagonal)
                    break;
            }

            // Si la pieza pl2 no tiene una de sus aristas igual a (d1, d2), pasamos a la siguiente pieza (si hay)
            if (!diagonal)
                continue;

            //std::cout<<"SI"<<std::endl;

            // Al ser (v21, v22) una arista de pl1 y pl2, comprobamos si podemos eliminarla y asi fusionar pl1 y pl2
            // Comprobamos si el punto d1 es convexo
            Point p = pl1->getPoint(v11);
            int pos = (v11 - 1 < 0) ? pl1->getNumPoints() - 1 : v11 - 1;
            Point p1 = pl1->getPoint(pos);
            pos = (v22 + 1) % pl2->getNumPoints();
            Point p2 = pl2->getPoint(pos);

            //std::cout<<"-> "<<p1.x<<", "<<p1.y<<" "<<p.x<<", "<<p.y<<" "<<p2.x<<", "<<p2.y<<std::endl;

            if (!isConvex(p1, p, p2)) {
                //std::cout<<"COONCAVO"<<std::endl;
                continue;
            }

            //std::cout<<"CONVEXO"<<std::endl;

            // Comprobamos si el punto d2 es convexo
            p = pl1->getPoint(v12);
            pos = (v12 + 1) % pl1->getNumPoints();
            p2 = pl1->getPoint(pos);
            pos = (v21 - 1 < 0) ? pl2->getNumPoints() - 1 : v21 - 1;
            p1 = pl2->getPoint(pos);

            //std::cout<<"-> "<<p1.x<<", "<<p1.y<<" "<<p.x<<", "<<p.y<<" "<<p2.x<<", "<<p2.y<<std::endl;

            if (!isConvex(p1, p, p2)) {
                //std::cout<<"COONCAVO"<<std::endl;
                continue;
            }

            //std::cout<<"CONVEXO"<<std::endl;

            //std::cout<<"-------------------------------------"<<std::endl;

            // Si ambos son convexos, creamos un nuevo poligono, que sustituira al pl1, y pl2 se eliminara
            Polygon newPolygon(pl1->getNumPoints() + pl2->getNumPoints() - 2, layer);
            int j = 0;

            for (int i = v12; i != v11; i = (i + 1) % pl1->getNumPoints()) {
                newPolygon[j] = pl1->getPoint(i);
                ++j;
            }

            for (int i = v22; i != v21; i = (i + 1) % pl2->getNumPoints()) {
                newPolygon[j] = pl2->getPoint(i);
                ++j;
            }

            // Eliminamos pl2
            triangles.erase(iter2);
            // Sustituimos pl1 por newPolygon
            *iter1 = newPolygon;
            pl1 = &(*iter1);
            v11 = -1;

            continue;
        }
    }

    /*std::cout<<"Resultado ConvexPartition:\n";
    for(auto i = triangles.begin(); i != triangles.end(); ++i){
        std::cout<<*i<<"\n";
    }*/


    // Pasamos el triangles a result
    for (auto p : triangles) {
        result->push_back(p);
    }

    return true;
}

void Partition::concaveToConvex(Polygon* p) {
    // Triangulamos el poligono
    std::vector<Polygon> result;

    triangulation(p, &result, p->getLayer());

    //for(size_t i = 0; i < result.size(); ++i){
    //    std::cout<<result[i]<<std::endl;
    //}

    result.clear();

    convexPartition(p, &result, p->getLayer());

    // Anyadimos las piezas convexas que sustituiran a la original para comprobar las colisiones
    for (size_t i = 0; i < result.size(); ++i) {
        p->addConvexPolygon(&result[i]);
    }
}

// Clase dxfFilter

dxfFilter::dxfFilter(const std::string& fileName) {
    file.open(fileName, std::ifstream::binary);
    std::string line;

    while (getline(file, line))
    {
        if(!line.empty())
            line.pop_back();

        if (line == "ENTITIES") {
#ifdef PRINTDEBUG
            std::cout << "FINNN" << std::endl;
#endif
            break;
        }
        //std::cout<<line<<std::endl;
    }
#ifdef PRINTDEBUG
    std::cout << "????????_" << std::endl;
#endif
}

dxfFilter::~dxfFilter() {
    file.close();
}

std::string dxfFilter::readLayer() {
    std::string line;

    while (getline(file, line)) {
        if (!line.empty())
            line.pop_back();

        if (line == "  8") {
            getline(file, line);
            //std::cout << "-> " << line << std::endl;

            if (!line.empty())
                line.pop_back();

            return line;
        }
    }

    return "";
}

void dxfFilter::addControlPoint(const DL_ControlPointData& data) {
    // Almacenamos los vértices de control (4) que conforman cada curva cúbica de Bezier 
    std::string layer = "";

    if (polys.size() == 0 || polys.back()->getNumPointsBC() % 4 == 0)
        layer = readLayer();
    else
        layer = polys.back()->getLayer();

    if (polys.size() == 0) {
#ifdef PRINTDEBUG
        std::cout << "NOOOOOOOO" << std::endl;
#endif
        polys.pushBack(Polygon(layer));
    }
    else {
        if (polys.back()->getLayer() != layer) {
            // Nuevo polígono
            // Calculamos centroide, proyecciones, ... del polígono anterior (ya completo)
            polys.back()->calculateCentroid();

            printPolygons();

            for (int i = 0; i < polys.back()->getNumPoints(); ++i) {
                polys.back()->setProjection(polys.back()->project(polys.back()->getAxisComponent(i)));
            }

            if (polys.back()->isConcavePolygon()) {
#ifdef PRINTDEBUG
                std::cout << "ES CONCAVO" << std::endl;
#endif
                Partition pp;
                pp.concaveToConvex(polys.back());
            }

            // Nuevo polígono
            polys.pushBack(Polygon(layer));
        }
    }

    // Ver quer vertices de la curva usamos para el convexhull

    // Insertamos el punto de control
    Point p(data.x, data.y);

    // Si el elemento anterior del polígono fue una arista (núm. de vértices del vector points es impar), insertamos 
    // el vertice p en points
    if (polys.back()->getNumPointsBC() == 0 || (polys.back()->getNumPointsBC() + 1) % 4 == 0) {
        Point p1;

        if (polys.back()->getNumPointsBC() == 0) {
            p1 = p;
        }
        else {
            p1 = polys.back()->getBCPoint(polys.back()->getNumPointsBC() - 3);
            double vX = p.x - p1.x;
            double vY = p.y - p1.y;

            polys.back()->setAxis(Point(-vY, vX));
        }

        //polys.back()->updateLowerVertice(p1);

        if (polys.back()->getNumPoints() == 0 || polys.back()->getPoint(polys.back()->getNumPoints() - 1) != p1) {
            polys.back()->insertVertice(p1);
            polys.back()->updateLowerVertice(p1);
        }
    }

    polys.back()->insertVerticeBC(p);
#ifdef PRINTDEBUG
    std::cout << "---> " << polys.back()->getNumPointsBC() << std::endl;
    
    std::cout <<layer<< " -> " << data.x << " " << data.y << " " << data.z << std::endl;
#endif
}

void dxfFilter::addLine(const DL_LineData& data) {
    std::string layer = readLayer();
    // Metemos los vertices que conforman la pieza en un vector

    // Si aun no hay ninguna pieza leida, la creamos
    if (polys.size() == 0) {
#ifdef PRINTDEBUG


        std::cout << "SIIIII" << std::endl;
#endif
        polys.pushBack(Polygon(layer));
        // polys.push_back(new Polygon(layer));
    }
    else {
        if (polys.back()->getLayer() != layer) {
            // Obtenemos el centroide del poligono
            polys.back()->calculateCentroid();

            // Obtenemos la proyeccion del poligono sobre el eje de proyeccion correspondiente a la ultima arista 
            // insertada
            for (int i = 0; i < polys.back()->getNumPoints(); ++i) {
                polys.back()->setProjection(polys.back()->project(polys.back()->getAxisComponent(i)));
            }

            // Comprobamos si la pieza concava, y si lo es, la descomponemos en piezas convexas
            if (polys.back()->isConcavePolygon()) {
#ifdef PRINTDEBUG
                std::cout<<"ES CONCAVO"<<std::endl;
#endif
                Partition pp;
                pp.concaveToConvex(polys.back());
            }
            //else{
            //    std::cout<<"ES CONVEXO"<<std::endl;
            //}

            // Obtenemos el número de polígonos convexos que forman la última pieza (si es que es concava) + en propio polígono
            polys.back()->setNumPolygons(polys.back()->getNumsubConvexPolygons() + 1);

            // Insertamos nueva pieza
            //polys.push_back(new Polygon(layer));
            polys.pushBack(Polygon(layer));
        }
    }

    // Introducimos el componente en la ultima pieza insertada
    Point p1(data.x1, data.y1);
    Point p2(data.x2, data.y2);

    // Obtenemos el eje de proyeccion correspondiente
    // Obtenemos el vector de los 2 puntos que forman la arista
    double vX = p2.x - p1.x;
    double vY = p2.y - p1.y;
    //std::cout<<vX<<" "<<vY<<std::endl;
    // Obtenemos el vector perpendicular a este (= eje de proyeccion)
    // NO estan normalizados
#ifdef PRINTDEBUG
    std::cout <<"----> " << Point(-vY, vX).x << " " << Point(-vY, vX).y << std::endl;
#endif
    polys.back()->setAxis(Point(-vY, vX));

    // Insertamos el vertice en el vector vertices
    Point newPoint(data.x1, data.y1);
#ifdef PRINTDEBUG
    std::cout << data.x1 << " * " << data.y1 << std::endl;
#endif

    // Si el vertice esta "mas abajo" que el vertice actual (lowerVertice), lo actualizamos (esto es para el
    // convexHull)
    polys.back()->updateLowerVertice(newPoint);

    polys.back()->insertVertice(newPoint);

    // Actualizamos el Bounding Box
    polys.back()->updateBB(newPoint);
}

void dxfFilter::readPolygons(const std::string& nameFile, DL_Dxf* dxf) {
    if (!dxf->in(nameFile, this)) {
        std::cerr << "NO se ha podido leer el archivo" << nameFile << std::endl;
        exit(1);
    }
#ifdef PRINTDEBUG
    std::cout << "??????" << std::endl;
    std::cout << polys.size() << std::endl;
#endif
        Polygon p = *polys.back();
#ifdef PRINTDEBUG
    for (int i = 0; i < p.getNumPoints(); ++i) {
        std::cout << p[i].x<<" "<<p[i].y << std::endl;
    }
#endif

    if (this->polys.size() == 0)
        return;

    // Hay que obtener la proyeccion de la ultima pieza leida
    for (int i = 0; i < this->getPolygon(this->getNumPolygons() - 1)->getNumPoints(); ++i) {
        this->getPolygon(this->getNumPolygons() - 1)->setProjection(
            this->getPolygon(this->getNumPolygons() - 1)->project(
                this->getPolygon(this->getNumPolygons() - 1)->getAxis(i)
            )
        );
    }

    // Obtenemos el centroide de la ultima pieza leida
    this->getPolygon(this->getNumPolygons() - 1)->calculateCentroid();

    if (this->getPolygon(this->getNumPolygons() - 1)->isConcavePolygon()) {
		Partition pp;
        pp.concaveToConvex(this->getPolygon(this->getNumPolygons() - 1));
    }

    // Obtenemos el número de polígonos convexos que forman la última pieza (si es que es concava) + en propio polígono
    polys.back()->setNumPolygons(polys.back()->getNumsubConvexPolygons() + 1);

    // Hacemos la copia del vector de polígonos leídos aquí
    polysCP = polys;
}

bool dxfFilter::collisionDetection1Polygon(void*** polygons, const int& numPolygons, const int& indexPolygon) {
    if (numPolygons == 1)
        return false;

    for (int i = 0; i < numPolygons; ++i) {
        if (i != indexPolygon) {
            int numPols1 = *(int*)(polygons[i][0]);
            int numPols2 = *(int*)(polygons[indexPolygon][0]);

            if (numPols1 == 1 && numPols2 == 1) {
                // Ambos son polígonos convexos
                if (checkCollision2(polygons[i], 0, polygons[indexPolygon], 0) && checkCollision2(polygons[indexPolygon], 0, polygons[i], 0)) {
                    return true;
                }
            }
            else if (numPols1 > 1 && numPols2 == 1) {
                for (int k = 0; k < numPols1; ++k) {
                    if (checkCollision2(polygons[i], 4 * k, polygons[indexPolygon], 0) && checkCollision2(polygons[indexPolygon], 0, polygons[i], 4 * k)) {
                        return true;
                    }
                }
            }
            else if (numPols1 == 1 && numPols2 > 1) {
                for (int k = 0; k < numPols1; ++k) {
                    if (checkCollision2(polygons[i], 0, polygons[indexPolygon], 4 * k) && checkCollision2(polygons[indexPolygon], 4 * k, polygons[i], 0)) {
                        return true;
                    }
                }
            }
            else {
                for (int k = 0; k < numPols1; ++k) {
                    for (int l = 0; l < numPols2; ++l) {
                        if (checkCollision2(polygons[i], 4 * k, polygons[indexPolygon], 4 * l) && checkCollision2(polygons[indexPolygon], 4 * l, polygons[i], 4 * k)) {
                            return true;
                        }
                    }
                }
            }
        }
    }

    return false;
}

bool dxfFilter::collisionDetectionBB(void **p1, void **p2) {
    // xmax ymax xmin ymin
    double* BB1 = (double*)(p1[3]);
    double* BB2 = (double*)(p2[3]);

    //std::cout << "(" << BB1[0] << " >= " << BB2[0] << " && " << BB1[2] << " <= " << BB2[0] << ") || (" << BB1[2] << " <= " << BB2[2] << " && " << BB1[0] <<
    //    ">= " << BB2[2] << ") || (" << BB1[0] << " <= " << BB2[0] << " && " << BB1[2] << " >= " << BB2[2] << ")\n";
    if ((BB1[0] >= BB2[0] && BB1[2] <= BB2[0]) || (BB1[2] <= BB2[2] && BB1[0] >= BB2[2]) || (BB1[0] <= BB2[0] && BB1[2] >= BB2[2])) {
        //std::cout << "(" << BB1[1] << " >= " << BB2[1] << " && " << BB1[3] << " <= " << BB2[1] << ") || (" << BB1[3] << " <= " << BB2[3] << " && " << BB1[1] <<
        //    " >= " << BB2[3] << ") || (" << BB1[1] << " <= " << BB2[1] << " && " << BB1[3] << " >= " << BB2[3] << ")\n";
        if ((BB1[1] >= BB2[1] && BB1[3] <= BB2[1]) || (BB1[3] <= BB2[3] && BB1[1] >= BB2[3]) ||(BB1[1] <= BB2[1] && BB1[3] >= BB2[3])) {
            return true;
        }
    }

    return false;
}

void printPolygon(void** polygon) {
    // Num pol + centroide + lowerPoint + BB + Num convex poly * (Num pts + puntos + axes + proyecciones)
    int numPol = *(int*)(polygon[0]);
    std::cout << "\nPoligono:\n";
    Point* cent = (Point*)(polygon[1]);
    std::cout << "Centroide: " << cent->x << " " << cent->y << std::endl;
    int index = 4;
    for (int i = 0; i < numPol; ++i) {
        int numPuntos = *(int*)(polygon[index++]);
        std::cout << "Numero puntos: " << numPuntos << "\nPuntos:\n";
        Point* points = (Point*)(polygon[index++]);
        for (int j = 0; j < numPuntos; ++j) {
            std::cout << points[j].x << ", " << points[j].y << std::endl;
        }
        std::cout << "Axes:\n";
        Point* axes = (Point*)(polygon[index++]);
        for (int j = 0; j < numPuntos; ++j) {
            std::cout << axes[j].x << ", " << axes[j].y << std::endl;
        }
        std::cout << "Projections:\n";
        Projection* proj = (Projection*)(polygon[index++]);
        for (int j = 0; j < numPuntos; ++j) {
            std::cout << proj[j].min << ", " << proj[j].max << std::endl;
        }
    }
    std::cout << "\n\n";
}

int dxfFilter::collisionDetection2(void*** polygons, const int &numPolygons) {
	int col = 0;
	if (numPolygons == 1)
		return col;

	for (int i = 0; i < numPolygons; ++i) {
		int j = i + 1;
		while (j < numPolygons) {
            if (collisionDetectionBB(polygons[i], polygons[j])) {
#ifdef PRINTDEBUG
                std::cout << "BB:\n";
                printPolygon(polygons[i]);
                printPolygon(polygons[j]);
                std::cout << std::endl;
#endif
                int numPols1 = *(int*)(polygons[i][0]);
                int numPols2 = *(int*)(polygons[j][0]);
                if (numPols1 == 1 && numPols2 == 1) {
                    // Ambos son polígonos convexos
                    if (checkCollision2(polygons[i], 0, polygons[j], 0) && checkCollision2(polygons[j], 0, polygons[i], 0)) {
                        col++;
                        //return true;
                    }
                }
                else if (numPols1 > 1 && numPols2 == 1) {
                    for (int k = 1; k < numPols1; ++k) {
                        if (checkCollision2(polygons[i], 4 * k, polygons[j], 0) && checkCollision2(polygons[j], 0, polygons[i], 4 * k)) {
                            col++;
                            //return true;
                        }
                    }
                }
                else if (numPols1 == 1 && numPols2 > 1) {
                    for (int k = 1; k < numPols2; ++k) {
                        if (checkCollision2(polygons[i], 0, polygons[j], 4 * k) && checkCollision2(polygons[j], 4 * k, polygons[i], 0)) {
                            col++;
                            //return true;
                        }
                    }
                }
                else {
                    for (int k = 1; k < numPols1; ++k) {
                        for (int l = 1; l < numPols2; ++l) {
                            if (checkCollision2(polygons[i], 4 * k, polygons[j], 4 * l) && checkCollision2(polygons[j], 4 * l, polygons[i], 4 * k)) {
                                col++;
                                //return true;
                            }
                        }
                    }
                }
            }
            
			++j;
		}
	}

	return col;
}

bool dxfFilter::collisionDetection() const {
    // Por ahora, probamos si existe interferencia entre todas las combinaciones de pares de piezas posibles
    // hasta que no queden o haya colision en uno de ellos
    bool collision = false;

    if (polys.size() <= 1)
        return false;

    for (size_t i = 0; i < polys.size(); ++i) {
        size_t j = i + 1;

        while (j < polys.size() && !collision) {
            if (polys.get(i)->getNumsubConvexPolygons() == 0 && polys.get(j)->getNumsubConvexPolygons() == 0) {
                if (polys.get(i)->checkCollision(polys.get(j)) && polys.get(j)->checkCollision(polys.get(i))) {
                    if (polys.get(i)->getNumPointsBC() > 0 && !polys.get(i)->checkCollisionBC(polys.get(j))) {
                        return false;
                    }
                    // Hay colision, por tanto, paramos
                    //std::cout << "COLISION" << std::endl;
                    return true;
                }
                else {
                    if (polys.get(i)->getNumPointsBC() > 0 && polys.get(i)->checkCollisionBC(polys.get(j)) ||
                        polys.get(j)->getNumPointsBC() > 0 && polys.get(j)->checkCollisionBC(polys.get(i))) {
                        return true;
                    }
                }
            }
            else if (polys.get(i)->getNumsubConvexPolygons() > 0 && polys.get(j)->getNumsubConvexPolygons() == 0) {
                for (int k = 0; k < polys.get(i)->getNumsubConvexPolygons(); ++k) {
                    if (polys.get(i)->getConvexShape(k).checkCollision(polys.get(j)) &&
                        polys.get(j)->checkCollision(&(polys.get(i)->getConvexShape(k)))) {
                        //std::cout << "COLISIOON" << std::endl;
                        return true;
                    }
                    else {
                        if (polys.get(i)->getNumPointsBC() > 0 && polys.get(i)->checkCollisionBC(polys.get(j)) ||
                            polys.get(j)->getNumPointsBC() > 0 && polys.get(j)->checkCollisionBC(polys.get(i))) {
                            return true;
                        }
                    }
                }
            }
            else if (polys.get(i)->getNumsubConvexPolygons() == 0 && polys.get(j)->getNumsubConvexPolygons() > 0) {
                for (int k = 0; k < polys.get(j)->getNumsubConvexPolygons(); ++k) {
                    if (polys.get(j)->getConvexShape(k).checkCollision(polys.get(i)) &&
                        polys.get(i)->checkCollision(&(polys.get(j)->getConvexShape(k)))) {
                        //std::cout << "COLISIOOoN" << std::endl;
                        return true;
                    }
                    else {
                        if (polys.get(i)->getNumPointsBC() > 0 && polys.get(i)->checkCollisionBC(polys.get(j)) ||
                            polys.get(j)->getNumPointsBC() > 0 && polys.get(j)->checkCollisionBC(polys.get(i))) {
                            return true;
                        }
                    }
                }
            }
            else {
                for (int k = 0; k < polys.get(i)->getNumsubConvexPolygons(); ++k) {
                    for (int l = 0; l < polys.get(j)->getNumsubConvexPolygons(); ++l) {
                        if (polys.get(i)->getConvexShape(k).checkCollision(&(polys.get(j)->getConvexShape(l))) &&
                            polys.get(j)->getConvexShape(l).checkCollision(&(polys.get(i)->getConvexShape(k)))) {
                            //std::cout << "COLISIOOooN" << std::endl;
                            return true;
                        }
                        else {
                            if (polys.get(i)->getNumPointsBC() > 0 && polys.get(i)->checkCollisionBC(polys.get(k)) ||
                                polys.get(k)->getNumPointsBC() > 0 && polys.get(k)->checkCollisionBC(polys.get(i))) {
                                return true;
                            }
                        }
                    }
                }
            }

            ++j;
        }
    }

    return false;
}

int orientation(Point p, Point q, Point r) {
    // Vemos que orientacion tienen los vertices p, q y r: 0 -> colineal, 1 -> sentido horario, 2 -> antihorario
    double val = (q.y - p.y) * (r.x - q.x) -
        (q.x - p.x) * (r.y - q.y);

    if (val == 0)
        return 0;

    return /*(val > 0) ? 1 : 2*/(val > 0) ? 2 : 1;
}

double distSq(Point p1, Point p2) {
    // Distancia entre p1 y p2
    return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
}

Point p0;

int compare(const void* vp1, const void* vp2) {
    // Esta funcion la usa std::qsort(...) para ordenar el array de puntos respecto a p0
    Point* p1 = (Point*)vp1;
    Point* p2 = (Point*)vp2;

    // Orientacion 
    int o = orientation(p0, *p1, *p2);

    if (o == 0) {
        // p0, p1 y p2 son colineales, por tanto, nos quedamos con el punto mas alejado de p0
        return (distSq(p0, *p2) >= distSq(p0, *p1)) ? -1 : 1;
    }

    return (o == 2) ? -1 : 1;
}

Point nextToTop(std::stack<Point>& s) {
    Point aux = s.top();
    s.pop();
    Point res = s.top();
    s.push(aux);

    return res;
}

double dxfFilter::convexHullArea(std::stack<Point> s) {
    double area = 0.0;
    int n = s.size();
    Point* points;
    points = new Point[n];
    // Point points[n];
    int i = 0;
#ifdef PRINTDEBUG
    std::cout<<"CONVEXHULL:\n";
#endif
    while (!s.empty()) {
        points[i] = s.top();
#ifdef PRINTDEBUG
		std::cout << points[i].x << ",, " << points[i].y << std::endl;
#endif
        // Por ahora metemos todos los putnos en el array points (como si todo fueran aristas)
        convexHullPolygon.insertVertice(s.top());
        //addVerticeConvexHull(s.top());
        //result.push_back(s.top());
        //std::cout<<s.top().x<<", "<<s.top().y<<" ";
        s.pop();

        ++i;
    }
    //std::cout<<std::endl;

    int j = n - 1;
    for (int i = 0; i < n; ++i) {
        area += (points[j].x + points[i].x) * (points[j].y - points[i].y);
        j = i;
    }

    delete points;

    return std::abs(area / 2.0);
}

Point dxfFilter::calculateLowerPoint() const {
    Point p = polys.front()->getLowerPoint();

    for (size_t i = 1; i < polys.size(); ++i) {
        if (polys.get(i)->getLowerPoint() < p)
            p = polys.get(i)->getLowerPoint();
    }

    //std::cout<<"LOwerPoint: "<<p.x<<", "<<p.y<<std::endl;

    return p;
}

double dxfFilter::convexHull() {
    //if (verticesConvexHull.size() > 0)
        //verticesConvexHull.clear();

    // Obtenemos el array con todos los puntos que forman el conjunto de piezas
    std::vector<Point> points;
    int lw = 0;

    if (polys.size() == 0)
        return 0;

    //for(auto p : polys){
    for (size_t p = 0; p < polys.size(); ++p) {
        for (int i = 0; i < polys.get(p)->getNumPoints(); ++i) {
            if (polys.get(p)->getPoint(i) == calculateLowerPoint())
                lw = (int)points.size();

            points.push_back(polys.get(p)->getPoint(i));
        }
    }
#ifdef PRINTDEBUG
    std::cout << "ConvexHUll points:\n";
    for (int i = 0; i < points.size(); ++i)
        std::cout << points[i].x << ", " << points[i].y << std::endl;

    std::cout << "lw " << points[lw].x<<", "<<points[lw].y << std::endl;
#endif
    // Ponemos el punto mas bajo en la primera posicion del array points
    std::swap(points[0], points[lw]);

    // Ordenamos el resto de puntos (numPoints - 1) por angulo polar en orden antihorario 
    // alrededor de points[0]. Es decir, un punto p1 ira antes que p2 en el array si p2 tiene un angulo 
    // polar mayor (en sentido contrario a las agujas del reloj) que p1
    // NOTA: la funcion de comparacion (compare) debe devolver un valor entero negativo 
    // si el primer argumento es menor que el segundo, un valor entero positivo si el primer argumento es mayor que 
    // el segundo y cero si los argumentos son equivalentes
    p0 = points[0];
    std::qsort(&points[1], points.size() - 1, sizeof(Point), compare);

    // Si dos o más puntos forman el mismo angulo con p0 (colineales), eliminamos todos estos, menos el que esta 
    // mas alejado de p0. En la clasificación anterior, nuestro criterio fue mantener el punto más alejado al final 
    // cuando mas de un punto tiene el mismo ángulo
    int m = 1;

    for (size_t i = 1; i < points.size(); ++i) {
        // Si hay 2 o mas puntos tienen el mismo angulo con p0, removemos todos menos el mas lejano (i + 1)
        while (i < points.size() - 1 && orientation(p0, points[i], points[i + 1]) == 0) {
            ++i;
        }

        points[m] = points[i];
        ++m;
    }

    // Si m < 3, no se puede obtener el convexHull
    if (m < 3)
        return false;

    // Creamos la pila y colocamos los 3 primeros puntos
    std::stack<Point> s;
    s.push(points[0]);
    s.push(points[1]);
    s.push(points[2]);

    for (int i = 3; i < m; ++i) {
        // Eliminamos de la pila el punto top de la pila, que el angulo formado con los puntos next-top (punto 
        // anterior al top de la pila), top point y point[i], NO sea en sentido antihorario
#ifdef PRINTDEBUG
        std::cout << s.top().x << ", " << s.top().y << std::endl;
#endif
        while (s.size() > 1 && orientation(nextToTop(s), s.top(), points[i]) != 2) {
            s.pop();
        }

        s.push(points[i]);
    }

    return convexHullArea(s);
}


void dxfFilter::freeData(void***& polygons) const {

	// Num pol + centroide + lowerPoint + BB + Num convex poly * (Num pts + puntos + axes + proyecciones)
	for (int i = 0; i < (int)polys.size(); i++) {
		for (int j = 0; j < *(int*)(polygons[i][0]); j++) {
			delete[](Point*)(polygons[i][4 * j + 5]);
			(polygons[i][4 * j + 5]) = NULL;
			delete[](Point*)(polygons[i][4 * j + 6]);
			polygons[i][4 * j + 6] = NULL;
			delete[] polygons[i][4 * j + 7];
			polygons[i][4 * j + 7] = NULL;
		}
		delete polygons[i][4 * 0 + 4];
		polygons[i][4 * 0 + 4] = NULL;

		delete polygons[i];
		polygons[i] = NULL;
	}

	delete polygons;
	polygons = NULL;
}

void dxfFilter::transferData(void*** &polygons) const {
	// Num pol + centroide + lowerPoint + BB + Num convex poly * (Num pts + puntos + axes + proyecciones)
	polygons = new void**[(int)polys.size()];
	int* numPols = new int[(int)polys.size()];
	Point* points;
	int* numPoints;
	Point* axesPoints;
	Projection* projections;
	Point* cent = new Point[(int)polys.size()];
	Point* lw = new Point[(int)polys.size()];

	int indexPolygon;
	int pols;
	for (int i = 0; i < (int)polys.size(); ++i) {
		indexPolygon = 0;
		pols = (polys.get(i)->getNumsubConvexPolygons() == 0) ? 1 : polys.get(i)->getNumsubConvexPolygons() + 1;
		polygons[i] = new void*[1 + 1 + 1 + 1 + pols * 4];
		numPoints = new int[pols];

		// Num pol
		numPols[i] = 1 + polys.get(i)->getNumsubConvexPolygons();
		polygons[i][indexPolygon++] = &(numPols[i]);
		// centroide
		cent[i] = polys.get(i)->getCentroid();
		polygons[i][indexPolygon++] = &(cent[i]);//polys.get(i)->getCentroid2();
		// lowerPoint
		lw[i] = polys.get(i)->getLowerPoint();
		polygons[i][indexPolygon++] = &(lw[i]);//polys.get(i)->getLowerPoint2();
        // Bounding Box
        polygons[i][indexPolygon++] = polys.get(i)->getBB();
#ifdef PRINTDEBUG
        for (int k = 0; k < 4; k++) {
            std::cout << ((double*)(polygons[i][indexPolygon - 1]))[k] << std::endl;
        }
#endif
		numPoints[0] = polys.get(i)->getNumPoints();
		points = new Point[numPoints[0]];
		for (int j = 0; j < polys.get(i)->getNumPoints(); ++j) {
			points[j] = polys.get(i)->getPoint(j);
		}
		axesPoints = new Point[numPoints[0]];
		for (int j = 0; j < polys.get(i)->getNumPoints(); ++j) {
			axesPoints[j] = polys.get(i)->getAxis(j);
		}
		projections = new Projection[numPoints[0]];
		for (int j = 0; j < polys.get(i)->getNumPoints(); ++j) {
			projections[j] = polys.get(i)->getProjection(j);
		}

		int j = 0;
		while (true) {
			// Num pts
			polygons[i][indexPolygon++] = &(numPoints[j]);
			// puntos
			polygons[i][indexPolygon++] = points;
			// axes
			polygons[i][indexPolygon++] = axesPoints;
			// projections
			polygons[i][indexPolygon++] = projections;

			if (++j >= pols)
				break;

			if (pols > 1) {
				numPoints[j] = polys.get(i)->getConvexShape2(j - 1)->getNumPoints();
				points = new Point[numPoints[j]];
				for (int k = 0; k < polys.get(i)->getConvexShape2(j - 1)->getNumPoints(); ++k) {
					points[k] = polys.get(i)->getConvexShape2(j - 1)->getPoint(k);
				}
				axesPoints = new Point[numPoints[j]];
				for (int k = 0; k < polys.get(i)->getConvexShape2(j - 1)->getNumPoints(); ++k) {
					axesPoints[k] = polys.get(i)->getConvexShape2(j - 1)->getAxis(k);
				}
				projections = new Projection[numPoints[j]];
				for (int k = 0; k < polys.get(i)->getConvexShape2(j - 1)->getNumPoints(); ++k) {
					projections[k] = polys.get(i)->getConvexShape2(j - 1)->getProjection(k);
				}
			}
		}
	}
}

/*__global__ double dxfFilter::evFunction2(void*** polygons, const int& thread, const double& x, const double& y, const double& r) {
    extern __shared__ void*** buffer;



    //for (int i = 0; i < numPolygons; ++i)
    printPolygon(polygons[thread]);
    // Movemos polígono
    move2(polygons[thread], Point(x, y));
    // Rotamos polígono
    rotate2(polygons[thread], r);

    // Copiar resultado en estructura shared (estructura similar al void***, pero compartida entre GPU y CPU)
    

    // Sincronizamos los hilos
    __syncthreads();



    // Detectar posibles colisiones y si no se encuentran, obtener convex hull
    if (collisionDetection2(polygons, numPolygons)) {
        std::cout << "Colision" << std::endl;
        cx = -1.0;
    }
    else {
        std::cout << "NO colision" << std::endl;
        cx = convexHull();
    }

    // Obtener dxf salida para visualizar la distribucion de las piezas y convex hull
    if (dxf)
        drawFile(dxf, cx, nameFile, area);

    // Restauramos los poligonos al estado inicial
    polys = polysCP;
    std::cout << std::endl;
    printPolygons();

    return cx;
}*/

double dxfFilter::evFunction(DL_Dxf* dxf, const std::string& nameFile, Area* area, const Point& p, const double& angle) {
    double cx = 0;
    // Las operatiocnes se realizan sobre el vector de copia, por tanto, lo limpiamos e inicializamos al original
    //polysCP = polys;

    printPolygons();

    //for (int i = 0; i < (int)polys.size(); ++i) {
        //polys.get(0)->move(p);
        //polys.get(0)->rotate(angle);
    //}

    polys.get(0)->move(p);

    printPolygons();

    /*if (collisionDetection())
        cx = -1.0;
    else
        cx = convexHull();

    if (dxf)
        drawFile(dxf, cx, nameFile, area);*/

    return cx;

    /*for (int i = 0; i < (int)polysCP.size(); ++i) {
        polysCP.get(i)->move(p);
        polysCP.get(i)->rotate(angle);
    }

    if (collisionDetection()) {
        cx = -1.0;
        if (dxf)
            drawFile(dxf, cx, nameFile, area);
    }
    else {
        cx = convexHull();
        if (area && !area->isInsideArea(verticesConvexHull)) {
            cx = -1.0;
        }

        if (dxf)
            drawFile(dxf, cx, nameFile, area);
    }

    return cx;*/
}

void dxfFilter::drawFile(DL_Dxf* dxf, const double& convexHullArea, const std::string& nameFile, Area* area) const {
    DL_Codes::version exportVersion = DL_Codes::AC1009;
    DL_WriterA* dw = dxf->out(nameFile.c_str(), exportVersion);

    if (polys.size() == 0)
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

    int numLayers = (int)polys.size();
    dw->tableLayers(numLayers);

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
    for (size_t i = 0; i < polys.size(); ++i) {
        // Dibujamos aristas
        for (int j = 0; j < polys.get(i)->getNumPoints(); ++j) {
            if(!isSpline(*polys.get(i), polys.get(i)->getPoint(j)))
                dxf->writeLine(
                    *dw,
                    DL_LineData(
                        polys.get(i)->getPoint(j).x,
                        polys.get(i)->getPoint(j).y,
                        0.0,
                        polys.get(i)->getPoint((j + 1) % polys.get(i)->getNumPoints()).x,
                        polys.get(i)->getPoint((j + 1) % polys.get(i)->getNumPoints()).y,
                        0.0),
                    DL_Attributes(
                        "Polygons",
                        256,
                        -1,
                        "BYLAYER"
                    )
                );
        }

        size_t j = 0;
        while (j < polys.get(i)->getNumPointsBC()) {
            // Dibujamos splines
            dxf->writeSpline(
                *dw,
                DL_SplineData(
                    3,
                    8,
                    4,
                    0
                ),
                DL_Attributes(
                    "Polygons",
                    256,
                    -1,
                    "BYLAYER"
                )
            );

            do {
                dxf->writeControlPoint(
                    *dw,
                    DL_ControlPointData(
                        polys.get(i)->getBCPoint(j).x,
                        polys.get(i)->getBCPoint(j).y,
                        0.0
                    )
                );

                ++j;
            } while (j % 4 != 0);
        }
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
#ifdef PRINTDEBUG
    std::cout << "Vertices CH:\n";
#endif
    /*for (size_t i = 0; i < verticesConvexHull.size(); ++i) {
        std::cout << verticesConvexHull.get(i)->x << ", " << verticesConvexHull.get(i)->y << std::endl;
    }*/
#ifdef PRINTDEBUG
    std::cout << convexHullPolygon << std::endl;
#endif
    if (convexHullArea >= 0) {
        // Area del convex Hull
        Point lowerPoint = calculateLowerPoint();

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
        for (int i = 0; i < convexHullPolygon.getNumPoints(); ++i) {
            dxf->writeLine(
                *dw,
                DL_LineData(
                    convexHullPolygon.getPoint(i).x,
                    convexHullPolygon.getPoint(i).y,
                    0.0,
                    convexHullPolygon.getPoint((i + 1) % convexHullPolygon.getNumPoints()).x,
                    convexHullPolygon.getPoint((i + 1) % convexHullPolygon.getNumPoints()).y,
                    0.0),
                DL_Attributes(
                    "ConvexHull",
                    256,
                    -1,
                    "BYLAYER"
                )
            );
        }

        // Dibujamos curvas
        /*int i = 0;
        while (i < convexHullPolygon.getNumPointsBC()) {
            dxf->writeSpline(
                *dw,
                DL_SplineData(
                    3,
                    8,
                    4,
                    0
                ),
                DL_Attributes(
                    "ConvexHull",
                    256,
                    -1,
                    "BYLAYER"
                )
            );

            do {
                dxf->writeControlPoint(
                    *dw,
                    DL_ControlPointData(
                        convexHullPolygon.getBCPoint(i).x,
                        convexHullPolygon.getBCPoint(i).y,
                        0.0
                    )
                );

                ++i;
            } while (i % 4 != 0);
        }*/
    }

    dw->sectionEnd();

    // Objects
    dxf->writeObjects(*dw);
    dxf->writeObjectsEnd(*dw);

    dw->dxfEOF();
    dw->close();

    delete dw;
}

void dxfFilter::printPolygons() const {
    std::cout << "Polygons:\n\t";

    for (size_t i = 0; i < polys.size(); ++i) {
        std::cout << *polys.get(i) << std::endl;
    }
}