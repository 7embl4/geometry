#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <climits>
#include <cassert>

template <typename T>
struct Base {
	T x, y;
	Base() : x(0), y(0) {}
	Base(T x, T y) : x(x), y(y) {}
	Base(const Base& other) : x(other.x), y(other.y) {}

	// polar angle in [0, 2*pi)
	double getPolarAngle() const {
		double angle = atan2(this->y, this->x);
		if (angle < 0) {
			angle += 2 * M_PI;
		}
		return angle;
	}

	// length of vector
	double len() const {
		return sqrt(x * x + y * y);
	}
};

// dot product
template <typename T>
T operator*(const Base<T>& v1, const Base<T>& v2) {
	return v1.x * v2.x + v1.y * v2.y;
}

// cross product
template <typename T>
T operator%(const Base<T>& v1, const Base<T>& v2) {
	return v1.x * v2.y - v2.x * v1.y;
}

template <typename T>
Base<T> operator+(const Base<T>& l, const Base<T>& r) {
	return Base<T>(l.x + r.x, l.y + r.y);
}

template <typename T>
Base<T> operator-(const Base<T>& l, const Base<T>& r) {
	return Base<T>(l.x - r.x, l.y - r.y);
}

template <typename T>
Base<T> operator*(const Base<T>& vec, int n) {
	return Base<T>(vec.x * n, vec.y * n);
}

template <typename T>
Base<T> operator/(const Base<T>& vec, int n) {
	assert(n != 0);
	return Base<T>(vec.x / n, vec.y / n);
}

template <typename T>
Base<T> operator-(const Base<T>& vec) {
	return Base<T>(-vec.x, -vec.y);
}

template <typename T>
bool operator==(const Base<T>& l, const Base<T>& r) {
	return l.x == r.x && l.y == r.y;
}
template <typename T>
bool operator!=(const Base<T>& l, const Base<T>& r) {
	return !(l == r);
}

template <typename T>
std::istream& operator>>(std::istream& in, Base<T>& b) {
	int x, y;
	in >> x >> y;
	b.x = x;
	b.y = y;
	return in;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const Base<T>& b) {
	out << b.x << ' ' << b.y;
	return out;
}

// useful usings
using vec2i   = Base<int>;
using point2i = Base<int>;
using vec2f   = Base<double>;
using point2f = Base<double>;

/* ----------------------------- solving problems ----------------------------- */

// partly solved
void interPointOfTwoLines(
	const point2i& begin1, const point2i& end1,
	const point2i& begin2, const point2i& end2
) {
	int a1 = begin1.y - end1.y;
	int b1 = end1.x - begin1.x;
	int c1 = begin1.x * end1.y - end1.x * begin1.y;

	int a2 = begin2.y - end2.y;
	int b2 = end2.x - begin2.x;
	int c2 = begin2.x * end2.y - end2.x * begin2.y;
	
	if (a1 * b2 - b1 * a2 == 0) {
		if (a1 * c2 - a2 * c1 == 0) {
			std::cout << 2;
		}
		else {
			std::cout << 0;
		}
	}
	else {
		double x = (b1 * c2 - c1 * b2) / (double)(a1 * b2 - b1 * a2);
		double y = (c1 * a2 - a1 * c2) / (double)(a1 * b2 - b1 * a2);
		std::cout << 1 << std::setprecision(15) << x << ' ' << y;
	}
}

std::string isPointsOnSameSide(const point2i& p1, const point2i& p2, int a, int b, int c) {
	int first_point = a * p1.x + b * p1.y + c;
	int second_point = a * p2.x + b * p2.y + c;
	if (first_point * second_point > 0) {
		return "YES";
	}
	return "NO";
}

double angleOfTwoVectors(const vec2i& v1, const vec2i& v2) {
	double angle = acos((v1 * v2) / (v1.len() * v2.len()));
	return angle < 0 ? M_PI - angle : angle;
}

std::vector<double> findPerpendicularThroughPoint(int a, int b, int c, const point2i& p) {
	double a_perp = -b;
	double b_perp = a;
	double c_perp = -(a_perp * p.x + b_perp * p.y);
	return { a_perp, b_perp, c_perp };
}

double calcTriangleSurface(const point2i& p1, const point2i& p2, const point2i& p3) {
	vec2i v1 = p3 - p1;
	vec2i v2 = p2 - p1;
	return fabs(v1 % v2) / 2;
}

double calcPolarDist(int r1, int angle1, int r2, int angle2) {\
	// find angle of two radius-vectors
	int angle_of_lines = abs(angle1 - angle2);
	if (angle_of_lines > 180) {
		angle_of_lines = 360 - angle_of_lines;
	}

	// find dist using cosinus theorem 
	double angle_in_radians = angle_of_lines * (M_PI / 180);
	return sqrt(r1 * r1 + r2 * r2 - 2 * r1 * r2 * cos(angle_in_radians));
}

size_t calcFine(const std::vector<point2f>& route) {
	if (route.size() < 3) {
		return 0;
	}

	size_t fine = 0;
	for (size_t i = 2; i != route.size(); ++i) {
		vec2f prev_vec = route[i - 1] - route[i - 2];
		vec2f curr_vec = route[i] - route[i - 1];
		if (prev_vec % curr_vec > 0) {
			++fine;
		}
	}
	return fine;
}

std::string isPointInSection(const point2i& point, const point2i& begin, const point2i& end) {
	if ((end - begin) % (begin - point) == 0 && (end - point) * (begin - point) <= 0) {
		return "YES";
	}
	return "NO";
}

std::vector<int> lineThroughTwoPoints(const point2i& p1, const point2i& p2) {
	int a = p1.y - p2.y;
	int b = p2.x - p1.x;
	int c = p1.x * p2.y - p2.x * p1.y;
	return { a, b, c };
}

void parallelLine(int a, int b, int c, int r) {
	std::cout << a << ' ' << b << ' ' << std::setprecision(10) << c - r * sqrt(a * a + b * b);
}

// partly solved
std::string isPointInConvexPolygon(const std::vector<point2i>& polygon, const point2i& point) {
	point2i p0 = polygon[0];
	size_t l = 1, r = polygon.size() - 1;

	// check if point is out of the polygon
	if (((polygon[l] - p0) % (point - p0)) * ((point - p0) % (polygon[r] - p0)) < 0) {
		return "NO";
	}
	// check if point is on an edge of polygon
	else if ((polygon[l] - p0) % (point - p0) == 0) {
		return isPointInSection(point, p0, polygon[l]);
	}
	else if ((point - p0) % (polygon[r] - p0) == 0) {
		return isPointInSection(point, p0, polygon[r]);
	}

	// find angle where point is
	while (l <= r) {
		size_t m = (l + r) / 2;
		if ((polygon[m] - p0) % (point - p0) < 0) {
			l = m + 1;
		}
		else if ((polygon[m] - p0) % (point - p0) > 0) {
			r = m - 1;
		}
		else {
			// point is on a vector
			// so check if point is in section (p0, polygon[m])
			return isPointInSection(point, p0, polygon[m]);
		}
	}

	// check if point is in triangle (p0, polygon[l], polygon[r])
	if (((polygon[l] - p0) % (point - p0)) * ((point - p0) % (polygon[r] - p0)) >= 0
		&& ((p0 - polygon[l]) % (point - polygon[l])) * ((point - polygon[l]) % (polygon[r] - polygon[l])) >= 0) {
		return "YES";
	}
	return "NO";
}

std::vector<int> lineWithPointAndNorm(const point2i& point, const vec2i& norm) {
	return { norm.x, norm.y, -(point.x * norm.x + point.y * norm.y) };
}

template <typename T>
double distFromPointToLine(const Base<T>& point, int a, int b, int c) {
	return fabs(a * point.x + b * point.y + c) / sqrt(a * a + b * b);
}

double distFromPointToSection(const point2i& point, const point2i& begin, const point2i& end) {
	if ((end - begin) * (point - begin) < 0) {
		return sqrt(pow(begin.x - point.x, 2) + pow(begin.y - point.y, 2));
	}
	else if ((begin - end) * (point - end) < 0) {
		return sqrt(pow(end.x - point.x, 2) + pow(end.y - point.y, 2));
	}
	return fabs((end - point) % (begin - point)) / (end - begin).len();
}

std::string isPointOnLine(const point2i& point, int a, int b, int c) {
	if (a * point.x + b * point.y + c == 0) {
		return "YES";
	}
	return "NO";
}

// partly solved
point2f heightsIntersectionPoint(const point2i& p1, const point2i& p2, const point2i& p3) {
	int a1 = p2.y - p3.y;
	int b1 = p3.x - p2.x;
	std::vector<int> first_line = lineWithPointAndNorm(p1, vec2i(-b1, a1));

	int a2 = p1.y - p2.y;
	int b2 = p2.x - p1.x;
	std::vector<int> second_line = lineWithPointAndNorm(p3, vec2i(-b2, a2));
	
	double x = (first_line[1] * second_line[2] - first_line[2] * second_line[1])
		/ (double)(first_line[0] * second_line[1] - first_line[1] * second_line[0]);
	double y = (first_line[2] * second_line[0] - first_line[0] * second_line[2])
		/ (double)(first_line[0] * second_line[1] - first_line[1] * second_line[0]);

	return point2f(x, y);
}

std::string isPointOnRay(const point2i& point, const point2i& begin, const point2i& end) {
	vec2i b_to_e = end - begin;
	vec2i b_to_p = point - begin;
	if (b_to_p * b_to_e >= 0 && b_to_p % b_to_e == 0) {
		return "YES";
	}
	return "NO";
}

point2f mediansIntersectionPoint(const point2i& p1, const point2i& p2, const point2i& p3) {
	double x = (p1.x + p2.x + p3.x) / 3.f;
	double y = (p1.y + p2.y + p3.y) / 3.f;
	return point2f(x, y);
}

double distFromPointToRay(const point2i& point, const point2i& begin, const point2i& end) {
	vec2i b_to_e = end - begin;
	vec2i b_to_p = point - begin;
	if (b_to_p * b_to_e < 0) {
		return sqrt(pow(point.x - begin.x, 2) + pow(point.y - begin.y, 2));
	}
	return abs((end - point) % (begin - point)) / (end - begin).len();
}

std::string cabinetThroughDoorway(size_t a, size_t b, size_t c, size_t x, size_t y) {
	size_t doorway_area = x * y;
	if (a * b <= doorway_area || a * c <= doorway_area || b * c <= doorway_area) {
		return "YES";
	}
	return "NO";
}

double rotatingDoor(size_t a, size_t b, size_t c) {
	std::vector<size_t> sizes = { a, b, c };
	std::sort(sizes.begin(), sizes.end());
	return sqrt(sizes[0] * sizes[0] + sizes[1] * sizes[1]);
}

std::vector<point2f> makeSquare(const point2i& a, const point2i& b) {
	point2f center((a.x + b.x) / 2.f, (a.y + b.y) / 2.f);
	vec2f norm(a.y - b.y, b.x - a.x);

	vec2f vec_from_a_to_center(center.x - a.x, center.y - a.y);
	point2f new_norm = norm / norm.len() * vec_from_a_to_center.len();
	if (new_norm != point2f(0, 0)) {
		norm = new_norm;
	}
	else {
		norm = -norm / norm.len() * vec_from_a_to_center.len();
	}
	std::cout << norm << '\n';

	std::vector<point2f> result(2);
	result[0] = point2f(norm.x + center.x, norm.y + center.y);
	result[1] = point2f(-norm.x + center.x, -norm.y + center.y);
	
	return result;
}

int sgn(float n) {
	if (n == 0) {
		return 0;
	}
	else if (n < 0) {
		return -1;
	}
	else {
		return 1;
	}
}

bool isConvexPolygon(const std::vector<point2i>& polygon) {
	int sign = 0;
	for (size_t i = 1; i != polygon.size() - 1; ++i) {
		point2i vec1 = polygon[i] - polygon[i - 1];
		point2i vec2 = polygon[i + 1] - polygon[i];
		if (sign != 0) {
			if (sgn(vec1 % vec2) != sign) {
				return false;
			}
		}
		else {
			if (sgn(vec1 % vec2) != 0) {
				sign = sgn(vec1 % vec2);
			}
		}
		
	}
	return true;
}

// partly solved
point2f interPointOfTwoLines(int a1, int b1, int c1, int a2, int b2, int c2) {
	point2f interpoint;
	interpoint.x = (b1 * c2 - b2 * c1) / (a1 * b2 - a2 * b1);
	if (b2 != 0) {
		interpoint.y = -(c2 + a2 * interpoint.x) / b2;
	}
	else {
		interpoint.y = -(c1 + a1 * interpoint.x) / b1;
	}
	return interpoint;
}

// partly solved
point2f makeReflection(int a, int b, int c, const point2i& point) {
	vec2f norm(a, b);
	double dist = distFromPointToLine(point, a, b, c);
	vec2f vec1(norm / norm.len() * 2 * dist);
	vec2f vec2(-norm / norm.len() * 2 * dist);
	
	point2f p1(point.x + vec1.x, point.y + vec1.y);
	point2f p2(point.x + vec2.x, point.y + vec2.y);

	if (distFromPointToLine(p1, a, b, c) == dist) {
		return p1;
	}
	return p2;
}

std::pair<point2f, point2f> changeCoordinateSystem(
	double a, double b, double c, double d, 
	double x1, double y1, double x2, double y2
) {
	double x_c = -a / (c - a);  // center of the first system in the second system
	double y_c = -b / (d - b);

	point2f from_first_in_second(x1 * (1 / (c - a)) + x_c, y1 * (1 / (d - b)) + y_c);
	point2f from_second_in_first((c - a) * x2 + a, (d - b) * y2 + b);

	return { from_first_in_second, from_second_in_first };
}

bool is_not_convex(const point2i& prelast_point, const point2i& last_point, const point2i& curr_point)
{
	vec2i v1 = last_point - prelast_point;
	vec2i v2 = curr_point - last_point;
	return v1 % v2 < 0;
}

static std::vector<point2i> convexHullGraham(std::vector<point2i> points)
{
	// find first point
	size_t p0_index = 0;
	point2i p0 = points[0];
	for (size_t i = 1; i != points.size(); ++i) {
		if (points[i].x < p0.x || points[i].x == p0.x && points[i].y < p0.y) {
			p0 = points[i];
			p0_index = i;
		}
	}

	// sort other points
	points.erase(points.begin() + p0_index);
	std::sort(
		points.begin(),
		points.end(),
		[&p0](const point2i& p1, const point2i& p2) {
			vec2i v1 = p1 - p0;
			vec2i v2 = p2 - p0;
			return v1 % v2;
		}
	);

	// find convex hull
	std::vector<point2i> convex_hull;
	convex_hull.push_back(p0);
	for (size_t i = 0; i != points.size(); ++i) {
		while (convex_hull.size() > 1 
			&& is_not_convex(convex_hull[convex_hull.size() - 2], convex_hull[convex_hull.size() - 1], points[i])) 
		{
			convex_hull.pop_back();
		}
		convex_hull.push_back(points[i]);
	}

	return convex_hull;
}
