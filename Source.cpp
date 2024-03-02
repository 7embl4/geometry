#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>

struct Base {
	int x, y;
	Base() : x(0), y(0) {}
	Base(int x, int y) : x(x), y(y) {}
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
int operator*(const Base& v1, const Base& v2) {
	return v1.x * v2.x + v1.y * v2.y;
}

// cross product
int operator%(const Base& v1, const Base& v2) {
	return v1.x * v2.y - v2.x * v1.y;
}

Base operator+(const Base& l, const Base& r) {
	return Base(l.x + r.x, l.y + r.y);
}

Base operator-(const Base& l, const Base& r) {
	return Base(l.x - r.x, l.y - r.y);
}

std::istream& operator>>(std::istream& in, Base& b) {
	int x, y;
	in >> x >> y;
	b.x = x;
	b.y = y;
	return in;
}

std::ostream& operator<<(std::ostream& out, const Base& b) {
	out << b.x << ' ' << b.y;
	return out;
}

using vec2 = Base;
using point2 = Base;

std::string isPointsOnSameSide(const point2& p1, const point2& p2, int a, int b, int c) {
	int first_point = a * p1.x + b * p1.y + c;
	int second_point = a * p2.x + b * p2.y + c;
	if (first_point * second_point > 0) {
		return "YES";
	}
	return "NO";
}

double angleOfTwoVectors(const vec2& v1, const vec2& v2) {
	double angle = acos((v1 * v2) / (v1.len() * v2.len()));
	return angle < 0 ? M_PI - angle : angle;
}

std::vector<double> findPerpendicularThroughPoint(int a, int b, int c, const point2& p) {
	double a_perp = -b;
	double b_perp = a;
	double c_perp = -(a_perp * p.x + b_perp * p.y);
	return { a_perp, b_perp, c_perp };
}

double calcTriangleSurface(const point2& p1, const point2& p2, const point2& p3) {
	vec2 v1 = p3 - p1;
	vec2 v2 = p2 - p1;
	return fabs(v1 % v2) / 2;
}

int main() {
	

	return 0;
}