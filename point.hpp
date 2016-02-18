#ifndef POINT_HPP
#define POINT_HPP

#include <cmath>

class Point
{
public:
	Point(double x, double y, double z)
		: x(x), y(y), z(z)
	{	}

	Point()
	{
		x=0; y=0; z=0;
	}

	double norm()
	{
		return sqrt(x*x + y*y + z*z);
	}

	void add(Point p)
	{
		x += p.x;
		y += p.y;
		z += p.z;
	}

	void sub(Point p)
	{
		x -= p.x;
		y -= p.y;
		z -= p.z;
	}

	void scale(double k)
	{
		x *= k;
		y *= k;
		z *= k;
	}

	Point operator-(const Point& b)
	{
		return Point(x - b.x, y - b.y, z - b.z);
	}

	double x;
	double y;
	double z;
};

#endif

