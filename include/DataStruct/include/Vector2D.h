/*!
* @file Vector2D.h
* @author Patrice Koehl
* @author Yossi Bokor Bleile
* @date April 2023
* @version 1
* @copyright BSD 3-Clause License.
*/

/*********************************************************************************
 *      The Vector2D class
 *********************************************************************************/

#ifndef VECTOR2D_H
#define VECTOR2D_H

#include "math.h"

namespace correa {
/*********************************************************************************
  class
 *********************************************************************************/

  class Vector2D {

	public:
		// initializes all components to zero
		Vector2D();

		// initializes with specified components
		Vector2D(double x, double y);

		// copy constructor
		Vector2D(const Vector2D& v);

		// access
		double& operator[](int index);

		// math
		Vector2D operator*(double s) const;
		Vector2D operator/(double s) const;
		Vector2D operator+(const Vector2D& v) const;
		Vector2D operator-(const Vector2D& v) const;

		Vector2D operator-() const;
		double operator*(const Vector2D& v) const;

		Vector2D& operator*=(double s);
		Vector2D& operator/=(double s);
		Vector2D& operator+=(const Vector2D& v);
		Vector2D& operator-=(const Vector2D& v);

		// returns Euclidean length
		double norm() const;

		// returns Euclidean length squared
		double norm2() const;

		// normalizes vector
		void normalize();

		// returns unit vector in the direction of this vector
		Vector2D unit() const;

		// members
		double x, y;
  };

  double dot(const Vector2D& u, const Vector2D& v);

/*********************************************************************************
  Constructor (initialize with 0)
 *********************************************************************************/

  Vector2D::Vector2D():
  x(0.0),
  y(0.0)
  {

  }

/*********************************************************************************
  Constructor (with specified coordinates)
 *********************************************************************************/

 Vector2D::Vector2D(double x_, double y_):
  x(x_),
  y(y_)
  {

  }

/*********************************************************************************
  Copy Constructor
 *********************************************************************************/

  Vector2D::Vector2D(const Vector2D& v):
  x(v.x),
  y(v.y)
  {

  }

/*********************************************************************************
  Access specific element
 *********************************************************************************/

 double& Vector2D::operator[](int index)
 {
	return (&x)[index];
 }

/*********************************************************************************
  Math: Multiply vector with a constant
 *********************************************************************************/

 Vector2D Vector2D::operator*(double s) const
 {
	return Vector2D(x*s, y*s);
 }

/*********************************************************************************
  Math: Divide vector with a constant
 *********************************************************************************/

 Vector2D Vector2D::operator/(double s) const
 {
	return Vector2D(x/s, y/s);
 }

/*********************************************************************************
  Math: Add to another vector
 *********************************************************************************/

  Vector2D Vector2D::operator+(const Vector2D& v) const
  {
	return Vector2D(x + v.x, y + v.y);
  }

/*********************************************************************************
  Math: substract another vector
 *********************************************************************************/

  Vector2D Vector2D::operator-(const Vector2D& v) const
  {
	return Vector2D(x - v.x, y - v.y);
  }

/*********************************************************************************
  Math: dot product with another vector
 *********************************************************************************/

  double Vector2D::operator*(const Vector2D& v) const
  {
	return x*v.x + y*v.y;
  }

/*********************************************************************************
  Math: take opposite of vector 
 *********************************************************************************/

 Vector2D Vector2D::operator-() const
 {
	return Vector2D(-x, -y);
 }

/*********************************************************************************
  Math: Multiply with a scalar (concise form)
 *********************************************************************************/

  Vector2D& Vector2D::operator*=(double s)
  {
	x *= s;
	y *= s;

	return *this;
  }

/*********************************************************************************
  Math: Divide with a scalar (concise form)
 *********************************************************************************/

  Vector2D& Vector2D::operator/=(double s)
  {
	x /= s;
	y /= s;

	return *this;
  }

/*********************************************************************************
  Math: Add another vector (concise form)
 *********************************************************************************/

  Vector2D& Vector2D::operator+=(const Vector2D& v)
  {
	x += v.x;
	y += v.y;

	return *this;
  }

/*********************************************************************************
  Math: substract another vector (concise form)
 *********************************************************************************/

  Vector2D& Vector2D::operator-=(const Vector2D& v)
  {
	x -= v.x;
	y -= v.y;

	return *this;
  }

/*********************************************************************************
  Math: norm of the vector
 *********************************************************************************/

  double Vector2D::norm() const
  {
	return std::sqrt(norm2());
  }

/*********************************************************************************
  Math: square of the norm of the vector
 *********************************************************************************/

  double Vector2D::norm2() const
  {
	return dot(*this, *this);
  }

/*********************************************************************************
  Math: Normalize vector
 *********************************************************************************/

  void Vector2D::normalize()
  {
	(*this) /= norm();
  }

/*********************************************************************************
  Math: generate corresponding unit vector
 *********************************************************************************/

  Vector2D Vector2D::unit() const
  {
	return (*this) / norm();
  }

/*********************************************************************************
  Math: dot product
 *********************************************************************************/

  double dot(const Vector2D& u, const Vector2D& v)
  {
	return u.x*v.x + u.y*v.y;
  }
} //end namespace correa
#endif
