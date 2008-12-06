/***********************************************************************
CylinderFitter - Functor plug-in to fit a cylinder to a set of points
using a Levenberg-Marquardt minimization algorithm.
Copyright (c) 2007-2008 Oliver Kreylos

This file is part of the LiDAR processing and analysis package.

The LiDAR processing and analysis package is free software; you can
redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

The LiDAR processing and analysis package is distributed in the hope
that it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with the LiDAR processing and analysis package; if not, write to the
Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
02111-1307 USA
***********************************************************************/

#ifndef CYLINDERFITTER_INCLUDED
#define CYLINDERFITTER_INCLUDED

#include <vector>
#include <Math/Math.h>
#include <Geometry/Point.h>
#include <Geometry/AffineCombiner.h>
#include <Geometry/Vector.h>

class CylinderFitter
	{
	/* Embedded classes: */
	public:
	typedef double Scalar; // Scalar type
	typedef Geometry::Point<Scalar,3> Point; // Type for target points
	typedef Geometry::Vector<Scalar,3> Vector; // Type for vectors
	static const int dimension=7; // Dimension of the optimization space
	typedef Geometry::ComponentArray<Scalar,dimension> Derivative; // Type for distance function derivatives
	
	/* Elements: */
	private:
	const std::vector<Point>& points; // Reference to vector containing target points
	Point center; // Current estimated point on cylinder axis
	Vector axis; // Current estimated normalized direction of cylinder axis
	Scalar radius; // Current estimated sphere radius
	Point centerSave; // Saved estimated sphere center
	Vector axisSave; // Saved estimated normalized direction of cylinder axis
	Scalar radiusSave; // Saved estimated sphere radius
	
	/* Constructors and destructors: */
	public:
	CylinderFitter(const std::vector<Point>& sPoints,int initialAxis) // Constructs cylinder fitter for given set of target points
		:points(sPoints)
		{
		/* Guess the initial state: */
		Point::AffineCombiner ac;
		for(std::vector<Point>::const_iterator pIt=points.begin();pIt!=points.end();++pIt)
			ac.addPoint(*pIt);
		center=ac.getPoint();
		axis=Vector::zero;
		axis[initialAxis]=Scalar(1);
		radius=Scalar(1);
		};
	
	/* Methods: */
	void setCenter(const Point& newCenter) // Sets the initial estimate for a point on the cylinder's axis
		{
		center=newCenter;
		};
	void setAxis(const Vector& newAxis) // Sets the initial estimate for the direction of the cylinder's axis
		{
		axis=newAxis;
		};
	void setRadius(Scalar newRadius) // Sets the initial estimate for the sphere's radius
		{
		radius=newRadius;
		};
	const Point& getCenter(void) const // Returns the estimated center
		{
		return center;
		};
	const Vector& getAxis(void) const // Returns the estimated axis
		{
		return axis;
		};
	Scalar getRadius(void) const // Returns the estimated radius
		{
		return radius;
		};
	void save(void) // Saves the current estimate
		{
		centerSave=center;
		axisSave=axis;
		radiusSave=radius;
		};
	void restore(void) // Restores the last saved estimate
		{
		center=centerSave;
		axis=axisSave;
		radius=radiusSave;
		};
	size_t getNumPoints(void) const // Returns the number of target points
		{
		return points.size();
		};
	Scalar calcDistance(size_t index) const // Calculates the distance value for the current estimate and the given target point
		{
		return Geometry::mag(Geometry::cross(axis,points[index]-center))-radius;
		};
	Derivative calcDistanceDerivative(size_t index) const // Calculates the derivative of the distance function for the current estimate and the given target point
		{
		Derivative result;
		Scalar dist=Geometry::mag(Geometry::cross(axis,points[index]-center));
		Scalar dist2=axis*(points[index]-center);
		if(dist!=Scalar(0))
			{
			for(int i=0;i<3;++i)
				result[i]=(axis[i]*dist2-(points[index][i]-center[i]))/dist;
			}
		else
			{
			for(int i=0;i<3;++i)
				result[i]=Math::sqrt(Scalar(1)-Math::sqr(axis[i]));
			}
		for(int i=0;i<3;++i)
			result[3+i]=dist2*result[i];
		result[6]=Scalar(-1);
		return result;
		};
	Scalar calcMag(void) const // Returns the magnitude of the current estimate
		{
		return Math::sqrt(Geometry::sqr(center)+Scalar(1)+Math::sqr(radius));
		};
	void increment(Derivative increment) // Increments the current estimate by the given difference vector
		{
		for(int i=0;i<3;++i)
			center[i]-=increment[i];
		for(int i=0;i<3;++i)
			axis[i]-=increment[3+i];
		radius-=increment[6];
		};
	void normalize(void) // Normalizes the current estimate
		{
		axis.normalize();
		center+=axis*((Point::origin-center)*axis);
		if(radius<0.0)
			radius=-radius;
		};
	};

#endif
