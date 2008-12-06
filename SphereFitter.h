/***********************************************************************
SphereFitter - Functor plug-in to fit a sphere to a set of points using
a Levenberg-Marquardt minimization algorithm.
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

#ifndef SPHEREFITTER_INCLUDED
#define SPHEREFITTER_INCLUDED

#include <vector>
#include <Math/Math.h>
#include <Geometry/Point.h>
#include <Geometry/AffineCombiner.h>

class SphereFitter
	{
	/* Embedded classes: */
	public:
	typedef double Scalar; // Scalar type
	typedef Geometry::Point<Scalar,3> Point; // Type for target points
	static const int dimension=4; // Dimension of the optimization space
	typedef Geometry::ComponentArray<Scalar,dimension> Derivative; // Type for distance function derivatives
	
	/* Elements: */
	private:
	const std::vector<Point>& points; // Reference to vector containing target points
	Point center; // Current estimated sphere center
	Scalar radius; // Current estimated sphere radius
	Point centerSave; // Saved estimated sphere center
	Scalar radiusSave; // Saved estimated sphere radius
	
	/* Constructors and destructors: */
	public:
	SphereFitter(const std::vector<Point>& sPoints) // Constructs sphere fitter for given set of target points
		:points(sPoints)
		{
		/* Guess the initial state: */
		Point::AffineCombiner ac;
		for(std::vector<Point>::const_iterator pIt=points.begin();pIt!=points.end();++pIt)
			ac.addPoint(*pIt);
		center=ac.getPoint();
		radius=Scalar(1);
		};
	
	/* Methods: */
	void setCenter(const Point& newCenter) // Sets the initial estimate for the sphere's center
		{
		center=newCenter;
		};
	void setRadius(Scalar newRadius) // Sets the initial estimate for the sphere's radius
		{
		radius=newRadius;
		};
	const Point& getCenter(void) const // Returns the estimated center
		{
		return center;
		};
	Scalar getRadius(void) const // Returns the estimated radius
		{
		return radius;
		};
	void save(void) // Saves the current estimate
		{
		centerSave=center;
		radiusSave=radius;
		};
	void restore(void) // Restores the last saved estimate
		{
		center=centerSave;
		radius=radiusSave;
		};
	size_t getNumPoints(void) const // Returns the number of target points
		{
		return points.size();
		};
	Scalar calcDistance(size_t index) const // Calculates the distance value for the current estimate and the given target point
		{
		return Geometry::dist(points[index],center)-radius;
		};
	Derivative calcDistanceDerivative(size_t index) const // Calculates the derivative of the distance function for the current estimate and the given target point
		{
		Derivative result;
		Scalar dist=Geometry::dist(points[index],center);
		for(int i=0;i<3;++i)
			result[i]=-(points[index][i]-center[i])/dist;
		result[3]=Scalar(-1);
		return result;
		};
	Scalar calcMag(void) const // Returns the magnitude of the current estimate
		{
		return Math::sqrt(Geometry::sqr(center)+Math::sqr(radius));
		};
	void increment(Derivative increment) // Increments the current estimate by the given difference vector
		{
		for(int i=0;i<3;++i)
			center[i]-=increment[i];
		radius-=increment[3];
		};
	void normalize(void) // Normalizes the current estimate
		{
		if(radius<0.0)
			radius=-radius;
		};
	};

#endif
