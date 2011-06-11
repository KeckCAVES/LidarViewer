/***********************************************************************
FallingSphereProcessor - Point processor class to drop a sphere onto a
LiDAR point cloud along the negative z direction.
Copyright (c) 2010 Oliver Kreylos

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

#ifndef FALLINGSPHEREPROCESSOR_INCLUDED
#define FALLINGSPHEREPROCESSOR_INCLUDED

#include <Math/Math.h>
#include <Math/Constants.h>

#include "LidarTypes.h"

class FallingSphereProcessor
	{
	/* Elements: */
	private:
	Point start; // Starting point of falling sphere
	Scalar radius,radius2; // Radius and squared radius of falling sphere
	Scalar minZ; // Currently lowest height the sphere can fall to
	
	/* Constructors and destructors: */
	public:
	FallingSphereProcessor(const Point& sStart,Scalar sRadius) // Initializes the processor
		:start(sStart),radius(sRadius),radius2(Math::sqr(radius)),
		 minZ(Math::Constants<Scalar>::min)
		{
		}
	
	/* Methods: */
	Box getBox(void) const // Returns the processor's region of interest
		{
		Box result;
		for(int i=0;i<2;++i)
			{
			result.min[i]=start[i]-radius;
			result.max[i]=start[i]+radius;
			}
		result.min[2]=Box::full.min[2];
		result.max[2]=start[2]+radius;
		return result;
		}
	void operator()(const LidarPoint& p) // Processes a LiDAR point
		{
		/* Calculate the LiDAR point's horizontal distance from the starting position: */
		Scalar xy2=Math::sqr(p[0]-start[0])+Math::sqr(p[1]-start[1]);
		
		/* Check if the sphere can touch the point: */
		if(xy2<=radius2)
			{
			/* Calculate the z value at which the sphere touches the point: */
			Scalar z=p[2]+Math::sqrt(radius2-xy2);
			
			/* Update the minimally possible z value: */
			if(minZ<z)
				minZ=z;
			}
		}
	Scalar getMinZ(void) const // Returns the lowest height to which the sphere can drop
		{
		return minZ;
		}
	};

#endif
