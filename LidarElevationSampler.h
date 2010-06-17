/***********************************************************************
LidarElevationSampler - Functor class to evaluate the implicit elevation
function of a 2.5D LiDAR point cloud.
Copyright (c) 2009-2010 Oliver Kreylos

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

#ifndef LIDARELEVATIONSAMPLER_INCLUDED
#define LIDARELEVATIONSAMPLER_INCLUDED

#include <Math/Math.h>
#include <Math/Constants.h>

#include "LidarTypes.h"

class LidarSeparableElevationSampler // Sampler using a separable Lanczos filter for Cartesian grid generation
	{
	/* Elements: */
	private:
	double pos[2]; // Sample point's (x, y) position
	double filterSize[2]; // Filter size, or cell size of target grid, in x and y
	double numLobes; // Number of lobes in Lanczos reconstruction filter
	double accumulator; // Accumulated convolution between point cloud and filter
	double weightSum; // Sum of weights for all accumulated points
	double absWeightSum; // Sum of absolute weights for all accumulated points
	
	/* Constructors and destructors: */
	public:
	LidarSeparableElevationSampler(const double sPos[2],const double sFilterSize[2],int sNumLobes) // Creates an empty sampler
		:numLobes(double(sNumLobes)),
		 accumulator(0.0),weightSum(0.0),absWeightSum(0.0)
		{
		for(int i=0;i<2;++i)
			{
			pos[i]=sPos[i];
			filterSize[i]=sFilterSize[i];
			}
		}
	
	/* Methods: */
	Box getBox(void) const
		{
		Box result;
		for(int i=0;i<2;++i)
			{
			result.min[i]=pos[i]-filterSize[i]*numLobes;
			result.max[i]=pos[i]+filterSize[i]*numLobes;
			}
		result.min[2]=Box::full.min[2];
		result.max[2]=Box::full.max[2];
		return result;
		}
	void operator()(const LidarPoint& p)
		{
		/* Calculate the Lanczos filter weight for the LiDAR point: */
		double weight=1.0;
		for(int i=0;i<2;++i)
			{
			double x=Math::Constants<double>::pi*(double(p[i])-pos[i])/filterSize[i];
			if(x!=0.0)
				{
				weight*=Math::sin(x)/x;
				x/=numLobes;
				weight*=Math::sin(x)/x;
				}
			}
		
		/* Accumulate the weighted point: */
		accumulator+=double(p[2])*weight;
		weightSum+=weight;
		absWeightSum+=Math::abs(weight);
		}
	double getWeightSum(void) const // Returns the sum of point weights
		{
		return weightSum;
		}
	double getAbsWeightSum(void) const // Returns the sum of absolute point weights
		{
		return absWeightSum;
		}
	double getValue(void) const // Returns the convolution result
		{
		/* Return the weighted average: */
		return accumulator/weightSum;
		}
	};

class LidarRadialElevationSampler // Sampler using a radial Lanczos filter
	{
	/* Elements: */
	private:
	double pos[2]; // Sample point's (x, y) position
	double filterRadius; // Filter radius
	double numLobes; // Number of lobes in Lanczos reconstruction filter
	double maxRadius2; // Squared maximum radius of Lanczos reconstruction filter
	double accumulator; // Accumulated convolution between point cloud and filter
	double weightSum; // Sum of weights for all accumulated points
	double absWeightSum; // Sum of absolute weights for all accumulated points
	
	/* Constructors and destructors: */
	public:
	LidarRadialElevationSampler(const double sPos[2],double sFilterRadius,int sNumLobes) // Creates an empty sampler
		:filterRadius(sFilterRadius),numLobes(double(sNumLobes)),maxRadius2(Math::sqr(filterRadius*numLobes)),
		 accumulator(0.0),weightSum(0.0),absWeightSum(0.0)
		{
		for(int i=0;i<2;++i)
			pos[i]=sPos[i];
		}
	
	/* Methods: */
	Box getBox(void) const
		{
		Box result;
		for(int i=0;i<2;++i)
			{
			result.min[i]=pos[i]-filterRadius*numLobes;
			result.max[i]=pos[i]+filterRadius*numLobes;
			}
		result.min[2]=Box::full.min[2];
		result.max[2]=Box::full.max[2];
		return result;
		}
	void operator()(const LidarPoint& p)
		{
		double r=Math::sqr(double(p[0])-pos[0])+Math::sqr(double(p[1])-pos[1]);
		if(r<maxRadius2)
			{
			/* Calculate the Lanczos filter weight for the LiDAR point: */
			r=Math::Constants<double>::pi*Math::sqrt(r)/filterRadius;
			double weight=1.0;
			if(r!=0.0)
				{
				weight*=Math::sin(r)/r;
				r/=numLobes;
				weight*=Math::sin(r)/r;
				}
			
			/* Accumulate the weighted point: */
			accumulator+=double(p[2])*weight;
			weightSum+=weight;
			absWeightSum+=Math::abs(weight);
			}
		}
	double getWeightSum(void) const // Returns the sum of point weights
		{
		return weightSum;
		}
	double getAbsWeightSum(void) const // Returns the sum of absolute point weights
		{
		return absWeightSum;
		}
	double getValue(void) const // Returns the convolution result
		{
		/* Return the weighted average: */
		return accumulator/weightSum;
		}
	};

#endif
