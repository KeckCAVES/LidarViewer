/***********************************************************************
PointBasedLightingShader - Class to maintain a GLSL point-based lighting
shader that tracks the current OpenGL lighting state.
Copyright (c) 2008-2013 Oliver Kreylos

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

#include "PointBasedLightingShader.h"

#include <string>
#include <iostream>
#include <Misc/PrintInteger.h>
#include <Misc/ThrowStdErr.h>
#include <GL/gl.h>
#include <GL/GLLightTracker.h>
#include <GL/GLClipPlaneTracker.h>
#include <GL/GLContextData.h>
#include <GL/Extensions/GLARBShaderObjects.h>
#include <GL/Extensions/GLARBVertexShader.h>
#include <GL/Extensions/GLARBGeometryShader4.h>
#include <GL/Extensions/GLARBFragmentShader.h>

/*****************************************
Methods of class PointBasedLightingShader:
*****************************************/

void PointBasedLightingShader::compileShader(void)
	{
	const GLLightTracker& lt=*(contextData.getLightTracker());
	const GLClipPlaneTracker& cpt=*(contextData.getClipPlaneTracker());
	
	std::string vertexShaderDefines;
	std::string vertexShaderFunctions;
	std::string vertexShaderMain;
	
	if(usePlaneDistance)
		{
		/* Create the plane distance mapping uniforms: */
		vertexShaderDefines+="\
			uniform vec4 planeDistancePlane;\n\
			uniform sampler1D planeDistanceMap;\n\
			\n";
		}
	
	/* Create the main vertex shader starting boilerplate: */
	vertexShaderMain+="\
		void main()\n\
			{\n\
			/* Compute the vertex position in eye coordinates: */\n\
			vec4 vertexEc=gl_ModelViewMatrix*gl_Vertex;\n\
			\n\
			/* Compute the normal vector in eye coordinates: */\n\
			vec3 normalEc=normalize(gl_NormalMatrix*gl_Normal);\n\
			\n\
			/* Let the normal vector always point towards the eye: */\n\
			normalEc=faceforward(normalEc,normalEc,vertexEc.xyz);\n\
			\n";
	
	/* Get the material components: */
	if(usePlaneDistance)
		{
		#ifdef LIDARVIEWER_VISUALIZE_WATER
		
		vertexShaderMain+="\
			/* Calculate the distance from the water surface: */\n\
			float planeDist=dot(planeDistancePlane,gl_Vertex);\n\
			vec4 ambient,diffuse;\n\
			if(planeDist<=0.5)\n\
				{\n\
				/* Get the material properties from the plane distance texture: */\n\
				ambient=texture1D(planeDistanceMap,planeDist);\n\
				diffuse=ambient;\n\
				}\n\
			else\n\
				{\n";
		
		if(usePointColors)
			{
			vertexShaderMain+="\
				/* Get the material properties from the current color: */\n\
				ambient=gl_Color;\n\
				diffuse=gl_Color;\n";
			}
		else
			{
			vertexShaderMain+="\
				/* Get the material properties from the material state: */\n\
				ambient=gl_FrontMaterial.ambient;\n\
				diffuse=gl_FrontMaterial.diffuse;\n";
			}
		
		vertexShaderMain+="\
				}\n";
		
		#else
		
		vertexShaderMain+="\
			/* Get the material properties from the plane distance texture: */\n\
			float planeDist=dot(planeDistancePlane,gl_Vertex);\n\
			vec4 ambient=texture1D(planeDistanceMap,planeDist);\n\
			vec4 diffuse=ambient;\n";
		
		#endif
		}
	else if(usePointColors)
		{
		vertexShaderMain+="\
			/* Get the material properties from the current color: */\n\
			vec4 ambient=gl_Color;\n\
			vec4 diffuse=gl_Color;\n";
		}
	else
		{
		vertexShaderMain+="\
			/* Get the material properties from the material state: */\n\
			vec4 ambient=gl_FrontMaterial.ambient;\n\
			vec4 diffuse=gl_FrontMaterial.diffuse;\n";
		}
	vertexShaderMain+="\
			vec4 specular=gl_FrontMaterial.specular;\n\
			float shininess=gl_FrontMaterial.shininess;\n\
			\n";
	
	/* Continue the main vertex shader: */
	vertexShaderMain+="\
			/* Calculate global ambient light term: */\n\
			vec4 ambientDiffuseAccum=gl_LightModel.ambient*ambient;\n\
			vec4 specularAccum=vec4(0.0,0.0,0.0,0.0);\n\
			\n\
			/* Accumulate all enabled light sources: */\n";
	
	/* Create light application functions for all enabled light sources: */
	for(int lightIndex=0;lightIndex<lt.getMaxNumLights();++lightIndex)
		if(lt.getLightState(lightIndex).isEnabled())
			{
			/* Create the light accumulation function: */
			vertexShaderFunctions+=lt.createAccumulateLightFunction(lightIndex);
			
			/* Call the light application function from the shader's main function: */
			vertexShaderMain+="\
				accumulateLight";
			char liBuffer[12];
			vertexShaderMain.append(Misc::print(lightIndex,liBuffer+11));
			vertexShaderMain+="(vertexEc,normalEc,ambient,diffuse,specular,shininess,ambientDiffuseAccum,specularAccum);\n";
			}
	
	/* Continue the main vertex shader: */
	vertexShaderMain+="\
			\n\
			/* Compute final vertex color: */\n\
			gl_FrontColor=ambientDiffuseAccum+specularAccum;\n\
			\n";
	
	/* Insert code to calculate the vertex' position relative to all user-specified clipping planes: */
	vertexShaderMain+=cpt.createCalcClipDistances("vertexEc");
	
	/* Finish the main vertex shader: */
	if(useSplatting)
		{
		/* Create the splatting varyings: */
		vertexShaderDefines+="\
			varying vec3 normal;\n\
			varying float splatSize;\n\
			\n";
		
		vertexShaderMain+="\
				/* Pass normal vector to geometry shader: */\n\
				normal=normalEc;\n\
				splatSize=length(gl_Normal);\n\
				\n\
				/* Pass eye coordinate vertex position to geometry shader: */\n\
				gl_Position=vertexEc;\n\
				}\n";
		}
	else
		{
		vertexShaderMain+="\
				/* Use standard vertex position: */\n\
				gl_Position=ftransform();\n\
				}\n";
		}
	
	/* Compile the vertex shader: */
	std::string vertexShaderSource=vertexShaderDefines+vertexShaderFunctions+vertexShaderMain;
	glCompileShaderFromString(vertexShader,vertexShaderSource.c_str());
	
	if(useSplatting)
		{
		if(!geometryShaderAttached)
			{
			/* Attach the geometry shader to the program object: */
			glAttachObjectARB(programObject,geometryShader);
			geometryShaderAttached=true;
			}
		
		/* Compile the surfel generation geometry shader: */
		const char* geometryShaderSource="\
			#version 120\n\
			#extension GL_ARB_geometry_shader4: enable\n\
			\n\
			uniform float surfelSize;\n\
			\n\
			varying in vec3 normal[];\n\
			varying in float splatSize[];\n\
			\n\
			void main()\n\
				{\n\
				/* Calculate quad base vectors based on the eye-coordinate vertex position and normal: */\n\
				vec3 x;\n\
				if(abs(normal[0].x)<abs(normal[0].y)&&abs(normal[0].x)<abs(normal[0].z))\n\
					x=normalize(vec3(0.0,normal[0].z,-normal[0].y));\n\
				else if(abs(normal[0].y)<abs(normal[0].z))\n\
					x=normalize(vec3(normal[0].z,0.0,-normal[0].x));\n\
				else\n\
					x=normalize(vec3(normal[0].y,-normal[0].x,0.0));\n\
				x*=splatSize[0]*surfelSize*1.41421356;\n\
				vec3 y=cross(normal[0],x);\n\
				\n\
				/* Emit the quad's four vertices: */\n\
				gl_TexCoord[0].st=vec2(-1.0,-1.0);\n\
				gl_FrontColor=gl_FrontColorIn[0];\n\
				gl_Position=gl_ProjectionMatrix*(gl_PositionIn[0]+vec4(x,0.0));\n\
				EmitVertex();\n\
				\n\
				gl_TexCoord[0].st=vec2(1.0,-1.0);\n\
				gl_FrontColor=gl_FrontColorIn[0];\n\
				gl_Position=gl_ProjectionMatrix*(gl_PositionIn[0]+vec4(y,0.0));\n\
				EmitVertex();\n\
				\n\
				gl_TexCoord[0].st=vec2(-1.0,1.0);\n\
				gl_FrontColor=gl_FrontColorIn[0];\n\
				gl_Position=gl_ProjectionMatrix*(gl_PositionIn[0]-vec4(y,0.0));\n\
				EmitVertex();\n\
				\n\
				gl_TexCoord[0].st=vec2(1.0,1.0);\n\
				gl_FrontColor=gl_FrontColorIn[0];\n\
				gl_Position=gl_ProjectionMatrix*(gl_PositionIn[0]-vec4(x,0.0));\n\
				EmitVertex();\n\
				}\n";
		glCompileShaderFromString(geometryShader,geometryShaderSource);
		
		/* Set the geometry shader's parameters: */
		glProgramParameteriARB(programObject,GL_GEOMETRY_VERTICES_OUT_ARB,4);
		glProgramParameteriARB(programObject,GL_GEOMETRY_INPUT_TYPE_ARB,GL_POINTS);
		glProgramParameteriARB(programObject,GL_GEOMETRY_OUTPUT_TYPE_ARB,GL_TRIANGLE_STRIP);
		
		/* Compile the surfel fragment shader: */
		const char* fragmentShaderSource=
			"\
			void main()\n\
				{\n\
				/* Discard fragments outside a unit-radius circle as defined by texture coordinates: */\n\
				if(dot(gl_TexCoord[0].xy,gl_TexCoord[0].xy)>1.0)\n\
					discard;\n\
				\n\
				gl_FragColor=gl_Color;\n\
				}\n";
		glCompileShaderFromString(fragmentShader,fragmentShaderSource);
		}
	else
		{
		if(geometryShaderAttached)
			{
			/* Detach the geometry shader from the program object: */
			glDetachObjectARB(programObject,geometryShader);
			geometryShaderAttached=false;
			}
		
		/* Compile the standard fragment shader: */
		const char* fragmentShaderSource=
			"\
			void main()\n\
				{\n\
				gl_FragColor=gl_Color;\n\
				}\n";
		glCompileShaderFromString(fragmentShader,fragmentShaderSource);
		}
	
	/* Link the program object: */
	glLinkProgramARB(programObject);
	
	/* Check if the program linked successfully: */
	GLint linkStatus;
	glGetObjectParameterivARB(programObject,GL_OBJECT_LINK_STATUS_ARB,&linkStatus);
	if(!linkStatus)
		{
		/* Get some more detailed information: */
		GLcharARB linkLogBuffer[2048];
		GLsizei linkLogSize;
		glGetInfoLogARB(programObject,sizeof(linkLogBuffer),&linkLogSize,linkLogBuffer);
		
		/* Signal an error: */
		Misc::throwStdErr("Error \"%s\" while linking shader program",linkLogBuffer);
		}
	
	if(useSplatting)
		{
		/* Get the locations of the uniform variables: */
		surfelSizeLocation=glGetUniformLocationARB(programObject,"surfelSize");
		}
	
	if(usePlaneDistance)
		{
		/* Get the locations of the uniform variables: */
		planeDistancePlaneLocation=glGetUniformLocationARB(programObject,"planeDistancePlane");
		planeDistanceMapLocation=glGetUniformLocationARB(programObject,"planeDistanceMap");
		}
	}

PointBasedLightingShader::PointBasedLightingShader(GLContextData& sContextData)
	:contextData(sContextData),
	 haveGeometryShaders(false),
	 lightStateVersion(0),clipPlaneStateVersion(0),shaderSettingsVersion(0),
	 settingsVersion(1),
	 usePlaneDistance(false),
	 usePointColors(false),
	 useSplatting(false),
	 vertexShader(0),geometryShader(0),fragmentShader(0),programObject(0),geometryShaderAttached(false)
	{
	/* Check for the required OpenGL extensions: */
	if(!GLARBShaderObjects::isSupported())
		Misc::throwStdErr("GLShader::GLShader: GL_ARB_shader_objects not supported");
	if(!GLARBVertexShader::isSupported())
		Misc::throwStdErr("GLShader::GLShader: GL_ARB_vertex_shader not supported");
	if(!GLARBFragmentShader::isSupported())
		Misc::throwStdErr("GLShader::GLShader: GL_ARB_fragment_shader not supported");
	
	/* Initialize the required extensions: */
	GLARBShaderObjects::initExtension();
	GLARBVertexShader::initExtension();
	GLARBFragmentShader::initExtension();
	
	/* Create the vertex and fragment shaders: */
	vertexShader=glCreateShaderObjectARB(GL_VERTEX_SHADER_ARB);
	fragmentShader=glCreateShaderObjectARB(GL_FRAGMENT_SHADER_ARB);
	
	/* Create the program object: */
	programObject=glCreateProgramObjectARB();
	glAttachObjectARB(programObject,vertexShader);
	glAttachObjectARB(programObject,fragmentShader);
	
	/* Check for the optional geometry shader extension: */
	if(GLARBGeometryShader4::isSupported())
		{
		/* Initialize the extension: */
		haveGeometryShaders=true;
		GLARBGeometryShader4::initExtension();
		
		/* Create the geometry shader: */
		geometryShader=glCreateShaderObjectARB(GL_GEOMETRY_SHADER_ARB);
		}
	}

PointBasedLightingShader::~PointBasedLightingShader(void)
	{
	glDeleteObjectARB(programObject);
	glDeleteObjectARB(vertexShader);
	if(haveGeometryShaders)
		glDeleteObjectARB(geometryShader);
	glDeleteObjectARB(fragmentShader);
	}

void PointBasedLightingShader::setUsePlaneDistance(bool newUsePlaneDistance)
	{
	if(usePlaneDistance!=newUsePlaneDistance)
		{
		/* Update the state: */
		usePlaneDistance=newUsePlaneDistance;
		++settingsVersion;
		}
	}

void PointBasedLightingShader::setUsePointColors(bool newUsePointColors)
	{
	if(usePointColors!=newUsePointColors)
		{
		usePointColors=newUsePointColors;
		++settingsVersion;
		}
	}

void PointBasedLightingShader::setUseSplatting(bool newUseSplatting)
	{
	/* Disable splatting if geometry shaders are not supported: */
	newUseSplatting=newUseSplatting&&haveGeometryShaders;
	
	if(useSplatting!=newUseSplatting)
		{
		useSplatting=newUseSplatting;
		++settingsVersion;
		}
	}

void PointBasedLightingShader::enable(void)
	{
	try
		{
		/* Re-compile the shader if it is out of line with current state: */
		const GLLightTracker& lt=*(contextData.getLightTracker());
		const GLClipPlaneTracker& cpt=*(contextData.getClipPlaneTracker());
		
		if(lightStateVersion!=lt.getVersion()||clipPlaneStateVersion!=cpt.getVersion()||shaderSettingsVersion!=settingsVersion)
			{
			/* Rebuild the shader: */
			compileShader();
			
			/* Mark the shader as up-to-date: */
			lightStateVersion=lt.getVersion();
			clipPlaneStateVersion=cpt.getVersion();
			shaderSettingsVersion=settingsVersion;
			}
		
		/* Enable the shader: */
		glUseProgramObjectARB(programObject);
		}
	catch(std::runtime_error err)
		{
		std::cerr<<"Disabling lighting shader due to exception "<<err.what()<<std::endl;
		}
	}

void PointBasedLightingShader::setSurfelSize(float surfelSize)
	{
	if(useSplatting)
		{
		/* Set the surfel size uniform variable: */
		glUniformARB(surfelSizeLocation,surfelSize);
		}
	}

void PointBasedLightingShader::setDistancePlane(int textureUnit,const PointBasedLightingShader::Plane& distancePlane,double distancePlaneScale) const
	{
	if(usePlaneDistance)
		{
		/* Set the plane equation variable: */
		GLfloat planeEq[4];
		for(int i=0;i<3;++i)
			planeEq[i]=GLfloat(distancePlane.getNormal()[i]/distancePlaneScale);
		planeEq[3]=GLfloat(0.5-distancePlane.getOffset()/distancePlaneScale);
		glUniformARB<4>(planeDistancePlaneLocation,1,planeEq);
		
		/* Set the texture unit variable: */
		glUniformARB(planeDistanceMapLocation,textureUnit);
		}
	}

void PointBasedLightingShader::disable(void)
	{
	/* Disable the shader: */
	glUseProgramObjectARB(0);
	}
