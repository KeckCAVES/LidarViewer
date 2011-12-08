/***********************************************************************
PointBasedLightingShader - Class to maintain a GLSL point-based lighting
shader that tracks the current OpenGL lighting state.
Copyright (c) 2008-2011 Oliver Kreylos

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

#include <string.h>
#include <stdio.h>
#include <iostream>
#include <Misc/ThrowStdErr.h>
#include <GL/gl.h>
#include <GL/GLExtensionManager.h>
#include <GL/Extensions/GLARBShaderObjects.h>
#include <GL/Extensions/GLARBVertexShader.h>
#include <GL/Extensions/GLARBFragmentShader.h>

/*****************************************
Methods of class PointBasedLightingShader:
*****************************************/

void PointBasedLightingShader::compileShader(void)
	{
	std::string vertexShaderFunctions;
	std::string vertexShaderMain;
	
	/* Create the main vertex shader starting boilerplate: */
	vertexShaderMain+=
		"\
		void main()\n\
			{\n\
			/* Compute the vertex position in eye coordinates: */\n\
			vec4 vertexEc=gl_ModelViewMatrix*gl_Vertex;\n\
			\n\
			/* Compute the normal vector in eye coordinates: */\n\
			vec3 normalEc=normalize(gl_NormalMatrix*gl_Normal);\n\
			\n\
			/* Let the normal vector always point towards the eye: */\n\
			if(dot(normalEc,vertexEc.xyz)>0.0)\n\
				normalEc=-normalEc;\n\
			\n";
	
	/* Get the material components: */
	if(usePointColors)
		{
		vertexShaderMain+=
			"\
			/* Get the material properties from the current color: */\n\
			vec4 ambient=gl_Color;\n\
			vec4 diffuse=gl_Color;\n\
			vec4 specular=gl_FrontMaterial.specular;\n\
			float shininess=gl_FrontMaterial.shininess;\n\
			\n";
		}
	else
		{
		vertexShaderMain+=
			"\
			/* Get the material properties from the material state: */\n\
			vec4 ambient=gl_FrontMaterial.ambient;\n\
			vec4 diffuse=gl_FrontMaterial.diffuse;\n\
			vec4 specular=gl_FrontMaterial.specular;\n\
			float shininess=gl_FrontMaterial.shininess;\n\
			\n";
		}
	
	/* Continue the main vertex shader: */
	vertexShaderMain+=
			"\
			/* Calculate global ambient light term: */\n\
			vec4 ambientDiffuseAccum=gl_LightModel.ambient*ambient;\n\
			vec4 specularAccum=vec4(0.0,0.0,0.0,0.0);\n\
			\n\
			/* Accumulate all enabled light sources: */\n";
	
	/* Create light application functions for all enabled light sources: */
	for(int lightIndex=0;lightIndex<lightTracker.getMaxNumLights();++lightIndex)
		if(lightTracker.getLightState(lightIndex).isEnabled())
			{
			/* Create the light accumulation function: */
			vertexShaderFunctions+=lightTracker.createAccumulateLightFunction(lightIndex);
			
			/* Call the light application function from the shader's main function: */
			char call[256];
			snprintf(call,sizeof(call),"\t\t\taccumulateLight%d(vertexEc,normalEc,ambient,diffuse,specular,shininess,ambientDiffuseAccum,specularAccum);\n",lightIndex);
			vertexShaderMain+=call;
			}
	
	/* Finish the main function: */
	vertexShaderMain+=
		"\
			\n\
			/* Compute final vertex color: */\n\
			gl_FrontColor=ambientDiffuseAccum+specularAccum;\n\
			\n\
			/* Use standard vertex position: */\n\
			gl_Position=ftransform();\n\
			}\n";
	
	/* Compile the vertex shader: */
	std::string vertexShaderSource=vertexShaderFunctions+vertexShaderMain;
	glCompileShaderFromString(vertexShader,vertexShaderSource.c_str());
	
	/* Compile the standard fragment shader: */
	const char* fragmentShaderSource=
		"\
		void main()\n\
			{\n\
			gl_FragColor=gl_Color;\n\
			}\n";
	glCompileShaderFromString(fragmentShader,fragmentShaderSource);
	
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
	}

PointBasedLightingShader::PointBasedLightingShader(const GLLightTracker& sLightTracker)
	:lightTracker(sLightTracker),lightStateVersion(0),
	 usePointColors(false),
	 vertexShader(0),fragmentShader(0),
	 programObject(0)
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
	}

PointBasedLightingShader::~PointBasedLightingShader(void)
	{
	glDeleteObjectARB(programObject);
	glDeleteObjectARB(vertexShader);
	glDeleteObjectARB(fragmentShader);
	}

void PointBasedLightingShader::setUsePointColors(bool newUsePointColors)
	{
	if(usePointColors!=newUsePointColors)
		{
		/* Invalidate the shader program: */
		lightStateVersion=0;
		}
	usePointColors=newUsePointColors;
	}

void PointBasedLightingShader::enable(void)
	{
	try
		{
		/* Re-compile the shader if it is out of line with the light tracker's lighting state: */
		if(lightStateVersion!=lightTracker.getVersion())
			{
			lightStateVersion=lightTracker.getVersion();
			compileShader();
			}
		
		/* Enable the shader: */
		glUseProgramObjectARB(programObject);
		}
	catch(std::runtime_error err)
		{
		std::cerr<<"Disabling lighting shader due to exception "<<err.what()<<std::endl;
		}
	}

void PointBasedLightingShader::disable(void)
	{
	/* Disable the shader: */
	glUseProgramObjectARB(0);
	}
