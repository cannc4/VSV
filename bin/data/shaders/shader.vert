#version 150

// these are for the programmable pipeline system
uniform mat4 modelViewProjectionMatrix;
in vec4 position;

void main()
{
	// finally set the pos to be that actual position rendered
	gl_Position = modelViewProjectionMatrix * position;
}