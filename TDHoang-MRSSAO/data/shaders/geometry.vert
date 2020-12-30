varying out vec4 pos;

void main()
{
  pos = gl_ModelViewMatrix * gl_Vertex;
  gl_Position = ftransform();
}