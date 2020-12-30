uniform mat4 gluOrtho;

void main()
{
  gl_Position = gluOrtho * gl_Vertex;
}