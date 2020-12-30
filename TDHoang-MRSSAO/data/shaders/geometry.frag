#version 140
#extension GL_ARB_texture_rectangle : require

in vec4 pos;
out vec4 Pos;
out vec4 Norm;

void main()
{    
  Pos = pos / pos.w;
    
  vec3 n = cross(dFdx(pos.xyz), dFdy(pos.xyz));    
  Norm = vec4(normalize(n), sign(abs(pos.z)));
}