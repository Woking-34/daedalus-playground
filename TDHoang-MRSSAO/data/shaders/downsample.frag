#version 140
#extension GL_ARB_texture_rectangle: enable

uniform sampler2DRect hiResNormTex;
uniform sampler2DRect hiResPosTex;

out vec4 Norm;
out vec4 Pos;

void main()
{
  vec2 xy2 = gl_FragCoord.xy * 2.0;

  vec4 pos[4];
  vec4 norm[4];

  pos[0] = texture2DRect(hiResPosTex, xy2 + vec2(-0.5, 0.5));
  pos[1] = texture2DRect(hiResPosTex, xy2 + vec2(0.5, 0.5));
  pos[2] = texture2DRect(hiResPosTex, xy2 + vec2(0.5, -0.5));
  pos[3] = texture2DRect(hiResPosTex, xy2 + vec2(-0.5, -0.5));

  norm[0] = texture2DRect(hiResNormTex, xy2 + vec2(-0.5, 0.5));
  norm[1] = texture2DRect(hiResNormTex, xy2 + vec2(0.5, 0.5));
  norm[2] = texture2DRect(hiResNormTex, xy2 + vec2(0.5, -0.5));
  norm[3] = texture2DRect(hiResNormTex, xy2 + vec2(-0.5, -0.5));

  float maxZ = max(max(pos[0].z, pos[1].z), max(pos[2].z, pos[3].z));
  float minZ = min(min(pos[0].z, pos[1].z), min(pos[2].z, pos[3].z));  

  int minPos, maxPos;

  for (int i = 0; i < 4; ++i)
  {
    if (pos[i].z == minZ)
      minPos = i;
    if (pos[i].z == maxZ)
      maxPos = i;
  }

  float d = distance(pos[minPos].xyz, pos[maxPos].xyz);  

  ivec2 median = ivec2(0, 0);
  int index = 0;

  for (int i = 0; i < 4 && index < 2; ++i)    
    if (i != minPos && i != maxPos)
      median[index++] = i;

  if (d < 1.0)
  {
    Pos = (pos[median.x] + pos[median.y]) / 2.0;
    Norm = (norm[median.x] + norm[median.y]) / 2.0;
  }
  else
  {
    Pos = pos[median.x];
    Norm = norm[median.x];
  }
}