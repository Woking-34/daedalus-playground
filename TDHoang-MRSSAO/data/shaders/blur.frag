#version 140
#extension GL_ARB_texture_rectangle: enable
#extension GL_EXT_gpu_shader4: enable

uniform sampler2DRect aoTex;
uniform sampler2DRect normTex;
uniform sampler2DRect posTex;

out vec3 AO;

void main()
{
  vec3 n = texture2DRect(normTex, gl_FragCoord.xy).xyz;
  vec3 p = texture2DRect(posTex, gl_FragCoord.xy).xyz;

  vec3 ss = vec3(0.0);

  float weight = 0.0;

  for (float i = -1.0; i <= 1.0; i += 1.0)
  {
    for (float j = -1.0; j <= 1.0; j += 1.0)
    {
      vec2 ij = vec2(i, j);
      vec3 t = texture2DRect(aoTex, gl_FragCoord.xy + ij).xyz;      
      vec3 norm = texture2DRect(normTex, gl_FragCoord.xy + ij).xyz;
      float depth = texture2DRect(posTex, gl_FragCoord.xy + ij).z;      
      
      float normWeight = (dot(norm, n) + 1.2) / 2.2;
      normWeight = pow(normWeight, 8.0);

      float depthWeight = 1.0 / (1.0 + abs(p.z - depth) * 0.2);
      depthWeight = pow(depthWeight, 16.0);

      float gaussianWeight = 1.0 / ((abs(i) + 1.0) * (abs(j) + 1.0));

      weight += normWeight * depthWeight * gaussianWeight;
      ss += t * normWeight * depthWeight * gaussianWeight;
    }
  }
  
  AO = vec3(ss / weight);
}