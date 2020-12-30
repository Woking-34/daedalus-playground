/* There is code duplication among the ao*.frag files. May be an include 
   mechanism for GLSL source file can resolve that... */

#version 140
#extension GL_ARB_texture_rectangle: enable
#extension GL_EXT_gpu_shader4: enable

uniform sampler2DRect normTex;
uniform sampler2DRect posTex;
uniform float dMax;
uniform float rMax;
uniform float r;

out vec4 AO;

vec3 n;
vec4 p;
float occlusion = 0.0;
float sampleCount = 0.0001;

void ComputeOcclusion(vec2 xy)
{
  vec4 samplePos = texture2DRect(posTex, xy);
  vec4 sampleNorm =  texture2DRect(normTex, xy);
  float d = distance(p.xyz, samplePos.xyz);  
  float t = min(1.0, (d * d) / (dMax * dMax));  
  t = 1.0 - t;  
  vec3 diff = normalize(samplePos.xyz - p.xyz);  
  float cosTheta = max(dot(n, diff), 0.0);  
  occlusion += t * cosTheta * sampleNorm.w;     
  sampleCount += 1.0;    
}

void main()
{
  n = texture2DRect(normTex, gl_FragCoord.xy).xyz;
  p = texture2DRect(posTex, gl_FragCoord.xy);

  float rangeMax = min(r / abs(p.z), rMax);    

  for (float x = 1.0; x <= rangeMax; x += 2.0)
  {
    for (float y = 1.0; y <= rangeMax; y += 2.0)
    {
      ComputeOcclusion(gl_FragCoord.xy + vec2(x, y));
      ComputeOcclusion(gl_FragCoord.xy + vec2(-x, y));
      ComputeOcclusion(gl_FragCoord.xy + vec2(-x, -y));
      ComputeOcclusion(gl_FragCoord.xy + vec2(x, -y));
    }
  }

  AO = vec4(occlusion / sampleCount, occlusion, sampleCount, 0.0);  
}