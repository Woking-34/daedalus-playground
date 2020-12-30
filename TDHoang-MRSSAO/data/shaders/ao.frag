/* There is code duplication among the ao*.frag files. May be an include 
   mechanism for GLSL source file can resolve that... */

#version 140
#extension GL_ARB_texture_rectangle: enable
#extension GL_EXT_gpu_shader4: enable

uniform sampler2DRect normTex;
uniform sampler2DRect posTex;
uniform sampler2DRect loResAOTex;
uniform sampler2DRect loResNormTex;
uniform sampler2DRect loResPosTex;
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

vec3 Upsample()
{
  vec2 loResCoord[4];
  loResCoord[0] = floor((gl_FragCoord.xy + vec2(-1.0, 1.0)) / 2.0) + vec2(0.5, 0.5);    
  loResCoord[1] = floor((gl_FragCoord.xy + vec2(1.0, 1.0)) / 2.0) + vec2(0.5, 0.5);
  loResCoord[2] = floor((gl_FragCoord.xy +vec2(-1.0, -1.0)) / 2.0) + vec2(0.5, 0.5);
  loResCoord[3] = floor((gl_FragCoord.xy + vec2(1.0, -1.0)) / 2.0) + vec2(0.5, 0.5);    
  vec3 loResAO[4];
  vec3 loResNorm[4];
  float loResDepth[4];
  for (int i = 0; i < 4; ++i)
  {
    loResNorm[i] = texture2DRect(loResNormTex, loResCoord[i]).xyz;
    loResDepth[i] = texture2DRect(loResPosTex, loResCoord[i]).z;
    loResAO[i] = texture2DRect(loResAOTex, loResCoord[i]).xyz;
  }      
  float normWeight[4];
  for (int i = 0; i < 4; ++i)
  {    
    normWeight[i] = (dot(loResNorm[i], n) + 1.1) / 2.1;
    normWeight[i] = pow(normWeight[i], 8.0);
  }    
  float depthWeight[4];
  for (int i = 0; i < 4; ++i)
  {
    depthWeight[i] = 1.0 / (1.0 + abs(p.z - loResDepth[i]) * 0.2);
    depthWeight[i] = pow(depthWeight[i], 16.0);
  }       
  float totalWeight = 0.0;
  vec3 combinedAO = vec3(0.0);
  for (int i = 0; i < 4; ++i)
  {
    float weight = normWeight[i] * depthWeight[i] * (9.0 / 16.0) /
      (abs((gl_FragCoord.x - loResCoord[i].x * 2.0) * (gl_FragCoord.y - loResCoord[i].y * 2.0)) * 4.0);    
    totalWeight += weight;
    combinedAO += loResAO[i] * weight;    
  }
  combinedAO /= totalWeight;  
  return combinedAO;
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

  vec3 upsample = Upsample();
  AO = vec4(max(upsample.x, occlusion / sampleCount), upsample.y + occlusion, upsample.z + sampleCount, 0.0);      
}