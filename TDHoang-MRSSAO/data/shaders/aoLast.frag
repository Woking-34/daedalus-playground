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
uniform sampler2DRect lastFrameAOTex;
uniform sampler2DRect lastFramePosTex;
uniform float dMax;
uniform float rMax;
uniform float r;
uniform float resolution;
uniform mat4 iMVMat;
uniform mat4 mVMat;
uniform mat4 projMat;
uniform float poissonDisk[32];

out vec4 AO;

vec3 n;
vec4 p;
float occlusion = 0.0;
float sampleCount = 0.0;

void ComputeOcclusion(vec2 xy)
{
  vec4 samplePos = texture2DRect(posTex, xy);
  vec4 sampleNorm =  texture2DRect(normTex, xy);
  float d = distance(p.xyz, samplePos.xyz);
  float t = min(1.0, (d * d) / (dMax * dMax));
  t = mix(1.0, 0.0, t);
  vec3 diff = normalize(samplePos.xyz - p.xyz);
  float cosTheta = max(dot(n, diff), 0.0);
  occlusion += t * cosTheta * sampleNorm.w;
  sampleCount += 1.0;
}

vec3 Upsample()
{
  vec2 loResCoord[4];
  loResCoord[0] = floor((gl_FragCoord.xy + vec2(-1.0,  1.0)) / 2.0) + vec2(0.5, 0.5);
  loResCoord[1] = floor((gl_FragCoord.xy + vec2( 1.0,  1.0)) / 2.0) + vec2(0.5, 0.5);
  loResCoord[2] = floor((gl_FragCoord.xy + vec2(-1.0, -1.0)) / 2.0) + vec2(0.5, 0.5);
  loResCoord[3] = floor((gl_FragCoord.xy + vec2( 1.0, -1.0)) / 2.0) + vec2(0.5, 0.5);
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

  for (int i = 0; i < 32; i += 2)
  {     
    ComputeOcclusion(gl_FragCoord.xy + vec2(poissonDisk[i], poissonDisk[i + 1]) * rangeMax);    
  }

  vec3 upsample = Upsample();
  float aoMax = max(upsample.x, occlusion / sampleCount);  
  float aoAverage = (upsample.y + occlusion) / (upsample.z + sampleCount);

  float currentFrameAO = (1.0 - aoMax) * (1.0 - aoAverage); 
  
  vec4 lastFrameEyePos = mVMat * (iMVMat * p);
  vec4 lastFrameScreenPos = projMat * lastFrameEyePos;
  vec2 lastFrameTexCoord = lastFrameScreenPos.xy / lastFrameScreenPos.w;
  lastFrameTexCoord = floor((lastFrameTexCoord + 1.0) * resolution / 2.0) + vec2(0.5, 0.5);
  float lastFrameZ = texture2DRect(lastFramePosTex, lastFrameTexCoord).z;
  float lastFrameWeight = 0.0;
  if (abs(1.0 - lastFrameZ / lastFrameEyePos.z) < 0.01)
    lastFrameWeight = 0.6;
  if (lastFrameTexCoord.x < 0.0 || lastFrameTexCoord.x > resolution ||
      lastFrameTexCoord.y < 0.0 || lastFrameTexCoord.y > resolution)
    lastFrameWeight = 0.0;
  float lastFrameAO = texture2DRect(lastFrameAOTex, lastFrameTexCoord).x;
  
  AO = vec4(lastFrameAO * lastFrameWeight + currentFrameAO * (1.0 - lastFrameWeight));
}