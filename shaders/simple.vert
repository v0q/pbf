#version 410 core

uniform mat4 u_MV;
uniform mat4 u_Projection;

layout(location = 0) in vec3 a_VertPosition;
layout(location = 2) in vec3 a_VertNormal;
layout(location = 1) in vec2 a_UV;

smooth out vec4 o_VertPosition;
smooth out vec3 o_VertNormal;

void main()
{
    // Transform the vertex normal by the inverse transpose modelview matrix
    o_VertNormal = normalize(a_VertNormal);

    // Compute the unprojected vertex position
    o_VertPosition = u_MV * vec4(a_VertPosition, 1.0f);

    gl_Position = u_Projection * o_VertPosition;
}
