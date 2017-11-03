#version 410 core

// Attributes passed on from the vertex shader
smooth in vec4 o_VertPosition;
smooth in vec3 o_VertNormal;

// Structure for holding light parameters
struct LightInfo {
    vec4 Position; // Light position in eye coords.
    vec3 La; // Ambient light intensity
    vec3 Ld; // Diffuse light intensity
    vec3 Ls; // Specular light intensity
};

// We'll have a single light in the scene
uniform LightInfo u_Light;
uniform LightInfo u_BackLight;
uniform vec4 u_Color;

// The material properties of our object
struct MaterialInfo {
    vec3 Ka; // Ambient reflectivity
    vec3 Kd; // Diffuse reflectivity
    vec3 Ks; // Specular reflectivity
    float Shininess; // Specular shininess factor
};
// The object has a material
uniform MaterialInfo u_Material;

// This is no longer a built-in variable
out vec4 o_FragColor;

void main() {
    // Calculate the normal
    vec3 n = normalize( o_VertNormal );

    // Calculate the light vector
    vec3 s = normalize( vec3(u_Light.Position) - vec3(o_VertPosition) );
    vec3 s2 = normalize( vec3(u_BackLight.Position) - vec3(o_VertPosition) );

    // Calculate the vertex position
    vec3 v = normalize(vec3(-o_VertPosition));

    // Reflect the light about the surface normal
    vec3 r = reflect( -s, n );
    vec3 r2 = reflect( -s2, n );

    // Compute the light from the ambient, diffuse and specular components
    vec3 lightColor = (
            u_Light.La * u_Material.Ka +
            u_Light.Ld * vec3(u_Color.rgb)/*u_Material.Kd*/ * max( dot(s, n), 0.0 ) +
            u_Light.Ls * u_Material.Ks * pow( max( dot(r,v), 0.0 ), u_Material.Shininess ));

    lightColor += (
            u_BackLight.La * u_Material.Ka +
            u_BackLight.Ld * vec3(u_Color.rgb)/*u_Material.Kd*/ * max( dot(s2, n), 0.0 ) +
            u_BackLight.Ls * u_Material.Ks * pow( max( dot(r2,v), 0.0 ), u_Material.Shininess ));

    // Set the output color of our current pixel
    o_FragColor = vec4(lightColor, 1.0);
}
