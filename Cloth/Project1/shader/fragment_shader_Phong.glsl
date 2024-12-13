#version 330 core

in vec3 Normal;
in vec3 FragPos;

uniform vec3 lightPos;
uniform vec3 eye;
uniform vec3 color;
uniform int shadow;
uniform int backFace;

out vec4 FragColor;

void main()
{
    if(shadow==1)
        FragColor = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    else
    {
        vec3 norm;
        if(backFace==0)
            norm = normalize(Normal);
        else
            norm = -normalize(Normal);
        vec3 lightDir = normalize(lightPos-FragPos);
        float diffuse = max(dot(norm,lightDir),0.0);

        vec3 viewDir = normalize(eye-FragPos);
        vec3 reflectDir = reflect(-lightDir, norm);

        float specular = pow(max(0.0, dot(norm,reflectDir)),32);

        FragColor = vec4(diffuse * color + specular * color + 0.1f * vec3(0.7f, 0.9f, 1.0f), 1.0f);
    }
    

} 