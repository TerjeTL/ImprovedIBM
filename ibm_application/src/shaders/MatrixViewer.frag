#version 460 core

out vec4 FragColor;

in vec2 TexCoord;

uniform uint size;
uniform vec3 color;

uniform sampler2D texture1;
uniform isampler2D texture2;

const float N = 20.0; // grid ratio
float gridTextureGradBox( in vec2 p, in vec2 ddx, in vec2 ddy )
{
	// filter kernel
    vec2 w = max(abs(ddx), abs(ddy)) + 0.01;

	// analytic (box) filtering
    vec2 a = p + 0.5*w;                        
    vec2 b = p - 0.5*w;           
    vec2 i = (floor(a)+min(fract(a)*N,1.0)-
              floor(b)-min(fract(b)*N,1.0))/(N*w);
    //pattern
    return (1.0-i.x)*(1.0-i.y);
}

const float smoothstep_factor = 0.1; // 10%
float GridNodes(in vec2 pt, in float n, in float radius)
{
    //scale up by n (representing number of nodes)
    pt = pt*n;
    
    // get the fractional part
    pt = fract(pt);

    const vec2 center = vec2(0.5); 
    vec2 len = pt-center;
    return 1.-smoothstep(radius - (radius*smoothstep_factor),
                         radius + (radius*smoothstep_factor),
                         dot(len,len)*4.0);
}

void main()
{
    vec2 coord = TexCoord * size + vec2(0.025, 0.025) - vec2(0.5, 0.5);
    vec3 material = vec3(1.0)*gridTextureGradBox( coord, dFdx(coord), dFdy(coord) );
    
    float test = texture(texture1, TexCoord).r;
    int node_flag = texture(texture2, TexCoord).r;
    float grid_nodes_int = GridNodes(TexCoord, size, 0.05);
    if (grid_nodes_int > 0.0)
    {
        if (node_flag == 0)
        {
            material.z = max(grid_nodes_int, material.x);
            material.y = min(1.-grid_nodes_int, material.y);
            material.x = min(1.-grid_nodes_int, material.z);
        }
        else if (node_flag == 2)
        {
            material.y = max(grid_nodes_int, material.x);
            material.x = min(1.-grid_nodes_int, material.y);
            material.z = min(1.-grid_nodes_int, material.z);
        }
        else
        {
            material.x = max(grid_nodes_int, material.x);
            material.y = min(1.-grid_nodes_int, material.y);
            material.z = min(1.-grid_nodes_int, material.z);
        }
    }

    // gamma correction	
	vec3 color = material;
	color = pow( color, vec3(0.85) );

    FragColor = vec4(color, 1.0);
    //FragColor = vec4(float(node_flag), 0.0, 0.0, 1.0); // view the node flag texture
};