//!HOOK POSTKERNEL
//!BIND HOOKED
//!BIND PREKERNEL

#define GetH(x,y)   HOOKED_texOff(vec2(x,y)).rgb
#define GetL(x,y)   PREKERNEL_tex(base + PREKERNEL_pt*vec2(x,y)).rgb
#define Luma(rgb)   ( dot(vec3(0.2126, 0.7152, 0.0722), rgb) )

vec4 hook() {
    vec2 fcoord = fract(PREKERNEL_pos * PREKERNEL_size - vec2(0.5));
    vec2 base = PREKERNEL_pos - fcoord * PREKERNEL_pt;

    float meanH = Luma((GetH(0, 0) + (GetH(-1, 0) + GetH(0, 1) + GetH(1, 0) + GetH(0, -1)) / 4.0) / 2.0);
    float meanL = Luma((GetL(0, 0) + (GetL(-1, 0) + GetL(0, 1) + GetL(1, 0) + GetL(0, -1)) / 4.0) / 2.0);

    // Get center four pixels
    vec3 c1 = GetL(0, 0);
    vec3 c2 = GetL(0, 1);
    vec3 c3 = GetL(1, 0);
    vec3 c4 = GetL(1, 1);

    vec3 lo = min(min(c1, c2), min(c3, c4));
    vec3 hi = max(max(c1, c2), max(c3, c4));

    vec4 color = HOOKED_texOff(0);

    float d = meanL/meanH;
    color.rgb = mix(clamp(color.rgb, lo, hi), color.rgb, mix(max(2.0-d, 0.0), d, step(d, 1.0)));

    return color;
}
